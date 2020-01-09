"""
Data assimilation example on the Lorenz 95 dynamical system.

To run from console, use:

    cd examples
    python3 -i lorenz95.py

The ``-i`` switch is important if you are using regular Python- it enables the interactive mode on the interpreter and
allows the figures that are generated at the end to stay visible after the script ends. Anaconda apparently does not
have this issue and the -i switch isn't needed.
"""


# Make sure the endas package can be found when running from the examples folder without having to install it first
import sys, os.path

import math

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from endas import cov, obs
from endas import ensemble
from endas import algorithms as alg
from endas.cov import DiagonalCovariance
from endas.localization import DomainLocalization, GenericStateSpace1d, taper

from scipy import linalg

from models import lorenz
from utils import make_data


# Avoid non-deterministic output for testing
np.random.seed(1234)

# ----------------------------------------------------------------------------------------------------------------------
# Problem setup
# ----------------------------------------------------------------------------------------------------------------------

# State space dimension
n = 40

# Ensemble size
N = 20

# Number of observations we will skip before we assimilate one. 1 step = 30min, 12 = 6h
obs_skip = 6

# Smoother lag as the number of time steps. Lag 0 disables smoothing (i.e. only the filtering solution is obtained)
lag = 50

# Model is the Lorenz-95 dynamic system
model = lorenz.Lorenz95(n)

# Set up localization implementation. For this synthetic case we will use the generic 1d partitioning of the state
# space and the Gaspari-Cohn function with correlation length 2 for distance-based observation covariance tapering.
ls = DomainLocalization(ssp=GenericStateSpace1d(n), taper_fn=taper.GaspariCohn(2))

# Filters/smoothers we are going to run

kf = alg.KalmanFilter(model=model.__call__, model_tl=model.dot, model_adj=model.adjdot, lag=lag, forgetting_factor=1.0)
enkfs = [
    #('EnKF-noloc', alg.EnsembleKalmanFilter(variant=alg.EnKF(), ensemble_size=N, lag=lag, forgetting_factor=1.0)),
    ('EnKF-loc', alg.EnsembleKalmanFilter(variant=alg.EnKF(), ensemble_size=N, lag=lag, forgetting_factor=1.0, loc_strategy=ls)),
    #('ETKF-noloc', alg.EnsembleKalmanFilter(variant=alg.ESTKF(), ensemble_size=N, lag=lag, forgetting_factor=1.0)),
    ('ETKF-loc', alg.EnsembleKalmanFilter(variant=alg.ESTKF(), ensemble_size=N, lag=lag, loc_strategy=ls))
]

# Observation operator: We will observe the last 3 variables in every 5 state vector variables, i.e. 24 out of 40
# in this case
k = 24
Hmat = np.zeros((k, n))
hi = np.arange(k, dtype=np.int32)
Hmat[hi, 5*(hi // 3) + (hi % 3) + 2] = 1
H = obs.MatrixObservationOp(Hmat)
observedStates = 5*(hi // 3) + (hi % 3) + 2

k = 40
H = obs.MatrixObservationOp(np.eye(n))
observedStates = np.arange(0, k)



# Initial error covariance
sig_clim = 3.6414723
P0 = cov.DiagonalCovariance(np.ones(n) * (0.5*sig_clim)**2)

# Model error covariance
Q = cov.DiagonalCovariance(np.ones(n) * (0.05*sig_clim)**2)

# Observation error covariance
R = cov.DiagonalCovariance(np.ones(k) * (0.15*sig_clim)**2)


# ----------------------------------------------------------------------------------------------------------------------
# Generate synthetic data for the twin experiment
# ----------------------------------------------------------------------------------------------------------------------

print("Generating true state and observations...")

dt = 0.025 / 6.0    # Model integration time step (equivalent to 30min)
n_steps_in_day = 48
n_steps = 1000

x0 = np.ones(n) * 8.0 # Initial state
x0[20] = 8.004        # Perturb 20-th coordinate

xt, yobs = make_data(model_fn=model.__call__, x0=x0, dt=dt, H=H, Q=Q, R=R, nsteps=n_steps, nspin=0)


# ----------------------------------------------------------------------------------------------------------------------
# Ready to run...
# ----------------------------------------------------------------------------------------------------------------------

xall = np.zeros((len(enkfs)+1, n_steps, n))
rmse = np.zeros((len(enkfs)+1, n_steps))

# Called when a smoother solution is ready. Here we will just store the result and RMSE to be plotted later
def on_result(x, A, t, args):
    kfi, is_ensemble = args
    xall[kfi, t, :] = x
    rmse[kfi, t] = math.sqrt(np.mean((x - xt[:, t])**2))


# Generate the initial system state and ensemble that all algorithms will start from. The initial guess is
# purposefully not very good
x0 = np.ones(n)
A0 = ensemble.generate(N, x0, P0)

# Run all Ensemble Kalman Filters/Smoothers

for kfi, (name, enkf) in enumerate(enkfs):
    print("Running {}...".format(name))

    # Avoid non-deterministic output for testing
    np.random.seed(1234)

    A = np.copy(A0)
    # Before the time stepping loop starts, let the smother know the initial system state. This way we will get
    # a smoother solution for it as well.
    enkf.smoother_begin(A, 0)

    # THIS IS THE MAIN "SIMULATION" TIME-STEPPING LOOP
    for t in range(1, n_steps):

        # Integrate the system state forward using the model. In this synthetic example we will use the forecast()
        # method to simplify the code
        A = enkf.forecast(model.__call__, A, Q, dt)

        # Assimilate observations. Please note that we need to call begin_analysis() and end_analysis() even if
        # we do not have any observations. If we didn't, we would not get the smoother solution for these time steps.
        enkf.begin_analysis(A, t)

        if t % obs_skip == 0: enkf.assimilate(z=yobs[:, t], H=H, R=R, z_coords=observedStates)

        A = enkf.end_analysis(on_smoother_result=on_result, result_args=(kfi, True))

    # Time-stepping loop completed. Call finish() to get any pending smoother solutions. Please note that this is
    # not needed if only the filtering solution is of interest
    enkf.smoother_finish(on_smoother_result=on_result, result_args=(kfi, True))


# Do an exact Kalman Filter/Smoother run for reference. The API is similar to the ensemble smoothing API, see the
# comments above for explanation

print("Running KF...")
np.random.seed(1234)

x = np.copy(x0)
P = P0.to_matrix(force_dense=True)

kf.smoother_begin(x, P, 0)

# THIS IS THE MAIN "SIMULATION" TIME-STEPPING LOOP
for t in range(1, n_steps):
    x, P = kf.forecast(x, P, Q, dt)

    kf.begin_analysis(x, P, t)
    if t % obs_skip == 0: kf.assimilate(z=yobs[:, t], H=H, R=R)
    x, P = kf.end_analysis(on_smoother_result=on_result, result_args=(-1, False))

kf.smoother_finish(on_smoother_result=on_result, result_args=(-1, False))




print("Running REFERENCE...")

np.random.seed(1234)

xall2 = np.zeros((n_steps, n))
rmse2 = np.zeros((n_steps, n))

x = np.copy(x0)
P = P0.to_matrix(force_dense=True)
QQ = Q.to_matrix(force_dense=True)

HH = H.to_matrix(force_dense=True)

for t in range(1, n_steps):

    # Forecast
    trj = model(x, dt)
    P = model.dot(trj, P)
    P = model.adjdot(trj, P)
    if Q is not None: P += QQ


    # Update
    if t % obs_skip == 0:
        F = HH.dot(P).dot(HH.T)

        if isinstance(R, np.ndarray):
            np.add(F, R, out=F)
        elif isinstance(R, cov.CovarianceOperator):
            R.add_to(F)

        z = yobs[:, t]

        # State update as xk + Cp*H'*F^-1*dz
        dz = z - HH.dot(x)
        x += P.dot(HH.T).dot(linalg.solve(F, dz, sym_pos=True, overwrite_b=True))
        #self._xa = self._xa.ravel()

        # Covariance estimate update as Cp - Cp*H'*F^-1*H*Cp
        P -= P.dot(HH.T).dot(linalg.solve(F, HH.dot(P), sym_pos=True, overwrite_a=True, overwrite_b=True))

    xall2[t, :] = x
    rmse2[t, :] = np.diagonal(P).ravel()


# ----------------------------------------------------------------------------------------------------------------------
# Data assimilation completed, print RMSEs and plots
# ----------------------------------------------------------------------------------------------------------------------

print("Done")

#RMSEskipStart = 1000
#for i in range(len(filtersToRun)):
#    logging.info("RMSE %s : %.4f", filtersToRun[i].name, np.mean(rmse[RMSEskipStart:,i]))

print("Plotting...")

line_styles = ['b-', 'g-', 'm-', 'c-', 'y-']

OBS2PLOT = 17
X2PLOT = observedStates[OBS2PLOT]


# Plot the trajectory of the 29-th state variable for the entire DA run
fig = plt.figure()
plt.plot(range(n_steps), xt[X2PLOT, 0:n_steps], 'k--', linewidth=1.3, alpha=0.5, label="truth")

for kfi, (name, enkf) in enumerate(enkfs):
    plt.plot(range(n_steps), xall[kfi, 0:n_steps, X2PLOT], line_styles[kfi], linewidth=1.0, label=name)
    #plt.plot(range(n_steps), rmse[kfi, 0:n_steps, X2PLOT], 'g-', linewidth=1.0, label=name + ' rmse')

plt.plot(range(n_steps), xall[-1, 0:n_steps, X2PLOT], 'k-', linewidth=1.0, label='KF')


plt.plot(range(n_steps), xall2[0:n_steps, X2PLOT], 'r--', linewidth=1.0, label='REFERENCE')
#plt.plot(range(n_steps), rmse2[0:n_steps, X2PLOT], 'r:', linewidth=1.0, label='REFERENCE rmse')

obs_t = np.arange(1, n_steps, obs_skip)
plt.plot(obs_t, yobs[OBS2PLOT, obs_t], 'kx', markersize=3)

plt.legend()
plt.grid(True)
plt.ylabel("x29")
plt.xlabel("t")
plt.show()

# RMSE plot
# fig = plt.figure()
#
# for kfi, (name, kf) in enumerate(kf_to_run):
#     plt.plot(range(n_steps), rmse[kfi, 0:n_steps], line_styles[kfi], linewidth=0.8, label=name)
#
# plt.legend()
# plt.grid(True)
# plt.ylabel("RMSE")
# plt.xlabel("t")
































