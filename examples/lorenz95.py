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

from endas.cov import DiagonalCovariance

sys.path.extend([os.path.abspath('..')])

import math

import numpy as np
import matplotlib.pyplot as plt

from endas import cov, obs
from endas import ensemble
from endas import algorithms as alg
from endas.localization import DomainLocalization, GenericStateSpace1d, taper

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
obs_skip = 12

# Smoother lag as the number of time steps. Lag 0 disables smoothing (i.e. only the filtering solution is obtained)
lag = 0

# Model is the Lorenz-95 dynamic system
model = lorenz.Lorenz95(n)

# Set up localization implementation. For this synthetic case we will use the generic 1d partitioning of the state
# space and the Gaspari-Cohn function with correlation length 2 for distance-based observation covariance tapering.
ls = None #DomainLocalization(ssp=GenericStateSpace1d(n), taper_fn=taper.GaspariCohn(2))

# Filters/smoothers we are going to run
kf_to_run = [

    ('KF', alg.KalmanFilter(model=model.__call__, model_tl=model.dot, model_adj=model.adjdot, lag=lag)),
    ('EnKF-noloc', alg.EnsembleKalmanFilter(variant=alg.EnKF(), ensemble_size=N, lag=lag))
#    ('EnKF-loc', alg.EnsembleKalmanFilter(variant=alg.EnKF(), ensemble_size=N, lag=lag, loc_strategy=ls))
]

# Observation operator: We will observe the last 3 variables in every 5 state vector variables, i.e. 24 out of 40
# in this case
k = 24
Hmat = np.zeros((k, n))
hi = np.arange(k, dtype=np.int32)
Hmat[hi, 5*(hi // 3) + (hi % 3) + 2] = 1
H = obs.MatrixObservationOp(Hmat)
observedStates = 5*(hi // 3) + (hi % 3) + 2


# Initial error covariance
sig_clim = 3.6414723
P0 = cov.DiagonalCovariance(np.ones(n) * (0.1*sig_clim)**2)

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
n_steps = 2000

x0 = np.ones(n) * 8.0 # Initial state
x0[20] = 8.004        # Perturb 20-th coordinate

xt, yobs = make_data(model_fn=model.__call__, x0=x0, dt=dt, H=H, Q=Q, R=R, nsteps=n_steps, nspin=0)


# ----------------------------------------------------------------------------------------------------------------------
# Ready to run...
# ----------------------------------------------------------------------------------------------------------------------

xall = np.zeros((len(kf_to_run), n_steps, n))
rmse = np.zeros((len(kf_to_run), n_steps))

# Called when a smoother solution is ready. Here we will just store the result and RMSE to be plotted later
def on_result(x, A, t, args):
    kfi, is_ensemble = args
    xall[kfi, t, :] = x
    rmse[kfi, t] = math.sqrt(np.mean((x - xt[:, t])**2))


# Generate the initial system state and ensemble that all algorithms will start from. The initial guess is
# purposefully not very good
x0 = np.ones(n)
A0 = ensemble.generate(N, x0, P0)

for kfi, (name, kf) in enumerate(kf_to_run):

    # Avoid non-deterministic output for testing
    np.random.seed(4321)

    is_ensemble = isinstance(kf, alg.EnsembleKalmanFilter)

    if is_ensemble:
        A = A0
    else:
        x = x0
        P = P0

    print("Running {}...".format(name))

    # Before the time stepping loop starts, let the smother know the initial system state. This way we will get
    # a smoother solution for it as well.
    if is_ensemble: kf.smoother_begin(A, 0)
    else: kf.smoother_begin(x, P, 0)

    # THIS IS THE MAIN "SIMULATION" TIME-STEPPING LOOP
    for t in range(1, n_steps):

        # Integrate the system state forward using the model. In this synthetic example we will use the kf.forecast()
        # method to simplify the code

        if is_ensemble: kf.forecast(model.__call__, A, Q, dt)
        else: kf.forecast(x, P, Q, dt)

        # Assimilate observations. Please note that we need to call begin_analysis() and end_analysis() even if
        # we do not have any observations. If we didn't, we would not get the smoother solution for these time steps.

        if is_ensemble: kf.begin_analysis(A, t)
        else: kf.begin_analysis(x, P, t)

        if t % obs_skip != 0:
            kf.assimilate(z=yobs[:, t], H=H, R=R)

        if is_ensemble: A = kf.end_analysis(on_smoother_result=on_result, result_args=(kfi, is_ensemble))
        else: x, P = kf.end_analysis(on_smoother_result=on_result, result_args=(kfi, is_ensemble))

    # Time-stepping loop completed. Call kf.finish() to get any pending smoother solutions. Please note that this is
    # not needed if only the filtering solution is of interest
    if is_ensemble: kf.finish(on_smoother_result=on_result, result_args=(kfi, is_ensemble))
    else: kf.finish(on_smoother_result=on_result, result_args=(kfi, is_ensemble))


# ----------------------------------------------------------------------------------------------------------------------
# Data assimilation completed, print RMSEs and plots
# ----------------------------------------------------------------------------------------------------------------------

print("Done")

#RMSEskipStart = 1000
#for i in range(len(filtersToRun)):
#    logging.info("RMSE %s : %.4f", filtersToRun[i].name, np.mean(rmse[RMSEskipStart:,i]))

print("Plotting...")
plt.ion()
line_styles = ['k-', 'b-', 'g-', 'm-', 'c-', 'y-']

# Plot the trajectory of the 29-th state variable for the entire DA run
fig = plt.figure()
plt.plot(range(n_steps), xt[29, 0:n_steps], 'k--', linewidth=1.3, alpha=0.5, label="truth")

for kfi, (name, kf) in enumerate(kf_to_run):
    plt.plot(range(n_steps), xall[kfi, 0:n_steps, 29],
             line_styles[kfi], linewidth=1.0, label=name)

obs_t = np.arange(1, n_steps, obs_skip)
plt.plot(obs_t, yobs[17, obs_t], 'kx', markersize=3)

plt.legend()
plt.grid(True)
plt.ylabel("$x_{29}$")
plt.xlabel("t")


# RMSE plot
fig = plt.figure()

for kfi, (name, kf) in enumerate(kf_to_run):
    plt.plot(range(n_steps), rmse[kfi, 0:n_steps], line_styles[kfi], linewidth=0.8, label=name)

plt.legend()
plt.grid(True)
plt.ylabel("RMSE")
plt.xlabel("t")

plt.show()































