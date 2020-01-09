#!/usr/bin/env python3
# -*- coding: utf-8; -*-
# Simple test of DFS

import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from endas import ensemble
from endas import algorithms as algs
from endas.obs import MatrixObservationOp
from endas.cov import DiagonalCovariance

# For deterministic output
np.random.seed(1234)

# EnKF ensemble size
N = 20

# Number of data assimilation steps. The model has a period of 12 time steps.
nsteps = 12 * 10

# Number of time steps for the lagged smoother
smoother_lag = 10

# Observation error standard deviation
obs_sigma = 0.4

# Initial system state
x0 = np.array([0, 0, 1], dtype=np.double)

# Model error covariance matrix
Q = DiagonalCovariance(diag=np.array([0.08, 0.01, 0.01])**2)

# Initial (background) error covariance matrix
P0 = DiagonalCovariance(diag=np.array([1., 1., 1.])**2)

# Observation operator. This is a 1x3 matrix representing a single observation
# of both the first and second state variable. The third state variable is unobserved.
H = MatrixObservationOp(np.array([[1., 1., 0.]]))

# Observation error covariance matrix
R = DiagonalCovariance(np.array([obs_sigma ** 2]))



# Simple dynamic model for propagating the system state forward in time.
# The model is a rotation of the state, around the first dimension, by pi/6. The model
# is therefore periodic with 12 steps needed to complete one cycle. The tangent linear
# and adjoint operations required by KalmanFilter are coded via the `dot()` and
# `adjdot()` methods.
class SimpleModel:
    def __init__(self):
        self.M = np.array([[1.0, 0.0, 0.0],
                           [0.0, np.cos(np.pi / 6.0), np.sin(np.pi / 6.0)],
                           [0.0, -np.sin(np.pi / 6.0), np.cos(np.pi / 6.0)]])

    # Call needs to update `x` in-place and can return any data needed for the tangent
    # linear and adjoint. With simple linear model neither the tangent linear nor adjoint
    # depends on `x`
    def __call__(self, x, dt):
        x[:] = self.M.dot(x)
        return None

    # Implements Mx
    def dot(self, _, x):
        return self.M.dot(x)

    # Implements xM'
    def adjdot(self, _, x):
        return x.dot(self.M.T)

model = SimpleModel()


# Before we start we generate some synthetic data that will serve as the "true" state
# and draw observations by applying the observation operator H and adding noise with
# covariance R.

x = x0
xtrue = np.zeros((nsteps, 3))
xtrue[0,:] = x
y = np.zeros((nsteps, 1))  # Array of observations

for i in range(1, nsteps):
    model(x, 1)
    x+= Q.random_multivariate_normal()
    xtrue[i,:] = x
    y[i,:] = H.dot(x) + R.random_multivariate_normal()



def onresult(x, A, t, result_args):
    est, err, is_enkf = result_args
    # print(k,t,x, np.std(E,axis=1))
    est[t, :] = x
    err[t, :] = np.std(A, axis=1) if is_enkf else np.diagonal(A)**0.5


# Collect results here
kf_est = np.zeros((nsteps, 3))
kf_std = np.zeros((nsteps, 3))
enkf_est = np.zeros((nsteps, 3))
enkf_std = np.zeros((nsteps, 3))

enkf = algs.EnsembleKalmanFilter(variant=algs.ESTKF(), ensemble_size=N, lag=smoother_lag)

kf = algs.KalmanFilter(model=model.__call__, model_tl=model.dot, model_adj=model.adjdot,
                       lag=smoother_lag)

# Run DA with EnKF
print("Running EnKF...")

result_args = (enkf_est, enkf_std, True)
A = ensemble.generate(N, x0, P0)
enkf.smoother_begin(A, 0)
x = x0
for t in range(1, nsteps):
    A = enkf.forecast(model.__call__, A, Q, 1)
    enkf.begin_analysis(A, t)
    enkf.assimilate(y[t], None, H, R)
    A = enkf.end_analysis(on_smoother_result=onresult, result_args=result_args)

enkf.smoother_finish(on_smoother_result=onresult, result_args=result_args)


# Run DA with KF
print("Running KF...")
result_args = (kf_est, kf_std, False)
kf.smoother_begin(x0, P0, 0)
x = x0
P = P0
for t in range(1, nsteps):
    x, P = kf.forecast(xb=x, Pb=P, Q=Q, dt=1)
    kf.begin_analysis(x, P, t)
    kf.assimilate(y[t], H, R)
    x, P = kf.end_analysis(on_smoother_result=onresult, result_args=result_args)

kf.smoother_finish(on_smoother_result=onresult, result_args=result_args)


# Plot results
print("Plotting...")

plt.figure(1, figsize=(10, 5))
plt.clf()
fig = plt.gcf()
# fig.set_size_inches(8, 5)
#res = np.flipud(np.vstack(resu))
#resstd = np.flipud(np.vstack(ress))
yhat = enkf_est[:,0] + enkf_est[:, 1]
ystd = np.sqrt(enkf_std[:, 0]**2 + enkf_std[:, 1]**2)

#t = np.arange(res.shape[0])
plt.plot(range(nsteps), yhat, "g-")
plt.fill_between(range(nsteps), yhat - 1.96 * ystd, yhat + 1.96 * ystd,
                 facecolor='green', alpha=0.3, label='enks')

# exact results
yhat = kf_est[:,0] + kf_est[:, 1]
ystd = np.sqrt(kf_std[:, 0]**2 + kf_std[:, 1]**2)

#t = np.arange(res.shape[0])
plt.plot(range(nsteps), yhat, "r-")
plt.fill_between(range(nsteps), yhat - 1.96 * ystd, yhat + 1.96 * ystd,
                 facecolor='red', alpha=0.3, label='ks')

# Truth
xt = xtrue[:,0] + xtrue[:, 1]
plt.plot(range(nsteps), xt, "k-")

#Observations
plt.plot(range(1, nsteps), y[1:], 'ks', markersize=3, label='Observation')


plt.grid(True)
plt.legend()
plt.title('Nens = {}'.format(N))
plt.show()

