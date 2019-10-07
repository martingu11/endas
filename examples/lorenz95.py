"""
Data assimilation example on the Lorenz 95 dynamical system.

To run from console, use:

    cd examples
    python3 -i lorenz95.py

The ``-i`` switch is important if you are using regular Python- it enables the interactive
mode on the interpreter and allows the figures that are generated at the end to stay visible
after the script ends. Anaconda apparently does not have this issue and the -i switch isn't
needed.
"""



# Make sure dfs package can be found when running from the examples folder
# without having to install the package first
import sys, os.path
sys.path.extend([os.path.abspath('..')])


import math
import logging

import numpy as np
import matplotlib.pyplot as plt

from endas import cov
from endas.localization import taper

from examples.models import lorenz
from endas import arraycache

from utils import makedata

# Avoid non-deterministic output for testing
np.random.seed(1234)


# Problem setup

n = 40              # State dimension


obsBatchSize = 100


# Model error covariance
sigClim = 3.6414723
q =  (0.05*sigClim)**2
Q = cov.DiagonalCovarianceOp(np.ones(n) * q)

# Observation operator: We will observe the last 3 variables in every 5 state
# vector variables, i.e. 24 out of 40 in this case
k = 24
Hmat = np.zeros((k, n))
hi = np.arange(k, dtype=np.int32)
Hmat[hi, 5*(hi // 3) + (hi % 3) + 2] = 1
H = observation.MatrixObservationOp(Hmat)

observedStates = 5*(hi // 3) + (hi % 3) + 2


# Observation error covariance
r = (0.15*sigClim)**2 #(0.15*sigClim)**2
R = cov.DiagonalCovarianceOp(np.ones(k) * r, blockSize=obsBatchSize)


# Model is the Lorenz-95 dynamic system
model = lorenz.Lorenz95(n)


# Generate synthetic data for the test

logging.info("Generating true state and observations...")

dt = 0.025 / 6.0    # Model integration time step (equivalent to 30min)
nStepsInDay = 48
nsteps = 2000


x0 = np.ones(n) * 8.0 # Initial state
x0[20] = 8.004        # Perturb 20-th coordinate

xt, yobs = makedata(model, x0=x0, dt=dt, H=H, Q=Q, R=R, nsteps=nsteps, nspin=0)



# Ready to run...

m = 20        # Ensemble size
obsskip = 12  # Number of observations we will skip before we assimilate one. 1 step = 30min, 12 = 6h

# Filters we are going to test
filtersToRun = [
  filter.EKF("EKF"),
  filter.EnKF(m, "EnKF"),
  filter.EnKF(m, "EnKF+loc"),
  filter.EnKF(m, "EnKS+loc"),
  #enkf_naive.EnKF(m, "EnKFnaive+loc")
  filter.EnSRF(m, "EnSRF+loc"),
  filter.EnSRF(m, "EnSRS+loc")
  #filters.VEnKF("VEnKF")
]

smoothers = ["EnKS", "EnSRS"]


xall = np.zeros((nsteps, n, len(filtersToRun)))
rmse = np.zeros((nsteps, len(filtersToRun)))

x0 = np.ones(n)  # Bad initial state guess
#x0 = xt[:,0] + Q.random_multivariate_normal()


dataCache = arraycache.ArrayCache()


for i in range(len(filtersToRun)):

  filter = filtersToRun[i]

  if "+loc" in filter.name:

    loc = localization.CovarianceLocalization(n, 1, taper.GaspariCohn(2))
    loc.set_obscoords(observedStates)
    filter.localize(loc)


    # localization.GaspariCohnTaperFn1d(2), observedStates)
    #filter.setCovInflationFactor(1.08)


  np.random.seed(4321)

  kf = filtersToRun[i]
  logging.info("Running %s..." % filter.name)

  x = x0
  P = Q

  dataProvider = ( (dt, yobs[:,i] if i % obsskip == 0 else None, H, Q, R) for i in range(nsteps) )

  def onResult(k, t, xa, Ea):
    xall[k-1, :, i] = xa
    rmse[k-1,i] = math.sqrt(np.mean((xa - xt[:,k-1])**2))


  # Running smoother or filter?
  if any(x in kf.name for x in smoothers):
    kf.smoother(model, x, P, dataProvider, dataCache, onResult, deflation=0.9)
  else:
    run.runfilter(kf, model, x, P, dataProvider, onResult)



# Done, print RMSEs and plots

logging.info("Done")


RMSEskipStart = 1000

for i in range(len(filtersToRun)):
  logging.info("RMSE %s : %.4f", filtersToRun[i].name, np.mean(rmse[RMSEskipStart:,i]))



logging.info("Runtime summary:")
runtimestats.print_summary(logging.getLogger())



logging.info("Plotting...")

filterNames = [f.name for f in filtersToRun]
lineStyles = ['k-', 'b-', 'g-', 'm-', 'c-', 'y-']



plt.ion()

# Trajectory of the 29-th state variable
fig = plt.figure()
plt.plot(range(nsteps), xt[29,0:nsteps], 'k--', linewidth=1.3, alpha=0.5)

for i in range(len(filtersToRun)):
  plt.plot(range(nsteps), xall[0:nsteps, 29,i], lineStyles[i], linewidth=1.0)

obsT = np.arange(0, nsteps, obsskip)
plt.plot(obsT, yobs[17, obsT], 'kx', markersize=3)

plt.legend(['truth'] + filterNames + ['observation'])
plt.grid(True)
plt.ylabel("$x_{29}$")
plt.xlabel("t")



# RMSE plot
fig = plt.figure()

for i in range(len(filtersToRun)):
  plt.plot(range(nsteps), rmse[0:nsteps, i], lineStyles[i], linewidth=0.8)

plt.legend(filterNames)
plt.grid(True)
plt.ylabel("RMSE")
plt.xlabel("t")

plt.show()































