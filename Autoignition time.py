###############################################################
#
# Autoignition of a methane air mixture at stoichiometry
# and atmospheric pressure, for different
# initial temperature
#
###############################################################

import sys
import numpy as np
from cantera import *
import cantera as ct
from pylab import *
import csv
#################################################################
# Mechanism used for the process

gas = Solution('gri30.cti')
# Initial temperature, Pressure and stoichiometry

# Specify the number of time steps and the time step size
nt = 100000
dt = 1.e-6  # s
Tmin = 0.65
Tmax = 0.85
npoints = 11
pmin= 0.8
pmax = 1
wmin=0
wmax=1.44
# Storage
# Temperature storage variables

Ti = np.zeros(npoints, 'd')
Ti2 = np.zeros(npoints, 'd')

# The initial storage variable become case dependent

tim = np.zeros(nt, 'd')
temp_cas = np.zeros(nt, 'd')
dtemp_cas = np.zeros(nt - 1, 'd')
# Additional storage variables are needed to differentiate each case

#Pressure storage variable

pi = np.zeros(npoints, 'd')

#stezenie
w = np.zeros(npoints,'d')


Autoignition_cas = np.zeros(npoints, 'd')
FinalTemp_cas = np.zeros(npoints, 'd')
mfrac_cas = np.zeros([npoints, gas.n_species], 'd')

# Loop over initial conditions
for l in range(npoints):
    w[l]= wmin + (wmax +wmin) * l / (npoints - 1)
    procent = (w[l]/(5.76+w[l]))*100
    for k in range(npoints):
        pi[k]= pmin +(pmax-pmin) * k / (npoints-1)
        one_atm= one_atm *pi[k]
        for j in range(npoints):
            Ti2[j] = Tmin + (Tmax - Tmin) * j / (npoints - 1)
            Ti[j] = 1000 / Ti2[j]
        # Set gas state, always at stoichiometry
            gas.TPX = Ti[j], one_atm, ('C2H6:1,O2:1,N2:3.76,H2O:'+str(w[l]))
        # Create the ideal batch reactor
            r = ct.IdealGasReactor(gas)
        # Now create a reactor network consisting of the single batch reactor
            sim = ct.ReactorNet([r])
    # Initial simulation time
            time = 0.0
    # Loop for nt time steps of dt seconds.
            for n in range(nt):
                time += dt
                sim.advance(time)
                tim[n] = time
                temp_cas[n] = r.T
            mfrac_cas[j][:] = r.thermo.Y
#################################################################
# Catch the autoignition timings
#################################################################
# Get autoignition timing
            Dtmax = [0, 0.0]
            for n in range(nt - 1):
                dtemp_cas[n] = (temp_cas[n + 1] - temp_cas[n]) / dt
                if (dtemp_cas[n] > Dtmax[1]):
                    Dtmax[0] = n
                    Dtmax[1] = dtemp_cas[n]
# Local print
            Autoignition     = ((tim[Dtmax[0]] + tim[Dtmax[0] + 1]) / 2.)
            print ('For ' + str(Ti[j]) +',p='+str(pi[k])+',st='+str(w[l])+ ', Autoignition time = (s) ' + str(Autoignition))
# Posterity
            Autoignition_cas[j] = Autoignition * 1000  # ms
            FinalTemp_cas[j] = temp_cas[nt - 1]
        plot(Ti2, Autoignition_cas, '^', color='orange')
        xlabel(r'Temp [1000/K]', fontsize=20)
        ylabel("Autoignition [ms]")
        title(r'Autoignition of $C_{2}H_{6}$ + Air mixture at $\Phi$ = 1, P ='+ str(pi[k]) +'bar and H2O-'+str(procent)+'of mixture', fontsize=22,horizontalalignment='center')
        axis([0.60, 0.90, 0.0, 100.0])
        grid()
        show()