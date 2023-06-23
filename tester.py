#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:41:45 2023

@author: Aaron
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi, arcsin, sqrt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import time


a = 0.246 # Graphene lattice constant in nm

""" Lorentzian as used in spectroscopy """
def lorentzian(x, b, c, x0, lw):
   f = b + c / pi * lw / ((x-x0)**2 + lw**2)
   return(f)

def jacobian(x, b, c, x0, lw):
    # Define partial derivatives
    Ldiffb = np.ones(x.shape)
    Ldiffc = lw / pi / ((x-x0)**2 + lw**2)
    Ldiffx0 = c / pi * 2 * lw * (x-x0) / ((x-x0)**2 + lw**2)**2
    Ldifflw = c / pi * (1 / ((x-x0)**2 + lw**2) -  2 * lw**2 / ((x-x0)**2 + lw**2)**2 )
    jac = np.stack((Ldiffb, Ldiffc, Ldiffx0, Ldifflw)).T
    return(jac)


x, y = 120, 56

p0 = [0, 250, 0, 2]

fitrange = (240, 290)

# This expression searches the xlist for the wave number values 
# closest to the chosen fitrange and retrieves their indexes
i_start = np.abs(xdata - fitrange[0]).argmin()
i_stop = np.abs(xdata - fitrange[1]).argmin()

# Define boundaries to speed up the fitting routine
lbounds = np.array((-np.inf, 0,      0,  0))
ubounds = np.array((np.inf,  np.inf, np.inf, np.inf))

# Retrieve the spectrum of the current scan data point
ydata = data[y,x]

# Slice ydata to keep only the relevant peak
# Slice xdata to keep only the relevant peak
xdata_fit = xdata[i_start:i_stop]
ydata_fit = ydata[i_start:i_stop]

# Determine highest data point, use position as x0 starting value
p0[2] = xdata_fit[np.argmax(ydata_fit)]

popt, pcov = curve_fit(lorentzian, xdata_fit, ydata_fit, p0, bounds=(lbounds, ubounds), jac=jacobian)
    
plt.scatter(xdata_fit, ydata_fit, zorder=1)
plt.plot(xdata_fit, lorentzian(xdata_fit, *p0))
plt.plot(xdata_fit, lorentzian(xdata_fit, *popt))
plt.xlabel(r"Raman shift (cm$^{-1})$")
plt.ylabel("Intensity (a.u.)")
plt.show()
print(popt)