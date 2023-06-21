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

a = 0.246 # Graphene lattice constant in nm

""" Lorentzian as used in spectroscopy """
def lorentzian(x, b, c, x0, lw):
   f = b + c / pi * lw / ((x-x0)**2 + lw**2)
   return(f)


test = (np.array((1,2,3)))
test_ma = np.ma.masked_array

# =============================================================================
# fitrange = (240,290)
# plotrange = (240,290)
# 
# i_start = min(range(len(xdata)), key=lambda i: abs(xdata[i]-fitrange[0]))
# i_stop = min(range(len(xdata)), key=lambda i: abs(xdata[i]-fitrange[1]))
# 
# x, y = 5, 18
# ydata = data[y, x]
# 
# # Slice data to keep only the relevant peak
# xdata_fit = xdata[i_start:i_stop]
# ydata_fit = ydata[i_start:i_stop]
# 
# # Optional: Replace x0 initial value by position of max value in the list
# x0 = xdata_fit[np.argmax(ydata_fit)]
# p0 = [10, 200, x0, 2]
# 
# fig = plt.figure()
# ax1 = fig.add_subplot(211)
# ax1.scatter(xdata_fit, ydata_fit, zorder=1)
# ax1.plot(xdata_fit, lorentzian(xdata_fit, *p0))
# ax1.set_xlabel(r"Raman shift (cm$^{-1})$")
# ax1.set_ylabel("Intensity (a.u.)")
# plt.tight_layout()
# 
# popt, pcov = curve_fit(lorentzian, xdata_fit, ydata_fit, p0)#, bounds=(lbounds, ubounds))
# print(popt)
# 
# fig = plt.figure()
# ax1 = fig.add_subplot(211)
# ax1.scatter(xdata_fit, ydata_fit, zorder=1)
# ax1.plot(xdata_fit, lorentzian(xdata_fit, *p0))
# ax1.plot(xdata_fit, lorentzian(xdata_fit, *popt))
# ax1.set_xlabel(r"Raman shift (cm$^{-1})$")
# ax1.set_ylabel("Intensity (a.u.)")
# plt.tight_layout()
# =============================================================================

