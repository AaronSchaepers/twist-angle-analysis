#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Aaron
@date: 2023/06/29
@github: https://github.com/AaronSchaepers/twist-angle-analysis
"""

###############################################################################
# Provide default values for the fit initial values and threshold intervals 
###############################################################################

# Starting parameters for Lorentzian peak fits in the order:
# Offset, intensity, linewidth.
# The starting value for the position is determined dynamically
initvalues_TA = [0, 250, 2]
initvalues_G = [0, 1E4, 16]
initvalues_LO = [0, 5E3, 4]
initvalues_2D = [0, 2E4, 15]

# Threshold parameters for excluding implausible fit results
# NOTE: These are not boundary values for the fitting routine!
#       Instead, the best fit parameters of each spectrum are compared to these  
#       threshold intervals. If any parameter lies outside of its threshold
#       interval, the corresponding data point is excluded from the map.

thresh_TA_c = [20, 1E4]      # Intensity
thresh_TA_x0 = [250, 275]    # Position
thresh_TA_lw = [0.2, 9]      # Linewidth

thresh_G_c = [1E3, 1E15]     # Intensity
thresh_G_x0 = [1577, 1587]   # Position
thresh_G_lw = [4, 25]        # Linewidth

thresh_LO_c = [300, 1E5]     # Intensity
thresh_LO_x0 = [1610, 1630]  # Position
thresh_LO_lw = [0.4, 15]     # Linewidth

thresh_2D_c = [1E3, 1E7]     # Intensity
thresh_2D_x0 = [2650, 2700]  # Position
thresh_2D_lw = [8, 30]       # Linewidth

max_gradient = 1 # °/µm, upper bound for the twist angle gradient map. Larger
                 #       values will be excluded from the map.

# Peak-specific spectral range for fitting (in wavenumbers)
fitrange_TA = (240, 290)
fitrange_G = (1500, 1610)
fitrange_LO = (1610, 1800)
fitrange_2D = (2600, 2800)

# Peak-specific spectral range for plotting (in wavenumbers)
plotrange_TA = (220, 310)
plotrange_G = (1500, 1800)
plotrange_LO = (1500, 1800)
plotrange_2D = (2400, 2900)



###############################################################################
# Functions that return a dicitonary with the required data
###############################################################################


def TA():
    dict_TA = {}
    # Insert a place holder value for the position starting value, which is 
    # determined dynamically in the fitting routine
    dict_TA["initvalues"] = [initvalues_TA[0], initvalues_TA[1], 0, initvalues_TA[2]]
    dict_TA["peakname"] = "TA"
    dict_TA["fitrange"] = fitrange_TA
    dict_TA["plotrange"] = plotrange_TA
    dict_TA["params_thresh"] = (thresh_TA_c, thresh_TA_x0, thresh_TA_lw)
    dict_TA["max_gradient"] = max_gradient
    return(dict_TA)

def G():
    dict_G = {}
    dict_G["initvalues"] = [initvalues_G[0], initvalues_G[1], 0, initvalues_G[2]]
    dict_G["peakname"] = "G"
    dict_G["fitrange"] = fitrange_G
    dict_G["plotrange"] = plotrange_G
    dict_G["params_thresh"] = (thresh_G_c, thresh_G_x0, thresh_G_lw)
    return(dict_G)

def LO():
    dict_LO = {}
    dict_LO["initvalues"] = [initvalues_LO[0], initvalues_LO[1], 0, initvalues_LO[2]]
    dict_LO["peakname"] = "LO"
    dict_LO["fitrange"] = fitrange_LO
    dict_LO["plotrange"] = plotrange_LO
    dict_LO["params_thresh"] = (thresh_LO_c, thresh_LO_x0, thresh_LO_lw)
    return(dict_LO)

def TwoD():
    dict_2D = {}
    dict_2D["initvalues"] = [initvalues_2D[0], initvalues_2D[1], 0, initvalues_2D[2]]
    dict_2D["peakname"] = "2D"
    dict_2D["fitrange"] = fitrange_2D
    dict_2D["plotrange"] = plotrange_2D
    dict_2D["params_thresh"] = (thresh_2D_c, thresh_2D_x0, thresh_2D_lw)
    return(dict_2D)
