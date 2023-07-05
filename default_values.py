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

# Starting parameters for Lorentzian peak fits in the order: Intensity, linewidth.
# The starting value for the position is determined dynamically, that for the 
# offset is hardcoded to zero
initvalues_TA = [250, 2]
initvalues_G = [1E4, 16]
initvalues_LO = [5E3, 4]
initvalues_2D = [2E4, 15]

# Threshold parameters for excluding implausible fit results
# NOTE: These are not boundary values for the fitting routine!
#       Instead, the best fit parameters of each spectrum are compared to these  
#       threshold intervals. If any parameter lies outside of its threshold
#       interval, the corresponding data point is excluded from the map.

thresh_TA_int = [20, 1E4]      # Intensity
thresh_TA_pos = [250, 275]    # Position
thresh_TA_lw = [0.2, 9]      # Linewidth

thresh_G_int = [1E3, 1E15]     # Intensity
thresh_G_pos = [1577, 1587]   # Position
thresh_G_lw = [4, 25]        # Linewidth

thresh_LO_int = [300, 1E5]     # Intensity
thresh_LO_pos = [1610, 1630]  # Position
thresh_LO_lw = [0.4, 15]     # Linewidth

thresh_2D_int = [1E3, 1E7]     # Intensity
thresh_2D_pos = [2650, 2700]  # Position
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
    # determined dynamically in the fitting routine.
    # Hardcode the offset initial value to 0.
    dict_TA["init_off"] = 0
    dict_TA["init_int"] = initvalues_TA[0]
    dict_TA["init_pos"] = 0
    dict_TA["init_lw"] = initvalues_TA[1]
    dict_TA["min_int"] = thresh_TA_int[0]
    dict_TA["max_int"] = thresh_TA_int[1]
    dict_TA["min_pos"] = thresh_TA_pos[0]
    dict_TA["max_pos"] = thresh_TA_pos[1]
    dict_TA["min_lw"] = thresh_TA_lw[0]
    dict_TA["max_lw"] = thresh_TA_lw[1]
    dict_TA["peakname"] = "TA"
    dict_TA["fitrange"] = fitrange_TA
    dict_TA["plotrange"] = plotrange_TA
    dict_TA["max_gradient"] = max_gradient
    return(dict_TA)

def G():
    dict_G = {}
    dict_G["init_off"] = 0
    dict_G["init_int"] = initvalues_G[0]
    dict_G["init_pos"] = 0
    dict_G["init_lw"] = initvalues_G[1]
    dict_G["min_int"] = thresh_G_int[0]
    dict_G["max_int"] = thresh_G_int[1]
    dict_G["min_pos"] = thresh_G_pos[0]
    dict_G["max_pos"] = thresh_G_pos[1]
    dict_G["min_lw"] = thresh_G_lw[0]
    dict_G["max_lw"] = thresh_G_lw[1]
    dict_G["peakname"] = "G"
    dict_G["fitrange"] = fitrange_G
    dict_G["plotrange"] = plotrange_G
    return(dict_G)

def LO():
    dict_LO = {}
    dict_LO["init_off"] = 0
    dict_LO["init_int"] = initvalues_LO[0]
    dict_LO["init_pos"] = 0
    dict_LO["init_lw"] = initvalues_LO[1]
    dict_LO["min_int"] = thresh_LO_int[0]
    dict_LO["max_int"] = thresh_LO_int[1]
    dict_LO["min_pos"] = thresh_LO_pos[0]
    dict_LO["max_pos"] = thresh_LO_pos[1]
    dict_LO["min_lw"] = thresh_LO_lw[0]
    dict_LO["max_lw"] = thresh_LO_lw[1]
    dict_LO["peakname"] = "LO"
    dict_LO["fitrange"] = fitrange_LO
    dict_LO["plotrange"] = plotrange_LO
    return(dict_LO)

def TwoD():
    dict_2D = {}
    dict_2D["init_off"] = 0
    dict_2D["init_int"] = initvalues_2D[0]
    dict_2D["init_pos"] = 0
    dict_2D["init_lw"] = initvalues_2D[1]
    dict_2D["min_int"] = thresh_2D_int[0]
    dict_2D["max_int"] = thresh_2D_int[1]
    dict_2D["min_pos"] = thresh_2D_pos[0]
    dict_2D["max_pos"] = thresh_2D_pos[1]
    dict_2D["min_lw"] = thresh_2D_lw[0]
    dict_2D["max_lw"] = thresh_2D_lw[1]
    dict_2D["peakname"] = "2D"
    dict_2D["fitrange"] = fitrange_2D
    dict_2D["plotrange"] = plotrange_2D
    dict_2D["params_thresh"] = (thresh_2D_int, thresh_2D_pos, thresh_2D_lw)
    return(dict_2D)
