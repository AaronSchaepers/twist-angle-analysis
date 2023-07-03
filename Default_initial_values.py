#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 19:22:07 2023

@author: Aaron
"""


# Starting parameters for Lorentzian peak fits in the order:
# Offset, intensity, linewidth.
# The starting value for the position is determined dynamically
startparams_TA = [0, 250, 2]
startparams_G = [0, 1E4, 16]
startparams_LO = [0, 5E3, 4]
startparams_2D = [0, 2E4, 15]

# Threshold parameters for excluding implausible fit results
# NOTE: These are not boundary values for the fitting routine!
#       Instead, the best fit parameters of each spectrum are compared to these  
#       threshold intervals. If any parameter lies outside of its threshold
#       interval, the corresponding data point is excluded from the map.

thresh_TA_c = [20, 6000]         # Intensity
thresh_TA_x0 = [250, 275]   # Position
thresh_TA_lw = [0.2, 9]         # Linewidth

thresh_G_c = [1E3, 1E15]       # Intensity
thresh_G_x0 = [1577, 1587]    # Position
thresh_G_lw = [4, 25]         # Linewidth

thresh_LO_c = [300, 1E5]      # Intensity
thresh_LO_x0 = [1610, 1630]   # Position
thresh_LO_lw = [0.4, 15]      # Linewidth

thresh_2D_c = [1E3, 1E7]      # Intensity
thresh_2D_x0 = [2650, 2700]   # Position
thresh_2D_lw = [8, 30]        # Linewidth