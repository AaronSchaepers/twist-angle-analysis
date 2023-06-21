#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Aaron
@date: 2023/06/19

PREPARATION

    Export Raman maps from ProjectFive as .txt by choosing "Table" in the
    export menu. Make sure to set the spectral unit to 1/cm.


HOW TO USE THIS CODE (SHORT VERSION)

    Enter the required information about your data in section 1 and run 
    this program to start the fitting and mapping procedure.
    
    
HOW TO USE THIS CODE (LONG VERSION)

    This code has three sections:
        1. User input section
        2. Advanced user input section
        3. Don't touch section
        
    In these sections, you need to do the following:
    1. User input section
        - Provide the required information about the scan (see comments)
        - Choose which peaks to map
        
    2. Advanced user input section
        - For each peak, set the starting parameters for the Lorentzian fit
        - For each fitting parameter and each peak, define a threshold interval
          within which the fit results are accepted as realistic
          
    3. Don't touch section
        - You know what to do with this one
"""

import lib
import numpy as np

###############################################################################
# 1. User input
###############################################################################

# Directory of the Raman data
folder = "/Users/Aaron/Desktop/Data_ETIRF04-100/Test_data" 

# Name of the .txt file containing the Raman data, given with suffix
file = "test_data.txt"  

# 1/cm, a spectral range without any features. This is used to calculate 
# the mean background noise which is subtracted from the data
spectral_mean_range = (320, 450) 

size_px = (100,100)       # Size of the Scan in pixels
size_um = (10, 10)        # Size of the Scan in µm

# What peaks shall be fitted?
b_fit_TA = False
b_fit_G = False
b_fit_LO = False
b_fit_2D = False

# What peaks shall be mapped?
b_map_TA = True
b_map_G = True
b_map_LO = True
b_map_2D = True


###############################################################################
# 2. Advanced user input
###############################################################################

# Starting parameters for Lorentzian peak fits in the order:
# Offset, intensity, linewidth
startparams_TA = [0, 250, 2]
startparams_G = [0, 1E11, 16]
startparams_LO = 0
startparams_2D = 0

# Threshold parameters for excluding implausible fit results
# NOTE: These are not boundary values for the fitting routine!
#       Instead, the best fit parameters of each spectrum are compared to these  
#       threshold intervals. If any parameter lies outside of its threshold
#       interval, the corresponding data point is excluded from the map.

thresh_TA_c = [40, 600]         # Intensity
thresh_TA_x0 = [250, 275]   # Position
thresh_TA_lw = [0.4, 9]         # Linewidth
thresh_TA_lw_std = [0,5]     # Covariance of linewidth

thresh_G_c = [1E3, 1E7]       # Intensity
thresh_G_x0 = [1577, 1587]    # Position
thresh_G_lw = [4, 25]         # Linewidth
thresh_G_lw_std = [0,20]      # Covariance of linewidth

thresh_LO_c = [1E3, 1E7]      # Intensity
thresh_LO_x0 = [1577, 1587]   # Position
thresh_LO_lw = [4, 25]        # Linewidth
thresh_LO_lw_std = [0,20]     # Covariance of linewidth

thresh_2D_c = [1E3, 1E7]      # Intensity
thresh_2D_x0 = [1577, 1587]   # Position
thresh_2D_lw = [4, 25]        # Linewidth
thresh_2D_lw_std = [0,20]     # Covariance of linewidth

max_gradient = 1 # °/µm, upper bound for the twist angle gradient map. Values
                 #       below will be excluded from the map.


###############################################################################
# 3. Don't touch 
###############################################################################

""" Assemble peak dictionaries with information for mapping and fitting """

# The first two lines check if the peak dictionary already exists and create
# it only if it doesn't. That way, stuff added later on won't be overwritten.
try: type(dict_TA)
except NameError: dict_TA = {}
# Insert a place holder value for the position starting value, which is 
# determined dynamically in the fitting routine
dict_TA["startparams"] = [startparams_TA[0], startparams_TA[1], 0, startparams_TA[2]]
dict_TA["size_px"] = size_px
dict_TA["size_um"] = size_um
dict_TA["peakname"] = "TA"
dict_TA["fitrange"] = (240, 290)             # Data range for fitting
dict_TA["plotrange"] = (240, 290)            # Data range for plotting
dict_TA["params_thresh"] = (thresh_TA_c, thresh_TA_x0, thresh_TA_lw, thresh_TA_lw_std)
dict_TA["max_gradient"] = max_gradient
lib.save_object(folder, dict_TA, "dict_TA")


# The first two lines check if the peak dictionary already exists and create
# it only if it doesn't. That way, stuff added later on won't be overwritten.
try: type(dict_G)
except NameError: dict_G = {}
dict_G["startparams"] = startparams_G
dict_G["size_px"] = size_px
dict_G["size_um"] = size_um
dict_G["peakname"] = "G"
dict_G["fitrange"] = (1500, 1610)             # Data range for fitting
dict_G["plotrange"] = (1500, 1610)            # Data range for plotting
dict_G["params_thresh"] = (thresh_G_c, thresh_G_x0, thresh_G_lw, thresh_G_lw_std)
lib.save_object(folder, dict_G, "dict_G")


""" Read data """

xdata, data = lib.read_raman_scan(folder, file, size_px, spectral_mean_range)
lib.save_object(folder, xdata, "xdata")

""" Do the fitting and mapping for each peak """

if b_fit_TA == True:
    fitresults_TA, fitresults_std_TA, fiterrors_TA = lib.fit_to_map(xdata, data, dict_TA)
      
if b_map_TA == True:
    lib.map_lorentz_parameters(fitresults_TA, fiterrors_TA, dict_TA, folder)
    lib.map_theta(fitresults_TA, fiterrors_TA, dict_TA, folder)
    
    
    
    




