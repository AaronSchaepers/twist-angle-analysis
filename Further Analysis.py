#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:31:29 2020

@author: Aaron
"""

import numpy as np
import numpy.ma as ma
from numpy import sqrt
import Library as lib 

folder = "//serveri2a/Transfer/Aaron/Twisted Bilayer Graphene/HiWi/LowFreq Measurements/7_ETI-02116_flipped/2020-10-17 TA-Map detailed3 (new Lorentzian)"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-09-30 TA-R'-Map nice (new Lorentzian)"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-10-16 TA-Map detailed2 (new Lorentzian)"


#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/6_ETI-0292_flipped"

#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-09-30 TA-R'-Map nice"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-10-02 TA-Map detailed"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-10-16 TA-Map detailed2"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-10-17 TA-Map detailed3"

#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/Samples/Samples flipped on SiO2/8_ETI-02128_flipped/Raman"

#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_JSHR05-46/JSHR05-46_8deg_5_nice"

#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-02128_TA_R'"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-02128_TA_R'_detailed"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-02128_TA_detailed/hBN+Si signal substracted"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-02128_TA_detailed2/Si+hBN substracted"

#folder = "//serveri2a/Transfer/Nikita/TDBG screenschots/8 degrees JSHR04146 (D4 Gate) 200505/Raman/LOW FREQUENCY"
#folder = "//serveri2a/Transfer/Nikita/TDBG screenschots/8 degrees JSHR04146 (D4 Gate) 200505/Raman/Low Frequency 100nm"


###############################################################################
"""Load data from the fits"""
###############################################################################
fitresults_1 = lib.load_obj("fitresults_1", folder)
fitresults_std_1 = lib.load_obj("fitresults_std_1", folder)
fiterrors_1 = lib.load_obj("fiterrors_1", folder)
dict_1 = lib.load_obj("dict_1", folder)
maps_images_1 = lib.load_obj("maps_images_1", folder)

try:
    theta_array_1 = lib.load_obj("theta_array_1", folder)
except:
    pass
    

# Rotate the images if necessary
# =============================================================================
# fitresults_1 = np.rot90(fitresults_1)
# fiterrors_1 = np.rot90(fiterrors_1)
# =============================================================================


# =============================================================================
# fitresults_2 = lib.load_obj("fitresults_2", folder)
# fitresults_std_2 = lib.load_obj("fitresults_std_2", folder)
# fiterrors_2 = lib.load_obj("fiterrors_2", folder)
# dict_2 = lib.load_obj("dict_2", folder)
# maps_images_2 = lib.load_obj("maps_images_2", folder)
# 
# fitresults_3 = lib.load_obj("fitresults_3", folder)
# fitresults_std_3 = lib.load_obj("fitresults_std_3", folder)
# fiterrors_3 = lib.load_obj("fiterrors_3", folder)
# dict_3 = lib.load_obj("dict_3", folder)
# maps_images_3 = lib.load_obj("maps_images_2", folder)
# =============================================================================

###############################################################################
"""Create histograms."""
"""Indexes: 0-Intensity, 1-Position, 2-Linewidth, 3-Height."""
###############################################################################
# =============================================================================
# # Histograms
# hist_data = fitresults_1[:,:,3]*2
# error_data = fiterrors_1
# bins = 50
# title = "Histogram TA FWHM"
# save = False
# lib.histogram(hist_data, error_data, bins, save, title, folder)
# =============================================================================


###############################################################################
"""Extract mean values of different quantities."""
"""Indices: 0 - offset, 1 - intensity, 2 - position, 3 - linewidth"""
###############################################################################
fitdata = fitresults_std_1      # Array from which to take the averaged quantity
errordata = fiterrors_1     # The corresponding fiterror array
mean_index = 2              # Index of the averaged quantity in the array
pdict = dict_1              # The dictionary of the corresponding peak
maps_images = maps_images_1 # Array containing the maps 

# With these parameters, the area that will be considered when calculating the mean can be adjusted
x_bounds = (0,100)              # Specify interval of relevant x coordinates
y_bounds = (0,100)              # Specify interval of relevant y coordinates
# =============================================================================
# x_bounds = (00,85)              # Specify interval of relevant x coordinates
# y_bounds = (57,72)              # Specify interval of relevant y coordinates
# 
# =============================================================================

# Define an additional threshold condition that a data point has to fulfill to be considered
cond_array = fitresults_std_1   # Define array from which to take a threshold parameter
cond_var = 2                # Index of the threshold parameter in the array (Intensity - 0, Position - 1, Linewidth - 2)
lb, ub = -np.inf, np.inf       # Upper and lower threshold

save = False               # Whether or not the result shall be saved
mean, std = lib.get_mean(fitdata, errordata, mean_index, pdict, maps_images, x_bounds, y_bounds, cond_array, cond_var, lb, ub, save, folder)


###############################################################################
"""Export position values for external use (linear correlation fit)."""
###############################################################################
# =============================================================================
# Rp_pos = fitresults_1[:,:,2]
# TA_pos = fitresults_2[:,:,2]
# np.savetxt(str(folder) + "/" + "Rp_pos.csv", Rp_pos, delimiter="\t", )
# np.savetxt(str(folder) + "/" + "TA_pos.csv", TA_pos, delimiter="\t", )
# =============================================================================


###############################################################################
"""Fancy scatter plots to reveal correlations between two quantities."""
"""Indexes: 0-Intensity, 1-Position, 2-Linewidth, 3-Height.    """
###############################################################################
# =============================================================================
# data1 = fitresults_1    #Array from which to take the first quantity
# index1 = 2             #Index of the first parameter
# data2 = fitresults_2    #Array from which to take the second quantity
# index2 = 2              #Index of the second parameter
# title = "TA over R' position"
# save = True         #Whether or not to save the scatter plot
# #lib.fancy_scatter_plot(data1, data2, index1, index2, dict_1, dict_2, fiterrors_2, title, save, folder=folder)
# =============================================================================


###############################################################################
"""Compute the gradient of any given fit parameter"""
###############################################################################
# =============================================================================
# # Settings for gradient data processing
# dict_grad = dict_1
# data_grad = fitresults_1[:,:,2]
# fiterrors_grad = fiterrors_1
# reach = 2
# stepsize = 70
# 
# gradient_abs_av, gradient_mask, exvalues = lib.averaged_gradient(data_grad, fiterrors_grad, reach, stepsize)
# =============================================================================




###############################################################################
"""Scatter plot to reveal correlations between two quantities, optionally color coded"""
"""Indexes for fitresults: 0 - offset, 1 - intensity, 2 - position, 3 - linewidth"""
"""Indexes for images: 0 - intensity, 1 - position, 2 - linewidth"""
###############################################################################
# =============================================================================
# # General ettings for the plot
# figsize = 3,4
# title = "Lw(TA) vs Grad(Theta) (2400gr_mm) 1"
# xlabel = r"Grad(Theta) (Degree/$\mu$m))"
# ylabel = "FWHM(TA) (1/cm)"
# save = False
# 
# # Color code settings
# color_code = False
# color_data = maps_images_1[:,:,0] # Array that serves for color coding the data points
# color_threshold = 0               # Threshold that can be applied to rule data points out according to their value in the color_data
# 
# # =============================================================================
# # # x-data
# # data_x = fitresults_2[:,:,2]
# # xlim = (1575,1587)
# # 
# # # y-data
# # data_y = fitresults_1[:,:,2]
# # ylim = (240,280)
# # 
# # =============================================================================
# 
# # x-data: Averaged gradient of position
# stepsize= 100*0.001 # From nm to mu
# reach = 1
# data_x, _ = lib.averaged_gradient(theta_array_1, fiterrors_1, reach, stepsize)
# xlim = (0, 0.1/stepsize)
# #data_x = np.rot90(data_x, 3)
# 
# # y-data: Averaged linewidth
# mask = True
# mask_array = fiterrors_1
# data_y = lib.moving_average(fitresults_1[:,:,3], reach, mask, mask_array)
# ylim = (1,7)
# 
# # Errorbar settings
# errorbars = False
# std_data = fitresults_std_1[:,:,3]
# 
# # Whether or not to include the theoretical curve in the plot
# withgradprediction = True
# 
# # =============================================================================
# # data_x = data_x[::2, ::2]
# # data_y = data_y[::2, ::2]
# # std_data = std_data[::2, ::2]
# # =============================================================================
# data_y*=2
# 
# lib.scatter_plot(data_x, data_y, figsize, xlim, ylim, title, xlabel, ylabel, color_code, color_data, color_threshold, errorbars, std_data, withgradprediction, save, folder)
# 
# =============================================================================



###############################################################################
"""Calculate int(R')/int(G) for every datapoint."""
"""Indexes: 0-Intensity, 1-Position, 2-Linewidth, 3-Height."""
###############################################################################
# =============================================================================
# data_1 = fitresults_1       # Data set from which to take the dividend
# data_2 = fitresults_2       # Data set from which to take the divisor
# ratio_index_1 = 0           # Index of the dividend
# ratio_index_2 = 0           # Index of the divisor
# height_threshold = 100      # Height that the "divident peak" must have in a given data point for the data point to be considered
# 
# #ratio_mean, ratio_std = lib.get_ratio(sample_size, data_1, data_2, ratio_index_1, ratio_index_2, height_threshold)
# =============================================================================
