#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:21:56 2020

@author: Aaron
"""

import Library as lib 
import numpy as np
import numpy.ma as ma
from numpy import sqrt
import os

###############################################################################
"""Set the working directory"""
###############################################################################

#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-09-30 TA-R'-Map nice (new Lorentzian)"
folder = "//serveri2a/Transfer/Aaron/Twisted Bilayer Graphene/HiWi/LowFreq Measurements/7_ETI-02116_flipped/2020-10-17 TA-Map detailed3 (new Lorentzian)"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/7_ETI-02116_flipped/2020-10-16 TA-Map detailed2 (new Lorentzian)"

#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/Samples/Samples flipped on SiO2/7,5_ETIHR02-120_flipped/Raman"

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
"""Create peak 1 maps."""
"""Indexes: 0-Intensity, 1-Position, 2-Linewidth, 3-Height."""
###############################################################################
fitresults_1 = lib.load_obj("fitresults_1", folder)
fitresults_std_1 = lib.load_obj("fitresults_std_1", folder)
fiterrors_1 = lib.load_obj("fiterrors_1", folder)
dict_1 = lib.load_obj("dict_1", folder)
sample_size = dict_1["sample_size"]


# =============================================================================
# fitresults_1 = lib.load_obj("fitresults_1_stdlw0,5", folder)
# fitresults_std_1 = lib.load_obj("fitresults_std_1_stdlw0,5", folder)
# fiterrors_1 = lib.load_obj("fiterrors_1_stdlw0,5", folder)
# dict_1 = lib.load_obj("dict_1_stdlw0,5", folder)
# maps_images_1 = lib.load_obj("maps_images_1_stdlw0,5", folder)
# =============================================================================


# =============================================================================
# print("std lw mean: ", np.mean(fitresults_std_1[:,:,3]))
# print("std lw median: ", np.median(fitresults_std_1[:,:,3]))
# =============================================================================

plot_size = (6,6)
save = False  # Whether or not the maps shall be saved
title = 0

# Get map data of c, x0, lw from fitresults
maps_data_1 = lib.create_mapsarray3D(fitresults_1, fiterrors_1, dict_1)

# Create maps from map data
maps_images_1 = lib.create_maps(dict_1, maps_data_1, fiterrors_1, plot_size, save, folder=folder, title=title)
lib.save_obj("maps_images_1", maps_images_1, folder)

# =============================================================================
# # For TA-peak: Create map of the twist angle
# theta_array_1 = lib.create_theta_map(dict_1, fitresults_1, fiterrors_1, plot_size, save, folder=folder, title=title)
# lib.save_obj("theta_array_1", theta_array_1, folder)
# # Check given datapoint in detail, inluding plot of the fit
# hx, hy = 72,115   # Coordinates of the datapoint
# map_index = 1       # 0 - intensity, 1 - position, 2 - linewidth
# #lib.check_peak(xlist, data, fitresults_1, fiterrors_1, sample_size, plot_size, dict_1, maps_data_1, map_index, "TA", hx, hy)
# 
# =============================================================================

###############################################################################
"""Create peak 2 maps."""
"""Indexes: 0-Intensity, 1-Position, 2-Linewidth, 3-Height."""
###############################################################################
# =============================================================================
# fitresults_2 = lib.load_obj("fitresults_2", folder)
# fitresults_std_2 = lib.load_obj("fitresults_std_2", folder)
# fiterrors_2 = lib.load_obj("fiterrors_2", folder)
# dict_2 = lib.load_obj("dict_2", folder)
# 
# plot_size = (4,8)
# save = False    #Whether or not the maps shall be saved
# 
# # Get map data of c, x0, lw from fitresults
# maps_data_2 = lib.create_mapsarray3D(fitresults_2, fiterrors_2, dict_2)
# 
# # Create maps from map data
# maps_images_2 = lib.create_maps(dict_2, maps_data_2, fiterrors_2, plot_size, save, folder=folder)
# #lib.save_obj("maps_images_2", maps_images_2, folder)
# 
# 
# #Check given datapoint in detail, inluding plot of the fit
# hx, hy = 30,100  #Coordinates of the datapoint
# map_index = 2       # 0 - intensity, 1 - position, 2 - linewidth
# #lib.check_peak(xlist, data, fitresults_2, fiterrors_2, sample_size, plot_size, dict_2, maps_data_2, map_index, "G", hx, hy)
# =============================================================================

###############################################################################
"""Create peak 3 maps."""
"""Indexes: 0-Intensity, 1-Position, 2-Linewidth, 3-Height."""
###############################################################################
# =============================================================================
# fitresults_3 = lib.load_obj("fitresults_3", folder)
# fitresults_std_3 = lib.load_obj("fitresults_std_3", folder)
# fiterrors_3 = lib.load_obj("fiterrors_3", folder)
# dict_3 = lib.load_obj("dict_3", folder)
# 
# plot_size = (8,4)
# save = False    # Whether or not the maps shall be saved
# 
# # Get map data of c, x0, lw from fitresults
# maps_data_3 = lib.create_mapsarray3D(fitresults_3, fiterrors_3, dict_3)
# 
# # Create maps from map data
# maps_images_3 = lib.create_maps(dict_3, maps_data_3, fiterrors_3, plot_size, save, folder=folder)
# lib.save_obj("maps_images_3", maps_images_3, folder)
# 
# 
# # Check given datapoint in detail, inluding plot of the fit
# hx, hy = 30,15  # Coordinates of the datapoint
# map_index = 2   # 0 - intensity, 1 - position, 2 - linewidth
# #lib.check_peak(xlist, data_substracted, fitresults_3, fiterrors_3, sample_size, plot_size, dict_3, maps_data_3, map_index, "R'", hx, hy)
# 
# =============================================================================

###############################################################################
"""Map the distance between two peaks"""
###############################################################################
# =============================================================================
# dict1 = dict_2
# fitdata1 = fitresults_2
# errordata1 = fiterrors_2
# 
# fitdata2 = fitresults_3
# errordata2 = fiterrors_3
# dict2 = dict_3
# 
# variable_index = 2  # 1 - intensity, 2 - position, 3 - linewidth
# save = True
# plot_size = (8,4)
# 
# masked_difference = lib.difference_map(dict1, dict2, fitdata1, fitdata2, errordata1, errordata2, variable_index, plot_size, save, folder) 
# lib.save_obj("masked_difference"+str(dict1["peak_name"])+str(dict2["peak_name"]), masked_difference, folder)
# =============================================================================


###############################################################################
"""Map the ratio between two peak parameters"""
###############################################################################
# =============================================================================
# save = True
# variable_index = 1  # 1 - intensity, 2 - position, 3 - linewidth
# plot_size = (8,6)
# masked_ratio = lib.ratio_map(dict_1, dict_2, fitresults_1, fitresults_2, fiterrors_1, fiterrors_2, variable_index, plot_size, save, folder)
# =============================================================================

###############################################################################
"""Map a fit parameter gradient"""
###############################################################################
# =============================================================================
# # Settings for data processing
# dict_grad = dict_1
# #data_grad = fitresults_1[:,:,2]
# data_grad = theta_array
# fiterrors_grad = fiterrors_1
# reach = 1
# stepsize = 100*0.001 # In micrometer
# 
# # Settings for the plot
# figsize = 4,4
# title = "Map Grad(Theta)"
# save = False
# 
# gradient_abs_av, exvalues = lib.averaged_gradient(data_grad, fiterrors_grad, reach, stepsize)
# lib.save_obj("gradient_abs_av", gradient_abs_av, folder)
# 
# maps_data_grad = lib.create_mapsarray2D(gradient_abs_av, dict_grad, exvalues)
# gradient_mask = gradient_abs_av.mask
# maps_images_grad = lib.create_maps(dict_grad, maps_data_grad, gradient_mask, figsize, save, map_index="none", title=title, exvalues=exvalues, folder=folder)
# 
# =============================================================================

###############################################################################
"""Nice maps for papers"""
###############################################################################
# =============================================================================
# pdict = dict_2
# mapdata = fitresults_2[:,:,2]
# #mapdata, _ = lib.averaged_gradient(data_grad, fiterrors_grad, reach, stepsize)
# mask = True
# errordata = fiterrors_2
# thresholds = (1578,1586)
# 
# mapdata = np.rot90(mapdata)
# errordata = np.rot90(errordata)
# 
# n_ticks = 5
# n_decimals = 0
# 
# scale2 = False
# data2 = fitresults_1[:,:,3]
# n_ticks2 = 5
# n_decimals2 = 0
# 
# plotsize = (6,3)
# fontsize = 12
# title = "Pos(G)"
# save = False
# 
# lib.papermap(pdict, mapdata, mask, errordata, thresholds, plotsize, n_ticks, n_decimals, fontsize, scale2, data2, n_ticks2, n_decimals2, save, title, folder)
# 
# =============================================================================

