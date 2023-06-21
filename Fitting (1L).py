#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:24:38 2020

@author: Aaron
"""

import Library as lib


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

#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-02128_TA_R'"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-2128_TA_R'_detailed"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-02128_TA_detailed/hBN+Si signal substracted"
#folder = "//serveri2a/Transfer/Aaron/Work on tBLG/LowFreq Measurements/8_ETI-2128_flipped/8_ETI-02128_TA_detailed2/Si+hBN substracted"

#folder = "//serveri2a/Transfer/Nikita/TDBG screenschots/8 degrees JSHR04146 (D4 Gate) 200505/Raman/LOW FREQUENCY"
#folder = "//serveri2a/Transfer/Nikita/TDBG screenschots/8 degrees JSHR04146 (D4 Gate) 200505/Raman/Low Frequency 100nm"


###############################################################################
"""information about the sample"""
###############################################################################
filename = "7_ETI-02116_TA-detailed3"  # Name of the file
sample_size = (100,100)       # Size of the Scan
meanrange = (320,450)    # Data range for finding zero-line, in 1/cm

#dict_1 = lib.load_obj("dict_1", folder)

###############################################################################
"""Peak 1 fit parameters (typical TA parameters)"""
###############################################################################
try: type(dict_1)               # These two lines check if the peak 1 dictionary already exists and
except NameError: dict_1 = {}   # create it only if it doesn't. That way, stuff added later on won't be overwritten

dict_1["sample_size"] = sample_size
dict_1["peak_name"] = "TA"
dict_1["variables"] = [ "Offset", "Intensity", "Position", "Linewidth"]
dict_1["units"] = ["CCD counts", "a.u.", "1/cm", "1/cm"]
dict_1["fitrange"] = (240, 290)               # Data range for fitting
dict_1["plotrange"] = (240, 290)              # Data range for plotting
dict_1["startparams"] = [0, 250, 9999, 2]   # Lorentzian start values b, c, x0, lw
dict_1["background"] = "constant"


# Thresholds for (fitting) parameters. These are no boundary values, instead they will be used to sort out fits that yielded weird results due to poor data quality
b_thresh_1 = [-50, 50]         # Offset
c_thresh_1 = [40, 600]      # Intensity
x0_thresh_1 = [250, 275]      # Position
lw_thresh_1 = [0.4, 9]           # Linewidth
additional_thresh_1 = [0, 5]     # Optional additional threshold that has to be specified in lib.perform_fit1L
dict_1["params_thresh"] = (b_thresh_1, c_thresh_1, x0_thresh_1, lw_thresh_1, additional_thresh_1)
lib.save_obj("dict_1", dict_1, folder)


###############################################################################
"""Peak 2 fit parameters (typical G parameters)"""
###############################################################################
# =============================================================================
# try: type(dict_2)               # These two lines check if the peak 1 dictionary already exists and
# except NameError: dict_2 = {}   # create it only if it doesn't. That way, stuff added later on won't be overwritten
# 
# dict_2["sample_size"] = sample_size
# dict_2["peak_name"] = "G"
# dict_2["variables"] = [ "Offset", "Intensity", "Position", "Linewidth"]
# dict_2["units"] = ["CCD counts", "a.u.", "1/cm", "1/cm"]
# dict_2["fitrange"] = (1500, 1610)               # Data range for fitting
# dict_2["plotrange"] = (1500, 1610)              # Data range for plotting
# dict_2["startparams"] = [0, 1E11, 1582, 16]   # Lorentzian start values b, c, x0, lw
# dict_2["background"] = "constant"
# 
# 
# # Thresholds for (fitting) parameters. These are no boundary values, instead they will be used to sort out fits that yielded weird results due to poor data quality
# b_thresh_2 = [-200, 200]         # Offset
# c_thresh_2 = [1E3, 1E7]      # Intensity
# x0_thresh_2 = [1577, 1587]      # Position
# lw_thresh_2 = [4, 25]           # Linewidth
# lw_std_thresh_2 = [0,20]     # Covariance of linewidth
# dict_2["params_thresh"] = (b_thresh_2, c_thresh_2, x0_thresh_2, lw_thresh_2, lw_std_thresh_2)
# lib.save_obj("dict_2", dict_2, folder)
# =============================================================================


###############################################################################
"""Peak 3 fit parameters (typical R' parameters)"""
###############################################################################
# =============================================================================
# try: type(dict_3)               #These two lines check if the peak 1 dictionary already exists and
# except NameError: dict_3 = {}   #create it only if it doesn't. That way, stuff added later on won't be overwritten
# 
# dict_3["sample_size"] = sample_size
# dict_3["peak_name"] = "R'"
# dict_3["variables"] = ["Offset", "Intensity", "Position", "Linewidth"]
# dict_3["units"] = ["CCD Counts", "a.u.", "1/cm", "1/cm"]
# dict_3["fitrange"] = (1610, 1660)               # Data range for fitting
# dict_3["startparams"] = [0, 9E9, 1622, 5]      # Lorentzian start values b, c, x0, lw
# dict_3["plotrange"] = (1610, 1660)              # Data range for plotting
# dict_3["background"] = "constant"
# 
# # Thresholds for (fitting) parameters. These are no boundary values, instead they will be used to sort out fits that yielded weird results due to poor data quality
# c_thresh_3 = [5E8, 1E13]
# x0_thresh_3 = [1615, 1630]
# lw_thresh_3 = [1, 9]
# b_thresh_3 = [-100, 100]
# lw_std_thresh_3 = [0,10]    #covariance of linewidth
# dict_3["params_thresh"] = (b_thresh_3, c_thresh_3, x0_thresh_3, lw_thresh_3, lw_std_thresh_3)
# lib.save_obj("dict_3", dict_3, folder)
# =============================================================================


###############################################################################
"""Read data"""
###############################################################################
#xlist, data, lines = lib.read_raman_data(filename, sample_size, meanrange, folder)
#np.savetxt(str(folder) + "/" + "xlist.csv", xlist, delimiter=",", fmt='%s')


###############################################################################
"""Fit peak 1 and save the results."""
###############################################################################
fitresults_1, fitresults_std_1, fiterrors_1, dict_1 = lib.perform_fit1L(xlist, data, sample_size, dict_1)

lib.save_obj("fitresults_1", fitresults_1, folder)
lib.save_obj("fitresults_std_1", fitresults_std_1, folder)
lib.save_obj("fiterrors_1", fiterrors_1, folder)
lib.save_obj("dict_1", dict_1, folder)


###############################################################################
"""Fit peak 2 and save the results."""
###############################################################################
# =============================================================================
# mask = False                # Mask the fit so that it is performed only where the first fit succeeded
# mask_array = 0    # Array to use for the mask
# fitresults_2, fitresults_std_2, fiterrors_2, dict_2 = lib.perform_fit1L(xlist, data, sample_size, dict_2, mask, mask_array)
# lib.save_obj("fitresults_2", fitresults_2, folder)
# lib.save_obj("fitresults_std_2", fitresults_std_2, folder)
# lib.save_obj("fiterrors_2", fiterrors_2, folder)
# lib.save_obj("dict_2", dict_2, folder)
# =============================================================================

# =============================================================================
# ###############################################################################
# """Substract given peak from all datapoints, set baseline to zero again."""
# ###############################################################################
# data_substracted = lib.substract_fit(sample_size, xlist, data, fitresults_2)
# lib.save_obj("data_substracted", data_substracted)
# 
# 
# ###############################################################################
# """Fit peak 3 and save the results."""
# ###############################################################################
# mask = False                # Mask the fit so that it is performed only where the first fit succeeded
# mask_array = fiterrors_1    # Array to use for the mask
# fitresults_3, fitresults_std_3, fiterrors_3, dict_3 = lib.perform_fit1L(xlist, data_substracted, sample_size, dict_3, mask, mask_array)
# lib.save_obj("fitresults_3", fitresults_3, folder)
# lib.save_obj("fitresults_std_3", fitresults_std_3, folder)
# lib.save_obj("fiterrors_3", fiterrors_3, folder)
# lib.save_obj("dict_3", dict_3, folder)
# =============================================================================
