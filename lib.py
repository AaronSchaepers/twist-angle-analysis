#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:10:18 2023

@author: Aaron
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi, arcsin, sqrt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

a = 0.246 # Graphene lattice constant in nm


###############################################################################
""" 1. Functions for importing and exporting data """
###############################################################################

""" Save an object using the pickle module """
def save_object(folder, file, name):

    with open(str(folder) + "/" + str(name)+".pkl", "wb") as f:
            pickle.dump(file, f, pickle.HIGHEST_PROTOCOL)
         
    return()


""" Load an object using the pickle module """
def load_object(folder, file):
    
    with open(str(folder)+"/"+str(file)+".pkl", "rb") as f:
            return(pickle.load(f))


""" Read Raman scan data from .txt file, set baseline to zero """
def read_raman_scan(folder, file, size_px, spectral_mean_range):
    nx, ny = size_px # Size of the scan in pixels
    
    # A tuple of strings with 1600 entries, i.e., one entry per CCD pixel, i.e., one entry
    # per data point in the Raman spectrum.
    # Each of the 1600 strings is a succession of numbers giving the CCD counts
    # on the corresponding pixel for all the spatial data points measured in the scan.
    lines = tuple(open(folder + "/"+ file, 'r')) 
    
    # A list containing 1600 sublists corresponding to the 1600 CCD pixels.
    # Each sublist contains nx*ny strings giving the the CCD counts on that 
    # pixel for all the spatial data points.
    data_str = [x.split('\t') for x in lines]
    
    # A list containing nx*ny+1 sublists with 1600 entries each. The first list 
    # contains the wave numbers, the other lists are the measured spectra.
    data_float = [[float(liste[i]) for liste in data_str ] for i in range(len(data_str[0]))]
    
    xlist = data_float.pop(0) # Extract first list which contains the x-data
    
    # A 3D array with the size nx*ny*1600, i.e., the first two axes correspond
    # to the x and y axes of the scan and the spectra will be stacked along the
    # third axis.
    data = np.zeros((nx,ny,len(data_str)))
    
    # Get indices of the averaging interval for baseline shifting
    i_xmin_mean = min(range(len(xlist)), key=lambda i: abs(xlist[i]-spectral_mean_range[0]))
    i_xmax_mean = min(range(len(xlist)), key=lambda i: abs(xlist[i]-spectral_mean_range[1]))
    
    # Scan over all spatial data points
    x,y = 0,0
    for i in range(nx*ny):
        # Calculate the average CCD counts in a range without peaks
        ymean = np.mean(data_float[i][i_xmin_mean:i_xmax_mean])   
        # Subtract that average from the entire spectrum and stack the spectrum
        # into the data array
        data[x,y] = [y-ymean for y in data_float[i]]        
        x = x+1
        if x == nx:
            y = y+1
            x = 0
            
    # Convert xlist to array for more efficient function evaluation
    xdata = np.asarray(xlist)
    
    return(xdata, data)

###############################################################################
""" 2. Functions for the fitting procedure """
###############################################################################

""" Lorentzian as used in spectroscopy """
def lorentzian(x, b, c, x0, lw):
   f = b + c / pi * lw / ((x-x0)**2 + lw**2)
   return(f)
    


""" Perform a Lorentzian fit across a map of spectra """
def fit_to_map(xdata, data, pdict):
    
    # Retrieve required variables from the dictionnary
    nx, ny = pdict["size_px"] # Scan size in pixels
    fitrange = pdict["fitrange"] # Spectral range for the fit
    p0 = pdict["startparams"] # Starting values
    (thresh_c, thresh_x0, thresh_lw, thresh_lw_std) = pdict["params_thresh"]
    
    # Create arrays to store the fit results
    fitresults = np.zeros((nx,ny,3))          # Array with (b, c, x0, lw) for each x/y-datapoint in which fit was successful
    fitresults_std = np.zeros((nx,ny,3))      # Array with the uncertainties of the fit results
    fiterrors = np.zeros((nx,ny), dtype= "object")    # Array with inormation if and why the fit failed for each x/y-datapoint
    
    # This expression searches the xlist for the wave number values 
    # closest to the chosen fitrange and retrieves their indexes
    i_start = min(range(len(xdata)), key=lambda i: abs(xdata[i]-fitrange[0]))
    i_stop = min(range(len(xdata)), key=lambda i: abs(xdata[i]-fitrange[1]))
    
    # Slice xdata to keep only the relevant peak
    xdata_fit = xdata[i_start:i_stop]
    
    # Define boundaries to exclude negative values
    lbounds = np.array((-np.inf, 0, 0, 0))
    ubounds = np.array((np.inf, np.inf, np.inf, np.inf))
    
    print("Fitting progress:")
    for x in range(nx):
        progress = np.round(x/nx*100)
        print(f"{progress}%")
        for y in range(ny):
            
            # Retrieve the spectrum of the current scan data point
            ydata = data[y,x]
            
            # Slice ydata to keep only the relevant peak
            ydata_fit = ydata[i_start:i_stop]
            
            # Determine highest data point, use position as x0 starting value
            p0[2] = xdata_fit[np.argmax(ydata_fit)]
            
            # Try to perform the fit.
            # If the fit fails, put a note in the fitresults array and proceed.
            try:
                popt, pcov = curve_fit(lorentzian, xdata_fit, ydata_fit, p0, bounds=(lbounds, ubounds))
                popt_std = np.sqrt(np.diag(pcov))
                
                # Check that the results do not contain infinities
                if None in popt:
                    fiterrors[y,x] = "Nan val"
                    continue
                if None in popt_std:
                    fiterrors[y,x] = "Nan std"
                    continue
                
                # Check if the fit results fall within their respective threshold
                # intervals.
                # If yes: Write fit results into result arrays
                # If not: Change the fiterrors entry to show which parameter(s) 
                # violated the threshold interval
                if not thresh_c[0] < popt[1] < thresh_c[1]:
                    fiterrors[y,x] = "c"
                    
                elif not thresh_x0[0] < popt[2] < thresh_x0[1]:
                    fiterrors[y,x] = "x0"
                    
                elif not thresh_lw[0] < popt[3] < thresh_lw[1]:
                    fiterrors[y,x] = "lw"
                    
                elif not thresh_lw_std[0] < popt_std[3] < thresh_lw_std[1]:
                    fiterrors[y,x] = "lw std"
                    
                else:
                    # Drop the first parameter (offset) as it is of no further
                    # interest for Raman maps
                    fitresults[y,x] = popt[1:]
                    fitresults_std[y,x] = popt_std[1:]
                    
            except: 
                fiterrors[y, x] = "fit fail"
                
    return(fitresults, fitresults_std, fiterrors)



###############################################################################
""" 3. Functions for the mapping """
###############################################################################

""" Map the intensity, position and linewidth of a given peak """
def map_lorentz_parameters(fitresults, fiterrors, pdict, folder):
    # Retrieve required variables from the dictionnary
    nx, ny = pdict["size_px"] # Scan size in pixels
    sx, sy = pdict["size_um"] # Scan size in microns
    peakname = pdict["peakname"] # Name of the peak
    
    # Mask the fitresults array in all points where the fit failed
    # Define the condition for masking
    mask_condition = fiterrors != 0
    # Create the mask, expand it to match the shape of fitresults
    expanded_mask = np.repeat(mask_condition[:, :, np.newaxis], fitresults.shape[2], axis=2)
    # Apply the mask to fitresults
    fitresults_ma = np.ma.array(fitresults, mask=expanded_mask)
        
    # The mapped quantities
    quantities = ["intensity", "position", "linewidth"]
    # Their units
    cbar_labels = ["Intensity (arb. u.)", r"$\omega$ (1/cm)", r"$\Gamma$ (1/cm)"]
    
    # In a loop, plot peak intensity, position and linewidth
    for i in range(3):
        im = plt.imshow(fitresults_ma[:,:,i], extent = [0, sx, 0, sy], cmap="gist_rainbow")
        plt.xlabel("µm")
        plt.ylabel("µm")
        plt.suptitle(peakname + " " + quantities[i])
        plt.colorbar(im, label=cbar_labels[i])
        plt.savefig(folder+"/"+peakname + "_" + quantities[i]+".pdf", format="pdf", dpi=300)
        plt.show()

    return()

""" Apply a moving average window [size (2s+1)*(2s+1)] to a masked 2D array """
def moving_2D_average(data, s):
    # Extract size of the data array
    ny,nx = data.shape
    # Create empty masked array to store the averaged data, with the mask set
    # to False (0) everywhere
    data_av = np.ma.array(np.zeros((ny,nx)), mask = np.zeros((ny,nx)))
    
    # Scan over the array, avoiding the array edges so far as to always have s 
    # neighbours around each data point in which a gradient is computed
    for x in range(nx-2*s):
        x = x+s
        
        for y in range(ny-2*s):
            y = y+s
            fitfails = 0
            
            # Scan over the viscinity of the current data point, considering a
            # square with side length 2s+1.
            dividend = 0
            divisor = 0
            for i in range(-s, s+1): # s+1 makes sure that i=(-1,0,1) if s=1
                for j in range(-s, s+1):
                    
                    # If the current viscinity data point is masked, raise the
                    # counter of fit fails
                    if data.mask[y+j, x+i] == True:
                        fitfails += 1
                        continue
                    
                    # If the current viscinity data point is not masked,
                    # include its value into the calculation of the local 
                    # average.
                    else:
                        value = data[y+j, x+i]
                        # The viscinity data points are weighted with
                        # 1 / their distance to the current data point (the one 
                        # the average is calculated for). The current data
                        # point itself gets a weight of 1.
                        if i==0 and j==0:
                            weight = 1
                        else:
                            weight = 1/sqrt(i**2 + j**2)
                                                    
                        dividend += value*weight
                        divisor += weight
                    
            # If at the end there have been <= 3 fitfails, calculate average
            if fitfails <= 3:
                data_av[y,x] = dividend/divisor
            
            # Otherwise, mask the current data point
            else:
                data_av.mask[y,x] = True
                
    # Apply mask around the edge where no average value was calculated
    data_av.mask[0:s,:] = 1
    data_av.mask[-s:,:] = 1
    data_av.mask[:,0:s] = 1
    data_av.mask[:,-s:] = 1
                
    return(data_av)     

    

""" Do the twist angle related maps """
def map_theta(fitresults, fiterrors, pdict, folder):
    
    ############################### Part 1 ####################################
    # Create a function that converts the TA peak position to the twist angle #
    ###########################################################################
    
    # Load phonon dispersion
    dispersion = np.loadtxt("phondis_graphene.dat")
    
    # Extract the dispersion of the TA branch
    # Use only the K-to-Gamma data (137:237)
    # Invert list to get Gamma-K instead of K-Gamma
    start, stop = 137, 237
    ta_branch = dispersion[start:stop,2][::-1]
        
    # This is the distance Gamma-K in units of the graphene lattice constant a
    gamma, k = 0, 4*pi/3/a
    # Create an array of crystal momenta from Gamma to K
    crystal_momentum = np.linspace(gamma, k, stop-start)
    
    # Convert crystal momenta to corresponding twist angles using the geometric
    # expression given in http://dx.doi.org/10.1021/nl201370m
    theta_deg = np.degrees(2*arcsin(sqrt(3)*a/8/pi*crystal_momentum))

    # Use scipy interpolation to create a function of the twist angle as a 
    # function of the phonon frequency, i.e., the TA peak position
    TA_position_to_theta = interp1d(ta_branch, theta_deg, kind="cubic")
    
    
    ############################### Part 2 ####################################
    #          Using the function from part 1, map the twist angle            #
    ###########################################################################
    
    # Retrieve required variables from the peak dictionary
    nx, ny = pdict["size_px"] # Scan size in pixels
    sx, sy = pdict["size_um"] # Scan size in microns
    
    # Calculate the twist angle array
    theta = TA_position_to_theta(fitresults[:,:,1])
    
    # Mask the twist angle array in all points where the fit failed
    # Define the condition for masking
    mask_condition = fiterrors != 0
    # Create the mask, expand it to match the shape of fitresults
    mask = np.array(mask_condition, dtype=bool)
    # Apply the mask to fitresults
    theta_ma = np.ma.array(theta, mask=mask)
    
    # Map the twist angle
    im = plt.imshow(theta_ma, extent = [0, sx, 0, sy], cmap="gist_rainbow")
    plt.xlabel("µm")
    plt.ylabel("µm")
    plt.suptitle("Twist angle")
    plt.colorbar(im, label=r"$\vartheta_{TA}$ (°)")
    plt.savefig(folder + "/" + "twist_angle.pdf", format="pdf", dpi=300)
    plt.show()
    
    
    ############################### Part 2 ####################################
    #               Calculate and map the twist angle gradient                #
    ###########################################################################
    
    # Retrieve required variables from the peak dictionary
    max_gradient = pdict["max_gradient"] # °/µm, maximum gradient value in the map
    
    # Calculate the twist angle gradient in x and y direction
    gradient_x, gradient_y = np.gradient(theta_ma, edge_order=1)
    
    # Calculate the moving average of gradient in x and y direction
    gradient_x_av = moving_2D_average(gradient_x, s=1)
    gradient_y_av = moving_2D_average(gradient_y, s=1)
    
    # Calculate absolute value of the gradient
    grad_theta_ma = np.sqrt(gradient_x_av**2 + gradient_y_av**2)
    
    # Divide numerical gradient by step size to convert in °/µm
    # sx/nx is the step size of the scan in pixel
    grad_theta_ma = grad_theta_ma/(sx/nx)

    # Mask all gradient values that are above the threshold set by the user
    updated_mask = np.logical_or(grad_theta_ma.mask, grad_theta_ma.data > max_gradient)
    grad_theta_ma = np.ma.masked_array(grad_theta_ma.data, updated_mask)
    
    # Map the twist angle gradient
    im = plt.imshow(grad_theta_ma, extent = [0, sx, 0, sy], cmap="gist_rainbow")
    plt.xlabel("µm")
    plt.ylabel("µm")
    plt.suptitle("Twist angle gradient")
    plt.colorbar(im, label=r"|$\nabla\vartheta_{TA}$| (°/µm)")
    plt.savefig(folder + "/" + "twist_angle_gradient.pdf", format="pdf", dpi=300)
    plt.show()

    return()
    
    
    


