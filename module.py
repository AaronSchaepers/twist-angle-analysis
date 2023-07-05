#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Aaron
@date: 2023/06/19
@github: https://github.com/AaronSchaepers/twist-angle-analysis


IDEAS
    - Use the threshold dynamically in the interactive plot
    - Make one interactive map per peak 
    - Disable intensity threshold, at least upper bound, because the max int 
      values depend heavily on the integration time
    - Export pdfs with dpi = 800 directly from interactive plot
    - Progress bars

PREPARATION

    Export Raman maps from ProjectFive as .txt by choosing "Table" in the
    export menu. Set the spectral unit to "rel. 1/cm".
    Also, make sure the file "phondis_graphene.dat" is located in the same 
    directory as this code.


HOW TO USE THIS CODE

    As you scroll down, you will see that this code has four sections:
        1. User input section
        2. Advanced user input section
        3. Don't touch section A (functions)
        4. Don't touch section B (code execution)
        
    Obviously, the last two sections are not meant to be edited as they contain
    the main body of this code (kind of the back end).
    As a user, here is what you have to do to analyse your Raman Scan:
    
    1. User input section:
        - Provide the required information about the scan (see comments)
        - Choose which peaks to fit and to map
        - Run this code to start the fitting and mapping, or proceed to the

    2. Advanced user input section:
        - For each peak, set the starting parameters for the Lorentzian fit
        - For each fitting parameter and each peak, define a threshold interval
          within which the fit results are accepted as realistic.
        - Run this code to start the fitting and mapping procedure
        
    The default values give good results for the TA, G, LO and 2D peak. You
    might want to save them by commenting them out if you change them so they
    don't get los

    If mapping is activated, the code will do two things. Firstly, it will 
    export the maps of all selected peaks to the same directory where the scan 
    data is stored. Secondly, it will open up an interactive window ONLY FOR
    THE LAST SELECTED PEAK. That means you can run the fitting and mapping for
    multiple peaks and you will find the results in the directory given in
    Section 1. But only for the last peak selected for mapping, the interactive
    window will show up.
    
    The interactive window shows the maps of intensity, position and linewidth
    and also the Raman spectrum and fit in any point you double click on any 
    of the maps. It is very helpful when it comes to finding the right
    threshold parameters that serve to exclude faulty spectra.
"""

import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from numpy import pi, arcsin, sqrt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

a = 0.246 # Graphene lattice constant in nm

###############################################################################
# 3.1 Functions for importing and exporting data 
###############################################################################

""" Save an object using the pickle module """
def save_object(folder, file, name):
    with open(str(folder) + "/" + str(name)+".pkl", "wb") as f:
            pickle.dump(file, f, pickle.HIGHEST_PROTOCOL)
    return()


""" Read Raman scan data from .txt file, set baseline to zero """
def import_raman_scan(filepath, size_x_px, size_y_px, min_mean_range, max_mean_range):
    
    # A tuple of strings with 1600 entries, i.e., one entry per CCD pixel, i.e., one entry
    # per data point in the Raman spectrum.
    # Each of the 1600 strings is a succession of numbers giving the CCD counts
    # on the corresponding pixel for all the spatial data points measured in the scan.
    lines = tuple(open(filepath, 'r')) 
    
    # A list containing 1600 sublists corresponding to the 1600 CCD pixels.
    # Each sublist contains size_x_px*size_y_px strings giving the the CCD counts on that 
    # pixel for all the spatial data points.
    data_str = [x.split('\t') for x in lines]
    
    # A list containing size_x_px*size_y_px+1 sublists with 1600 entries each. The first list 
    # contains the wave numbers, the other lists are the measured spectra.
    data_float = np.array(data_str, dtype=float).T
    
    xlist = data_float[0].copy()  # Extract first row which contains the x-data
    data_float = np.delete(data_float, 0, axis=0)  # Remove the first row from data_float
    
    # A 3D array with the size size_x_px*size_y_px*1600, i.e., the first two axes correspond
    # to the x and y axes of the scan and the spectra will be stacked along the
    # third axis.
    data = np.zeros((size_y_px, size_x_px, len(data_str)))
    
    # Get indices of the averaging interval for baseline shifting
    i_xmin_mean = np.abs(xlist - min_mean_range).argmin()
    i_xmax_mean = np.abs(xlist - max_mean_range).argmin()
    
    # Scan over all spatial data points
    x,y = 0,0
    for i in range(size_x_px*size_y_px):
        # Calculate the average CCD counts in a range without peaks
        ymean = np.mean(data_float[i, i_xmin_mean:i_xmax_mean])
        # Subtract that average from the entire spectrum and stack the spectrum
        # into the data array
        data[y,x] = [y-ymean for y in data_float[i]]          

        x = x+1
        if x == size_x_px:
            y = y+1
            x = 0
            
    # Convert xlist to array for more efficient function evaluation
    xdata = np.asarray(xlist)
    
    return(xdata, data)


###############################################################################
# 3.2 Functions for the fitting procedure 
###############################################################################

""" Lorentzian as used in spectroscopy """
def lorentzian(x, b, c, x0, lw):
   f = b + c / pi * lw / ((x-x0)**2 + lw**2)
   return(f)
    

""" Jacobian matrix of the Lorentzian as defined above """
def jacobian(x, b, c, x0, lw):
    # Define partial derivatives
    Ldiffb = np.ones(x.shape)
    Ldiffc = lw / pi / ((x-x0)**2 + lw**2)
    Ldiffx0 = c / pi * 2 * lw * (x-x0) / ((x-x0)**2 + lw**2)**2
    Ldifflw = c / pi * (1 / ((x-x0)**2 + lw**2) -  2 * lw**2 / ((x-x0)**2 + lw**2)**2 )
    jac = np.stack((Ldiffb, Ldiffc, Ldiffx0, Ldifflw)).T
    return(jac)


""" Perform a Lorentzian fit across a map of spectra """
def fit_to_map(xdata, data, pdict):
    
    # Retrieve required variables from the dictionnary
    size_x_px = pdict["size_x_px"] # Size of the scan in pixels
    size_y_px = pdict["size_y_px"]
    fitrange = pdict["fitrange"] # Spectral range for the fit
    p0 = [pdict["init_off"], pdict["init_int"], pdict["init_pos"], pdict["init_lw"]] # Initial values
    
    # Array with (b, c, x0, lw) for each x/y-datapoint in which fit was successful
    fitresults = np.zeros((size_y_px, size_x_px, 4))
    
    # This expression searches the xlist for the wave number values 
    # closest to the chosen fitrange and retrieves their indexes
    i_start = np.abs(xdata - fitrange[0]).argmin()
    i_stop = np.abs(xdata - fitrange[1]).argmin()
    
    # Slice xdata to keep only the relevant peak
    xdata_fit = xdata[i_start:i_stop]
    
    # Define boundaries to speed up the fitting routine
    #                   b        c       x0      lw
    lbounds = np.array((-np.inf, 0,      0,      0))
    ubounds = np.array((np.inf,  np.inf, np.inf, np.inf))
    
    print("Fitting progress:")
    for x in range(size_x_px):
        progress = np.round(x/size_x_px*100, 1)
        print(f"{progress}%")
        for y in range(size_y_px):
            
            # Retrieve the spectrum of the current scan data point
            ydata = data[y,x]
            
            # Slice ydata to keep only the relevant peak
            ydata_fit = ydata[i_start:i_stop]
            
            # Determine highest data point, use its wavenumber as initial value
            # for the position
            p0[2] = xdata_fit[np.argmax(ydata_fit)]
            
            # Try to perform the fit
            try:
                popt, pcov = curve_fit(lorentzian, xdata_fit, ydata_fit, p0, bounds=(lbounds, ubounds), jac=jacobian)
                popt_std = np.sqrt(np.diag(pcov))
                
                # Check if any fit results contain NaN values
                if np.isnan(popt).any() or np.isnan(popt_std).any():
                    continue
                # If not: write best fit parameters into the fitresults
                else:
                    fitresults[y,x] = popt
                    
            # If the fit fails, don't change the fitresults entry in that point
            # so that it remains [0,0,0,0]               
            except: 
                continue
                
    return(fitresults)



###############################################################################
# 3.3 Functions for the mapping 
###############################################################################

""" Map the intensity, position and linewidth of a given peak """
def map_lorentz_parameters(fitresults, fiterrors, pdict, folder):
    # Retrieve required variables from the dictionnary
    nx, ny = pdict["size_px"] # Scan size in pixels
    sx, sy = pdict["size_um"] # Scan size in microns
    peakname = pdict["peakname"] # Name of the peak
    (thresh_c, thresh_x0, thresh_lw) = pdict["params_thresh"]    
    
    # Conditions to check if the best fit parameters fall into their threshold interval
    conditions = [
        (thresh_c[0] > fitresults[:, :, 1]),
        (fitresults[:, :, 1] > thresh_c[1]),
        (thresh_x0[0] > fitresults[:, :, 2]),
        (fitresults[:, :, 2] > thresh_x0[1]),
        (thresh_lw[0] > fitresults[:, :, 3]),
        (fitresults[:, :, 3] > thresh_lw[1]),
        (fiterrors != 0)
    ]
    
    # Create the mask in 2D
    mask = np.logical_or.reduce(conditions)
    
    # Apply the mask to fitresults
    fitresults_ma = np.ma.masked_array(fitresults, mask=np.broadcast_to(mask[:, :, np.newaxis], fitresults.shape))

    # The mapped quantities
    quantities = ["intensity", "position", "linewidth"]
    # Their units
    cbar_labels = ["Intensity (arb. u.)", r"$\omega$ (1/cm)", r"$\Gamma$ (1/cm)"]
    
    # In a loop, plot peak intensity, position and linewidth
    for i in range(3):
        # The first plt.close() avoids funny interactions with the interactive plot
        plt.close()
        # In this line, i = i+1 to skip the first parameter which is the offset
        im = plt.imshow(fitresults_ma[:,:,i+1], extent = [0, sx, 0, sy], cmap="gist_rainbow")
        plt.xlabel("µm")
        plt.ylabel("µm")
        plt.suptitle(peakname + " " + quantities[i])
        plt.colorbar(im, label=cbar_labels[i])
        plt.savefig(folder+"/"+peakname + "_" + quantities[i]+".pdf", format="pdf", dpi=300)
        plt.close()

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
    (thresh_c, thresh_x0, thresh_lw) = pdict["params_thresh"]    
    
    # Conditions to check if the best fit parameters fall into their threshold interval
    conditions = [
        (thresh_c[0] > fitresults[:, :, 1]),
        (fitresults[:, :, 1] > thresh_c[1]),
        (thresh_x0[0] > fitresults[:, :, 2]),
        (fitresults[:, :, 2] > thresh_x0[1]),
        (thresh_lw[0] > fitresults[:, :, 3]),
        (fitresults[:, :, 3] > thresh_lw[1]),
        (fiterrors != 0)
    ]
    
    # Create the mask in 2D
    mask = np.logical_or.reduce(conditions)
    
    # Apply mask to fitresults and fill masked entries with 0 because interpl1d
    # does not accept masked arrays as an input
    fitresults_0 = np.ma.masked_array(fitresults, mask[..., np.newaxis], fill_value=0)
    
    # Calculate the twist angle array
    theta = TA_position_to_theta(fitresults_0[:,:,2])
    
    # Apply the mask to theta
    theta_ma = np.ma.masked_array(theta, mask=mask)
    
    # Map the twist angle
    im = plt.imshow(theta_ma, extent = [0, sx, 0, sy], cmap="gist_rainbow")
    plt.xlabel("µm")
    plt.ylabel("µm")
    plt.suptitle("Twist angle")
    plt.colorbar(im, label=r"$\vartheta_{TA}$ (°)")
    plt.savefig(folder + "/" + "twist_angle.pdf", format="pdf", dpi=300)
    plt.close()
    
    
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
    plt.close()

    return()