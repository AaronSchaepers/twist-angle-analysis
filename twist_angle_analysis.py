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
""" 1. User input """
###############################################################################

# Directory of the Raman data
folder = "/Users/Aaron/Desktop/ETIRF04-100/G+2D+LO" 

# Name of the .txt file containing the Raman data, given with suffix
file = "ETIRF04-100_G_LO_2D.txt" 

# In 1/cm, a spectral range without any features. This is used to calculate 
# the mean background noise which is subtracted from the data
spectral_mean_range = (2000, 2100) 

size_px = (50, 40)    # Size of the Scan in pixels
size_um = (87, 59)        # Size of the Scan in µm

# What peaks shall be fitted?
b_fit_TA = False
b_fit_G = False
b_fit_LO = False
b_fit_2D = False

# What peaks shall be mapped?
b_map_TA = False
b_map_G = True
b_map_LO = True
b_map_2D = True



###############################################################################
""" 2. Advanced user input """
###############################################################################

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

thresh_TA_c = [20, 6000]      # Intensity
thresh_TA_x0 = [250, 275]     # Position
thresh_TA_lw = [0.2, 9]       # Linewidth

thresh_G_c = [1E3, 1E15]      # Intensity
thresh_G_x0 = [1577, 1587]    # Position
thresh_G_lw = [4, 25]         # Linewidth

thresh_LO_c = [300, 1E5]      # Intensity
thresh_LO_x0 = [1610, 1630]   # Position
thresh_LO_lw = [0.4, 15]      # Linewidth

thresh_2D_c = [1E3, 1E7]      # Intensity
thresh_2D_x0 = [2650, 2700]   # Position
thresh_2D_lw = [8, 30]        # Linewidth

max_gradient = 1 # °/µm, upper bound for the twist angle gradient map. Larger
                 #       values will be excluded from the map.


###############################################################################
""" 3. Don't touch section A (functions) """
###############################################################################


###############################################################################
# 3.1 Functions for importing and exporting data 
###############################################################################

""" Save an object using the pickle module """
def save_object(folder, file, name):
    with open(str(folder) + "/" + str(name)+".pkl", "wb") as f:
            pickle.dump(file, f, pickle.HIGHEST_PROTOCOL)
    return()


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
    data_float = np.array(data_str, dtype=float).T
    
    xlist = data_float[0].copy()  # Extract first row which contains the x-data
    data_float = np.delete(data_float, 0, axis=0)  # Remove the first row from data_float
    
    # A 3D array with the size nx*ny*1600, i.e., the first two axes correspond
    # to the x and y axes of the scan and the spectra will be stacked along the
    # third axis.
    data = np.zeros((ny,nx,len(data_str)))
    
    # Get indices of the averaging interval for baseline shifting
    i_xmin_mean = np.abs(xlist - spectral_mean_range[0]).argmin()
    i_xmax_mean = np.abs(xlist - spectral_mean_range[1]).argmin()
    
    # Scan over all spatial data points
    x,y = 0,0
    for i in range(nx*ny):
        # Calculate the average CCD counts in a range without peaks
        ymean = np.mean(data_float[i, i_xmin_mean:i_xmax_mean])
        # Subtract that average from the entire spectrum and stack the spectrum
        # into the data array
        data[y,x] = [y-ymean for y in data_float[i]]          

        x = x+1
        if x == nx:
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
    nx, ny = pdict["size_px"] # Scan size in pixels
    fitrange = pdict["fitrange"] # Spectral range for the fit
    p0 = pdict["startparams"] # Starting values
    
    # Create arrays to store the fit results
    fitresults = np.zeros((ny,nx,4))          # Array with (b, c, x0, lw) for each x/y-datapoint in which fit was successful
    fitresults_std = np.zeros((ny,nx,4))      # Array with the uncertainties of the fit results
    fiterrors = np.zeros((ny,nx), dtype= "object")    # Array with inormation if and why the fit failed for each x/y-datapoint
    
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
    for x in range(nx):
        progress = np.round(x/nx*100, 1)
        print(f"{progress}%")
        for y in range(ny):
            
            # Retrieve the spectrum of the current scan data point
            ydata = data[y,x]
            
            # Slice ydata to keep only the relevant peak
            ydata_fit = ydata[i_start:i_stop]
            
            # Determine highest data point, use position as x0 starting value
            p0[2] = xdata_fit[np.argmax(ydata_fit)]
            
            # Try to perform the fit
            try:
                popt, pcov = curve_fit(lorentzian, xdata_fit, ydata_fit, p0, bounds=(lbounds, ubounds), jac=jacobian)
                popt_std = np.sqrt(np.diag(pcov))
                
                # Check if any fit results contain NaN values
                if np.isnan(popt).any() or np.isnan(popt_std).any():
                    fiterrors[y, x] = "NaN"
                    continue
                # If not: write best fit parameters into the fitresults
                else:
                    fitresults[y,x] = popt
                    fitresults_std[y,x] = popt_std
                    
            # If the fit fails, put a note in the fitresults array and proceed.                    
            except: 
                fiterrors[y, x] = "fit fail"
                
    return(fitresults, fitresults_std, fiterrors)



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


###############################################################################
# 3.4 Experimental: interactive stuff 
###############################################################################

# This function is activated when double clicking in tne interactive map and 
# plots the spectrum + fit in the point that was clicked
def onclick(event, fitresults, fiterrors, pdict):
    if event.dblclick:
        if event.button == 1:
            if isinstance(event.ydata, float): # Check if the doubleclick event happened inside the plot and returns correct values
                # Exact location of click in axis units
                valx, valy = event.xdata, event.ydata
                
                # Find index coordinates of click position
                x_map = np.abs(xaxis - valx).argmin()
                y_map = np.abs(yaxis - valy).argmin()
                
                # Extract required variables from the pdict
                plotrange = pdict["plotrange"]
                fitrange = pdict["fitrange"]
                (thresh_c, thresh_x0, thresh_lw) = pdict["params_thresh"]    
                
                # Find start and stop indices of plotrange and fitrange in spectrum
                i_plotstart = np.abs(xdata - plotrange[0]).argmin()
                i_plotstop = np.abs(xdata - plotrange[1]).argmin()
                i_fitstart = np.abs(xdata - fitrange[0]).argmin()
                i_fitstop = np.abs(xdata - fitrange[1]).argmin()

                # Create x and y array in the plotrange
                xdata_plot = xdata[i_plotstart:i_plotstop]
                ydata_plot = data[-y_map, x_map, i_plotstart:i_plotstop]
                # Create x array in the fitrange
                xdata_fit = xdata[i_fitstart:i_fitstop]

                # Clear the previous spectrum + fit from the plot
                ax3.cla()
                
                # Scatter the spectrum
                ax3.scatter(xdata_plot, ydata_plot, s=5, zorder=1)
                
                # Plot the fit if it succeeded in that point
                if fiterrors[-y_map, x_map] == 0:
                    ax3.plot(xdata_fit, lorentzian(xdata_fit, *fitresults[-y_map, x_map]), color="tab:orange")
                    ax3.set_xlim(plotrange)
                    
                    # Assemble text string with best fit parameters
                    text_I = r"I = %.0f " %(fitresults[-y_map, x_map][1])
                    text_x0 = "\n" + "$\omega$ = %.2f rel. 1/cm" %(fitresults[-y_map, x_map][2])
                    text_lw = "\n" + "$\Gamma$ = %.2f 1/cm" %(fitresults[-y_map, x_map][3])
                    
                    # Define conditions to check if the fitresults violate threshold boundaries
                    cond_c = (thresh_c[0] > fitresults[-y_map, x_map, 1]) or (fitresults[-y_map, x_map, 1] > thresh_c[1])
                    cond_x0 = (thresh_x0[0] > fitresults[-y_map, x_map, 2]) or (fitresults[-y_map, x_map, 2] > thresh_x0[1])
                    cond_lw = (thresh_lw[0] > fitresults[-y_map, x_map, 3])  or (fitresults[-y_map, x_map, 3] > thresh_lw[1])

                    # Add information about threshold violations to the string
                    if cond_c:
                        ax3.text(0.01, 0.95, text_I, ha='left', va='top', transform=ax3.transAxes, color="tab:red")
                    else:
                        ax3.text(0.01, 0.95, text_I, ha='left', va='top', transform=ax3.transAxes)
                    
                    if cond_x0:
                        ax3.text(0.01, 0.94, text_x0, ha='left', va='top', transform=ax3.transAxes, color="tab:red")
                    else:
                        ax3.text(0.01, 0.94, text_x0, ha='left', va='top', transform=ax3.transAxes)
                        
                    if cond_lw:
                        ax3.text(0.01, 0.85, text_lw, ha='left', va='top', transform=ax3.transAxes, color="tab:red")
                    else:
                        ax3.text(0.01, 0.85, text_lw, ha='left', va='top', transform=ax3.transAxes)
                        
                    # Plot the text
                    #ax3.text(0.01, 0.95, textstr, ha='left', va='top', transform=ax3.transAxes)

                # Update the plot in the window
                ax3.set_xlabel("Raman shift (rel. 1/cm)")
                ax3.set_ylabel("CCD counts")
                ax3.figure.canvas.draw()
                
            else:
                print('Please double click inside a map!')
                return None


def make_figure(fitresults, fiterrors, pdict):
    
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
    
    # Build the figure with grid and subplots
    fig = plt.figure(figsize = (15, 10))
    grid = mpl.gridspec.GridSpec(6, 6, figure = fig, hspace=0.6, wspace=0.8)
    ax0 = fig.add_subplot(grid[:3, 0:2])
    ax1 = fig.add_subplot(grid[:3, 2:4])
    ax2 = fig.add_subplot(grid[:3, 4:6])
    ax3 = fig.add_subplot(grid[3:6, :])
    
    # Plot the maps
    im0 = ax0.imshow(fitresults_ma[:,:,1], extent = [0, sx, 0, sy], cmap="gist_rainbow") # Intensity
    im1 = ax1.imshow(fitresults_ma[:,:,2], extent = [0, sx, 0, sy], cmap="gist_rainbow") # Position
    im2 = ax2.imshow(fitresults_ma[:,:,3], extent = [0, sx, 0, sy], cmap="gist_rainbow") # Linewidth
    
    # Label axes
    ax0.set_xlabel("µm")
    ax0.set_ylabel("µm")
    ax1.set_xlabel("µm")
    ax1.set_ylabel("µm")
    ax2.set_xlabel("µm")
    ax2.set_ylabel("µm")
    ax3.set_xlabel("Raman shift (rel. 1/cm)")
    ax3.set_ylabel("CCD counts")
    
    # Add colorbars
    plt.colorbar(im0, ax=ax0, label="Intensity (arb. u.)")   
    plt.colorbar(im1, ax=ax1, label="Position (rel. 1/cm)") 
    plt.colorbar(im2, ax=ax2, label="Linewidth (1/cm)") 
    
    # Set titles
    ax0.set_title(peakname + " intensity")
    ax1.set_title(peakname + " position")
    ax2.set_title(peakname + " linewidth")
   
    fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, fitresults, fiterrors, pdict))
    
    return(ax0, ax1, ax2, ax3)


###############################################################################
""" 4. Don't touch section B (executing the code) """
###############################################################################


###############################################################################
# 4.1 Assemble peak dictionaries with information for mapping and fitting
###############################################################################

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
dict_TA["plotrange"] = (220, 310)            # Data range for plotting
dict_TA["params_thresh"] = (thresh_TA_c, thresh_TA_x0, thresh_TA_lw)
dict_TA["max_gradient"] = max_gradient
save_object(folder, dict_TA, "dict_TA")

# The first two lines check if the peak dictionary already exists and create
# it only if it doesn't. That way, stuff added later on won't be overwritten.
try: type(dict_G)
except NameError: dict_G = {}
dict_G["startparams"] = [startparams_G[0], startparams_G[1], 0, startparams_G[2]]
dict_G["size_px"] = size_px
dict_G["size_um"] = size_um
dict_G["peakname"] = "G"
dict_G["fitrange"] = (1500, 1610)             # Data range for fitting
dict_G["plotrange"] = (1500, 1800)            # Data range for plotting
dict_G["params_thresh"] = (thresh_G_c, thresh_G_x0, thresh_G_lw)
save_object(folder, dict_G, "dict_G")

# The first two lines check if the peak dictionary already exists and create
# it only if it doesn't. That way, stuff added later on won't be overwritten.
try: type(dict_LO)
except NameError: dict_LO = {}
dict_LO["startparams"] = [startparams_LO[0], startparams_LO[1], 0, startparams_LO[2]]
dict_LO["size_px"] = size_px
dict_LO["size_um"] = size_um
dict_LO["peakname"] = "LO"
dict_LO["fitrange"] = (1610, 1800)             # Data range for fitting
dict_LO["plotrange"] = (1500, 1800)            # Data range for plotting
dict_LO["params_thresh"] = (thresh_LO_c, thresh_LO_x0, thresh_LO_lw)
save_object(folder, dict_LO, "dict_LO")

# The first two lines check if the peak dictionary already exists and create
# it only if it doesn't. That way, stuff added later on won't be overwritten.
try: type(dict_2D)
except NameError: dict_2D = {}
dict_2D["startparams"] = [startparams_2D[0], startparams_2D[1], 0, startparams_2D[2]]
dict_2D["size_px"] = size_px
dict_2D["size_um"] = size_um
dict_2D["peakname"] = "2D"
dict_2D["fitrange"] = (2600, 2800)             # Data range for fitting
dict_2D["plotrange"] = (2400, 2900)            # Data range for plotting
dict_2D["params_thresh"] = (thresh_2D_c, thresh_2D_x0, thresh_2D_lw)
save_object(folder, dict_2D, "dict_2D")


###############################################################################
# 4.3 Read the Raman scan data
###############################################################################

xdata, data = read_raman_scan(folder, file, size_px, spectral_mean_range)

###############################################################################
# 4.4 Perform the actual fitting and mapping
###############################################################################

# Create lists that are needed to locate the index coordinates of a click event
xaxis = np.linspace(0, size_um[0], size_px[0])
yaxis = np.linspace(0, size_um[1], size_px[1]) 

# Fitting
if b_fit_TA == True:
    fitresults_TA, fitresults_std_TA, fiterrors_TA = fit_to_map(xdata, data, dict_TA)
if b_fit_G == True:
    fitresults_G, fitresults_std_G, fiterrors_G = fit_to_map(xdata, data, dict_G)
if b_fit_LO == True:
    fitresults_LO, fitresults_std_LO, fiterrors_LO = fit_to_map(xdata, data, dict_LO)
if b_fit_2D == True:
    fitresults_2D, fitresults_std_2D, fiterrors_2D = fit_to_map(xdata, data, dict_2D)
     
# Save maps, open interactive figures
if b_map_TA == True:
    map_lorentz_parameters(fitresults_TA, fiterrors_TA, dict_TA, folder)
    map_theta(fitresults_TA, fiterrors_TA, dict_TA, folder)
    ax0, ax1, ax2, ax3 = make_figure(fitresults_TA, fiterrors_TA, dict_TA)
if b_map_G == True:
    map_lorentz_parameters(fitresults_G, fiterrors_G, dict_G, folder)
    ax0, ax1, ax2, ax3 = make_figure(fitresults_G, fiterrors_G, dict_G)
if b_map_LO == True:
    map_lorentz_parameters(fitresults_LO, fiterrors_LO, dict_LO, folder)
    ax0, ax1, ax2, ax3 = make_figure(fitresults_LO, fiterrors_LO, dict_LO)
if b_map_2D == True:
    map_lorentz_parameters(fitresults_2D, fiterrors_2D, dict_2D, folder)
    ax0, ax1, ax2, ax3 = make_figure(fitresults_2D, fiterrors_2D, dict_2D)
    
    