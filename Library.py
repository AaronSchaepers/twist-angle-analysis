#!/usr/bin/env python3
# -*- coding: utf-8 -*
"""
Created on Fri Oct 11 16:56:36 2019
@author: Aaron
"""

import numpy as np
from numpy import sqrt, arcsin, pi
import numpy.ma as ma

from scipy.interpolate import interp1d
import math
from math import floor
import lmfit as lm
import pickle

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)
mpl.rcParams.update({'font.size': 12})

mpl.rc('axes', linewidth=1)
mpl.rcParams["ytick.major.pad"] = 5
mpl.rcParams["xtick.major.pad"] = 5
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.direction"] = "in"
ticklength = 5
tickwidth = 1



###############################################################################
"""Read measurement data from table, set baseline to zero"""
###############################################################################

def read_raman_data(filename, sample_size, mean_range, folder=0):
    nx, ny = sample_size
    if folder==0:
        lines = tuple(open(filename + ".txt", 'r'))
    else:
        lines = tuple(open(folder + "/"+ filename + ".txt", 'r'))
    data_str = [x.split('\t') for x in lines]
    data_float = [[float(liste[i]) for liste in data_str ] for i in range(len(data_str[0]))]
    #extract first list which contains the x-data
    xlist = data_float.pop(0)

    # Required if the "base line" of the data should be set to zero
    xmin_mean = min(range(len(xlist)), key=lambda i: abs(xlist[i]-mean_range[0]))
    xmax_mean = min(range(len(xlist)), key=lambda i: abs(xlist[i]-mean_range[1]))

    #sort y-data into 3D-array
    data = np.zeros((nx,ny,len(data_str)))
    
    x,y = 0,0
    for i in range(nx*ny):
        # Shift baseline to zero
        ymean = np.mean(data_float[i][xmin_mean:xmax_mean])        
        data[x,y] = [y-ymean for y in data_float[i]]
        
        #data[x,y] = data_float[i]
        x = x+1

        if x == nx:
            y = y+1
            x = 0

    return(xlist, data, lines)


###############################################################################
"""Read Raman data of signle spectrum"""
###############################################################################

def read_single_spectrum(filename, folder=0):
    if folder==0:
        lines = tuple(open(filename, 'r'))
    else:
        lines = tuple(open(folder + "/"+ filename, 'r'))
        
    lines = tuple(open(folder + "/"+ filename, 'r'))
    data_str = [x.split('\t') for x in lines]
    data_float = [[float(liste[i]) for liste in data_str ] for i in range(len(data_str[0]))]
    
    # Required if the "base line" of the data should be set to zero
    # =============================================================================
    # xmin_mean = min(range(len(xlist)), key=lambda i: abs(xlist[i]-mean_range[0]))
    # xmax_mean = min(range(len(xlist)), key=lambda i: abs(xlist[i]-mean_range[1]))
    # =============================================================================

    # Extract first list which contains the x-data
    xlist = data_float.pop(0)
    ylist = data_float[0]

    return(xlist, ylist)


###############################################################################
"""Plot examplary spectrum"""
###############################################################################

def plot_spectrum(xlist, ylist, plotsize, save=False, title=0, folder=0):
# =============================================================================
#     # Find indices of plotrange values
#     start = min(range(len(xlist)), key=lambda i: abs(xlist[i]-plotrange[0]))
#     stop = min(range(len(xlist)), key=lambda i: abs(xlist[i]-plotrange[1]))
# =============================================================================

    fig = plt.figure(figsize=plotsize)
    ax1 = fig.add_subplot(211)
    ax1.plot(xlist, ylist, zorder=1)
    ax1.set_xlabel(r"Raman shift (cm$^{-1})$")
    ax1.set_ylabel("Intensity (a.u.)")
    plt.tight_layout()
    if save == True:
        if folder == 0:
            plt.savefig(str(title)+".png", format="png", dpi=900)
        else:
            plt.savefig(str(folder)+"/"+str(title)+".png", format="png", dpi=900)           

        
###############################################################################
"""Histogram of array"""
###############################################################################
def histogram(data, errordata, bins, save=False, title=0, folder=0):
    # Mask the arrays to exclude areas where fit failed
    masked_data = np.ma.masked_where(errordata != 0, data)
    # Flatten data to 1D
    masked_data = np.ndarray.flatten(masked_data)
    fig = plt.figure(figsize=(4,4))
    ax1 = fig.add_subplot(111)
    ax1.hist(masked_data, bins, color="blue")
    if save == True:
        plt.suptitle(str(title))
        if folder == 0:
            plt.savefig(str(title)+".png", format="png", dpi=900)
        else:
            plt.savefig(str(folder)+"/"+str(title)+".png", format="png", dpi=900)           

    

###############################################################################
"""Save and load objects using the pickle module"""
###############################################################################

def save_obj(name, file, folder=0):
    if folder == 0:
        with open(str(name)+".pkl", "wb") as f:
            pickle.dump(file, f, pickle.HIGHEST_PROTOCOL)
    else:        
        with open(str(folder) + "/" + str(name)+".pkl", "wb") as f:
            pickle.dump(file, f, pickle.HIGHEST_PROTOCOL)
            
            
def load_obj(name, folder=0,fullpath=False):
#The kwag folder is required if the file is located in another folder
#than the script from which the function is called

    if folder == 0 and not fullpath:
        with open(str(name)+".pkl", "rb") as f:
            return(pickle.load(f))
    elif fullpath:
        
        with open(str(name), "rb") as f:
            return(pickle.load(f))
        
    else:
        with open(str(folder)+"/"+str(name)+".pkl", "rb") as f:
            return(pickle.load(f))


############################################################################### 
"""Mathematical functions"""
###############################################################################

# =============================================================================
# # Relativistic Lorentzian (difference of squares)
# def lorentzian(x, c, x0, lw):
#    f = c/((x**2-x0**2)**2 + (x0**2)*(lw**2))
#    return(f)
# =============================================================================

# Lorentzian as used in spectroscopy
def lorentzian(x, c, x0, lw):
   f = c / pi * lw / ((x-x0)**2 + lw**2)
   return(f)

def constant(x, b):
    return(b)

def linear(x,m,b):
    return(m*x+b)

# Carozo model: Twist angle theta as function of Q_intra
a = 0.246 # Graphene lattice ceonstant in nm
def theta(q):
    theta_rad = 2*arcsin(sqrt(3)*a/8/pi*q)
    #print(theta_rad)
    theta_deg = math.degrees(theta_rad)
    #print(theta_deg)
    return(theta_deg)   

# Linear dependency Theta(Pos(TA)), holds for up to 8.5 deg. Values from linear fit,
# see Deviation Carozo vs Linear for TA.py in //serveri2a/Transfer/Aaron/Code/Carozo Model
def theta_of_posTA(pos):
    m = 0.027247566645635615
    b = -0.15295066791136236
    theta = m*pos + b
    return(theta)


###############################################################################
"""Perform Lorentzian fit with linear background on a given spectrum"""
###############################################################################

def fit1L_linear_bg(xdata, ydata, p0, printresults=False, plot=False, returnPlotValues=False):    

    model = lm.Model(linear) + lm.Model(lorentzian)   # Make composite model by adding components
    b, c, x0, lw = p0                                   # Read out start values
    params = model.make_params()                        # Make parameters for the full model
    params['m'] = lm.parameter.Parameter(name="m", value=0)
    params['b'] = lm.parameter.Parameter(name="b", value=b)
    params['c'] = lm.parameter.Parameter(name="c", value=c, min=0)
    params['x0'] = lm.parameter.Parameter(name="x0", value=x0, min=0)
    params['lw'] = lm.parameter.Parameter(name="lw", value=lw, min=0)
    
    # Optional: Replace x0 initial value by position of max value in the list
    x0_new = xdata[np.argmax(ydata)]
    params['x0'] = lm.parameter.Parameter(name="x0", value=x0_new, min=0)
   
    
    # Do a fit over the selected data range
    result = model.fit(ydata, params, x=xdata)  # This object has the type of a class, tricky to work with
    final = result.best_fit                     # This object is an array, easy to work with

    # The lmfit routine returns fit results in a dictionary, write them into a list for easier handling:
    popt = []
    popt_std = []
    for key in result.params:
        if key != "m":
            popt.append(abs(result.params[key].value))
            popt_std.append(result.params[key].stderr)
    
    #print(result.fit_report()) #debug

    # Print the results if required
    if printresults == True:
        print(result.fit_report())

    # Plot the results if required
    if plot == True:
        fig = plt.figure(figsize=(8,8))
        ax1 = fig.add_subplot(211)
        ax1.plot(xdata, final, color="red", zorder=1)
        ax1.scatter(xdata, ydata, s=20, zorder=0)
        ax1.set_xlabel(r"Raman shift (cm$^{-1})$")
        ax1.set_ylabel("Intensity (a.u.)")
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if returnPlotValues:
        return final
    else:
        return(popt, popt_std)


###############################################################################
"""Perform Lorentzian fit on a given spectrum"""
###############################################################################

def fit1L(xdata, ydata, p0, printresults=False, plot=False, returnPlotValues=False):    

    model = lm.Model(constant) + lm.Model(lorentzian)   # Make composite model by adding components
    b, c, x0, lw = p0                                   # Read out start values
    params = model.make_params()                        # Make parameters for the full model
    params['b'] = lm.parameter.Parameter(name="b", value=b)
    params['c'] = lm.parameter.Parameter(name="c", value=c, min=0)
    params['x0'] = lm.parameter.Parameter(name="x0", value=x0, min=0)
    params['lw'] = lm.parameter.Parameter(name="lw", value=lw, min=0)
    
    # Optional: Replace x0 initial value by position of max value in the list
    x0_new = xdata[np.argmax(ydata)]
    params['x0'] = lm.parameter.Parameter(name="x0", value=x0_new, min=0)
    
    # Do a fit over the selected data range
    result = model.fit(ydata, params, x=xdata)  # This object has the type of a class, tricky to work with
    final = result.best_fit                     # This object is an array, easy to work with

    # The lmfit routine returns fit results in a dictionary, write them into a list for easier handling:
    popt = []
    popt_std = []
    for key in result.params:
        #popt.append(abs(result.params[key].value))
        popt.append(result.params[key].value)
        popt_std.append(result.params[key].stderr)
    
    #print(result.fit_report()) #debug

    # Print the results if required
    if printresults == True:
        print(result.fit_report())

    # Plot the results if required
    if plot == True:
        xlist_fit = np.arange(xdata[0], xdata[-1], 0.1)
        ylist_fit = [popt[0] + lorentzian(x, popt[1], popt[2], popt[3]) for x in xlist_fit]
        fig = plt.figure(figsize=(8,8))
        ax1 = fig.add_subplot(211)
        ax1.plot(xlist_fit, ylist_fit, color="red", zorder=1)
        ax1.scatter(xdata, ydata, s=20, zorder=0)
        ax1.set_xlabel(r"Rel. Raman shift (cm$^{-1})$")
        ax1.set_ylabel("Intensity (a.u.)")
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if returnPlotValues:
        return final
    else:
        return(popt, popt_std)


###############################################################################
"""Perform Lorentzian fit to every x/y-point of the sample"""
###############################################################################

def perform_fit1L(xlist, data, sample_size, pdict, masked=False, mask=0):
    
    nx, ny = sample_size
    background = pdict["background"]
    
    # Extract variables, starting values and thresholds from dictionary
    fitrange = pdict["fitrange"]
    p0 = pdict["startparams"] # Starting values
    (b_thresh, c_thresh, x0_thresh, lw_thresh, additional_thresh) = pdict["params_thresh"] # Threshold values
    b_lb, b_ub = b_thresh
    c_lb, c_ub = c_thresh # lb/ub = lower/upper boundary
    x0_lb, x0_ub = x0_thresh
    lw_lb, lw_ub = lw_thresh
    additional_threshold_lb, additional_threshold_ub = additional_thresh
    
    # Initiate variables to save overall extreme values of the fits
    max_x0, max_lw, max_c = 1,1,1
    min_x0, min_lw, min_c = 1E20, 1E20, 1E20

    fitresults = np.zeros((nx,ny,len(p0)))          # Array with (b, c, x0, lw) for each x/y-datapoint in which fit was successful
    fitresults_std = np.zeros((nx,ny,len(p0)))      # Array with the uncertainties of the fit results
    fiterrors = np.zeros((nx,ny), dtype= "object")  # Array with inormation if and why the fit failed for each x/y-datapoint

    # Scan through the data set and perform fit in every datapoint
    for x in range(nx):
        for y in range(ny):
            #print(x,y) # debug
            
            # In case a mask is applied: skip fit if it failed in this data point before
            if masked==True and mask[x,y]!=0:
                fiterrors[x,y] = "1st fit failed"
                continue

            ylist = data[x,y]

            # This expression searches the xlist for the values closest to the chosen fit_range and retrieves their indexes    
            start = min(range(len(xlist)), key=lambda i: abs(xlist[i]-fitrange[0]))
            stop = min(range(len(xlist)), key=lambda i: abs(xlist[i]-fitrange[1]))

            # Slice data to keep only the relevant peak
            ylist_fit = ylist[start:stop]
            xlist_fit = xlist[start:stop]
            
            # Do the fitting
            if background == "linear":
                popt, popt_std = fit1L_linear_bg(xlist_fit, ylist_fit, p0)
                b, c, x0, lw = popt
                b_std, c_std, x0_std, lw_std = popt_std
                try: additional_value = c_std/c
                except: additional_value = 1E20
                #print(lw_std)
            else:
                popt, popt_std = fit1L(xlist_fit, ylist_fit, p0)
                b, c, x0, lw = popt
                b_std, c_std, x0_std, lw_std = popt_std
                try: additional_value = lw_std
                except: additional_value = 1E20
           
            # Check the results do not contain infinities
            if None in popt:
                fiterrors[x,y] = "None val"
                continue
            if None in popt_std:
                fiterrors[x,y] = "None std"
                continue
            
            # Check if fit results are within the respective thresholds            
            if b_lb < b < b_ub:
                if c_lb < c < c_ub:    
                    if x0_lb < x0 < x0_ub:
                        if lw_lb < lw < lw_ub:
                            if additional_threshold_lb < additional_value < additional_threshold_ub:
                                    
                                    # If we got till here it worked, then save parameters
                                    fitresults[x,y] = [b, c, x0, lw]
                                    fitresults_std[x,y] = [b_std, c_std, x0_std, lw_std]
                                    # Check if we hit any min or max value in the parameters
                                    if x0 < min_x0: min_x0 = x0
                                    if lw < min_lw: min_lw = lw
                                    if c < min_c: min_c = c
                                    if x0 > max_x0: max_x0 = x0
                                    if lw > max_lw: max_lw = lw
                                    if c > max_c: max_c = c

                            else: fiterrors[x,y] = "ad. cond."  # Additional condition outside threshold
                        else: fiterrors[x,y] = "lw"             # Linewidth outside threshold
                    else: fiterrors[x,y] = "x0"                 # Position outside threshold
                else: fiterrors[x,y] = "c"                      # Intensity outside threshold 
            else: fiterrors[x,y] = "b"                          # Offset of the linear background outside threshold

    

    # Mask fitresults where the fit failed
    for i in range(len(fitresults[0,0])):
        fitresults[:,:,i] = np.ma.masked_where(fiterrors != 0, fitresults[:,:,i])

    # Save the extremal values of this peak on this sample
    maxvalues = (max_c, max_x0, max_lw)
    minvalues = (min_c, min_x0, min_lw)
    pdict["minvalues"] = minvalues
    pdict["maxvalues"] = maxvalues

    return(fitresults, fitresults_std, fiterrors, pdict)



###############################################################################
"""Return proper color-value for datapoint, given min/max values and current value"""
###############################################################################

def color(x, xmin, xmax):
    return(abs(255*(x-xmin)/(xmin-xmax)))


###############################################################################
"""Create image arrays from 3D array (intensity, position and linewidth)"""
###############################################################################

def create_mapsarray3D(fitresults, fiterrors, pdict):
    nx, ny = pdict["sample_size"]
    minvalues = pdict["minvalues"] 
    maxvalues = pdict["maxvalues"]

    # Create 3D array with map data, one layer per mapped parameter
    maps_array = np.zeros((nx,ny,3),dtype=np.uint8)
    
    x,y = 0,0
    # Loop over all three parameters
    for i in range(3): # Three parameters to be mapped: c, x0, lw
        for y in range(ny):    
            for x in range(nx):
                #print(fiterrors[x,y]) #debug
                
                # Leave out datapoints where the fit failed
                if fiterrors[x,y] == 0:
                    value = fitresults[x,y,i+1] # Fitresults contain: 0-b, 1-c, 2-x0, 3-lw. Require i = 1,2,3.
                    #print(value)
                    maps_array[x,y,i] = color(value, minvalues[i], maxvalues[i])

    return(maps_array)

###############################################################################
"""Create image array from 2D array"""
###############################################################################

def create_mapsarray2D(fitresults, pdict, exvalues):
    nx, ny = pdict["sample_size"]
    minvalue, maxvalue = exvalues
    
    # Create 3D array with map data, one layer per mapped parameter
    maps_array = np.zeros((nx,ny,1),dtype=np.uint8)
    
    x,y = 0,0
    for y in range(ny):    
        for x in range(nx):
            
            # Leave out datapoints where the fit failed
            if fitresults.mask[x,y] == False:
                value = fitresults[x,y] # Fitresults contain: 0-b, 1-c, 2-x0, 3-lw. Require i = 1,2,3.
                #print(value)
                maps_array[x,y,0] = color(value, minvalue, maxvalue)

    return(maps_array)

###############################################################################
"""Create image arrays of twist angle from fitresults"""
###############################################################################

def create_theta_map(pdict, fitresults, fiterrors, plot_size, save=False, folder=0, title=0):
    # Load phonon data
    dispersion = np.loadtxt("phondis_graphene.dat")
    
    # Save each phonon band into a separate list
    # Use only the K-to-Gamma-data (137:237)
    start, stop = 137, 237
    ta_branch = dispersion[start:stop,2]
    la_branch = dispersion[start:stop,3]

    # Invert lists to get Gamma-K instead of K-Gamma
    phon_branch = ta_branch[::-1]
    
    # Create proper y list containing twist angles
    gamma, k = 0, 4*pi/3/a #This is the distance Gamma-K in units of the graphene lattice constant a
    ylist = np.linspace(gamma, k, stop-start)
    # Convert ylist with crystal momenta Q_intra to corresponding twist angles
    thetas = [theta(q) for q in ylist]
    
    
    # Use scipy interpolation to create functions of twist angle as function
    # of phonon frequency
    order = "cubic"
    f_phon = interp1d(phon_branch, thetas, kind=order)
    
    # Load fit results
    nx, ny = pdict["sample_size"] 
    minvalues = pdict["minvalues"] 
    maxvalues = pdict["maxvalues"]
    min_theta, max_theta = f_phon(minvalues[1]), f_phon(maxvalues[1])
    #print(min_theta, max_theta)
    
    # Create 3D array with map data, one layer per mapped parameter
    map_array = np.zeros((nx,ny),dtype=np.uint8)
    theta_array = np.zeros((nx,ny))
    
    x,y = 0,0
    for y in range(ny):    
        for x in range(nx):
            #print(fiterrors[x,y]) #debug
            
            # Leave out datapoints where the fit failed
            if fiterrors[x,y] == 0:
                
                value = f_phon(fitresults[x,y,2]) # Fitresults contain: 0-b, 1-c, 2-x0, 3-lw. Require i = 1,2,3.
                theta_array[x,y] = value
                map_array[x,y] = color(value, min_theta, max_theta)

    
    # Set values and location for colormap labels
    ticklocs = np.linspace(round(np.amin(map_array),1), round(np.amax(map_array),1), 5)
    numbers = np.linspace(min_theta, max_theta, 5)
    ticklabels = ["{0:.1f}".format(n) for n in numbers]
    
    # Mask the array to exclude areas where fit failed
    maskedimage = np.ma.masked_where(fiterrors != 0, map_array)

    # Rotate the arrays to get them horizontal
    maskedimage = np.rot90(maskedimage)   
    
    # Do the plot
    fig = plt.figure(figsize=plot_size)
    if title==0:
        fig.suptitle("Twist angle (°) using TA branch")
    else:
        fig.suptitle("Twist angle (°) using TA branch" + title)
    ax = fig.add_subplot(111)
    ax.set_xticks([])
    ax.set_yticks([])
    cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
    im = ax.imshow(maskedimage, cmap=cmap, origin="lower")
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax, orientation="horizontal", cmap=cmap, ticks=ticklocs)
    cax.set_xticklabels(ticklabels, size=12) 

    if save == True:
        if folder == 0:
            if title==0:
                plt.savefig("Map_theta"+".png", format="png", dpi=900)
            else:
                plt.savefig("Map_theta "+str(title)+".png", format="png", dpi=900)
        else:
            if title==0:
                plt.savefig(str(folder)+"/"+"Map_theta"+".png", format="png", dpi=900)
            else:
                plt.savefig(str(folder)+"/"+"Map_theta "+str(title)+".png", format="png", dpi=900)
    return(theta_array)

###############################################################################   
"""Process image array to map"""
###############################################################################
# Creates a new colormap from a given one
def truncate_colormap(cmap, minval, maxval, n):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def create_maps(pdict, maps_data, fiterrors, plot_size, save, map_index="none", title=0, exvalues=0, folder=0, hx=0, hy=0, highlight=False):  
    nx, ny = pdict["sample_size"]
    minvalues = pdict["minvalues"]
    maxvalues = pdict["maxvalues"]
    peak_name = pdict["peak_name"]
    
    # Unless specified differently, build maps for all parameters
    if map_index == "none":
        print(0)
        n = len(maps_data[0,0,:])
        # Create array to save image arrays in, rotate them because the images 
        # are rotated by 90deg with respect to the data arrays  
        maps_images = np.zeros((nx,ny,n),dtype=np.uint8)
        maps_images = np.rot90(maps_images)
        
        for i in range(n):
            image = maps_data[:,:,i]   
            variable = pdict["variables"][i+1]
            ticklocs = np.linspace(0, np.amax(image), 5)
            
            # Format of the ticklabels depends on the variable (exp. or float)
            if exvalues!=0:
                minvalue, maxvalue = round(exvalues[0],1), round(exvalues[1],1)
                numbers = np.linspace(minvalue, maxvalue, 5)
                ticklabels = ["{0:.1f}".format(n) for n in numbers]
                
            elif variable == "Intensity" or variable =="Position":
                numbers = np.linspace(round(minvalues[i]), round(maxvalues[i]), 5)
                ticklabels = ["{0:.0f}".format(n) for n in numbers]
            else:
                numbers = np.linspace(round(minvalues[i],1), round(maxvalues[i],1), 5)
                ticklabels = ["{0:.1f}".format(n) for n in numbers]           
                
            # Mask the arrays to exclude areas where fit failed
            maskedimage = np.ma.masked_where(fiterrors != 0, image)

            # Rotate the arrays to get them horizontal
            maskedimage = np.rot90(maskedimage)   
    
            # Save final image data to maps_images
            maps_images[:,:,i] = maskedimage
            
            # Do the plots
            fig = plt.figure(figsize=plot_size)
            if title==0:
                fig.suptitle(str(peak_name)+"-Peak " + str(variable) + " (" + str(pdict["units"][i+1]) + ")")
            else:
                fig.suptitle(title)
                #fig.suptitle(str(peak_name)+"-Peak " + str(variable) + " (" + str(pdict["units"][i+1]) + ")" + title)
                
            ax = fig.add_subplot(111)
            ax.set_xticks([])
            ax.set_yticks([])
            cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
            im = ax.imshow(maskedimage, cmap=cmap, origin="lower")
            # Colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("bottom", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax, orientation="horizontal", cmap=cmap, ticks=ticklocs)
            cax.set_xticklabels(ticklabels, size=12)
        
            if save == True:
                if folder == 0:
                    if title==0:
                        plt.savefig(str(peak_name)+"_map_"+str(variable)+".png", format="png", dpi=900)
                    else:
                        plt.savefig(str(peak_name)+"_map_"+str(variable)+str(title)+".png", format="png", dpi=900)
                else:
                    if title==0:
                        plt.savefig(str(folder)+"/"+str(peak_name)+"_map_"+str(variable)+".png", format="png", dpi=900)
                    else:
                        plt.savefig(str(folder)+"/"+str(peak_name)+"_map_"+str(variable)+str(title)+".png", format="png", dpi=900)
                        #plt.savefig(str(folder)+"/"+str(peak_name)+"_map_"+str(variable)+str(title)+".png", format="png", dpi=900)
    
        # Return the image data stored in maps_images for further use
        return(maps_images)
            
            
    # If a single fitting parameter is specified using the map_index keyword, create a map only for this parameter
    else:
        
        image = maps_data[:,:,map_index]   
        variable = pdict["variables"][map_index+1]
        ticklocs = np.linspace(0, np.amax(image), 5)
        
        # Format of the ticklabels depends on the variable (exp. or float)
        if variable == "Intensity":
            numbers = np.linspace(minvalues[map_index], maxvalues[map_index], 5)
            ticklabels = ["{0:.2E}".format(n) for n in numbers]
        else:
            numbers = np.linspace(minvalues[map_index], maxvalues[map_index], 5)
            ticklabels = ["{0:.2f}".format(n) for n in numbers]
            
        # Mask the arrays to exclude areas where fit failed
        maskedimage = np.ma.masked_where(fiterrors != 0, image)

        # Rotate the arrays to get them horizontal
        maskedimage = np.rot90(maskedimage)   
        
        # Do the plots
        fig = plt.figure(figsize=plot_size)
        if title==0:
            fig.suptitle(str(peak_name)+"-Peak " + str(variable) + " (" + str(pdict["units"][map_index+1]) + ")")
        else:
            fig.suptitle(title)
        ax = fig.add_subplot(111)
        ax.set_xticks([])
        ax.set_yticks([])
        cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
        im = ax.imshow(maskedimage, cmap=cmap, origin="lower")
        # Colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax, orientation="horizontal", cmap=cmap, ticks=ticklocs)
        cax.set_xticklabels(ticklabels, size=12)

        # highlight a given position if required
        if highlight == True:
            position = np.zeros((nx, ny))
            position[:,hy] = 1
            position[hx,:] = 1
            maskedposition = np.ma.masked_where(position!=1, position)
            maskedposition = np.rot90(maskedposition)
            ax.imshow(maskedposition, cmap="binary_r", alpha=0.2, origin="lower")
    
        if save == True:
            if folder == 0:
                if title==0:
                    plt.savefig(str(peak_name)+"_map_"+str(variable)+".png", format="png", dpi=900)
                else:
                    plt.savefig(str(peak_name)+"_map_"+str(variable)+str(title)+".png", format="png", dpi=900)
            else:
                if title==0:
                    plt.savefig(str(folder)+"/"+str(peak_name)+"_map_"+str(variable)+".png", format="png", dpi=900)
                else:
                    plt.savefig(str(folder)+"/"+str(peak_name)+"_map_"+str(variable)+str(title)+".png", format="png", dpi=900)  


###############################################################################
"""For paper figures: Create maps directly from fitresults, more flexibility"""
###############################################################################
def papermap(pdict, data, mask, fiterrors, thresholds, plotsize, n_ticks, n_decimals, fontsize, scale2=False, data2=0, n_ticks2=0, n_decimals2=0, save=False, title=0, folder=0):
    nx, ny = data.shape
    lb, ub = thresholds
    
    # Mask data according to thresholds and errors, if desired
    masked_data = np.ma.masked_where(data > ub, data)
    masked_data = np.ma.masked_where(masked_data < lb, masked_data)
    if mask == True:
        masked_data = np.ma.masked_where(fiterrors != 0, masked_data)
    
    # Get min/maxvalues
    minvalue, maxvalue = np.amin(masked_data), np.amax(masked_data)
    
    # Create 3D array with map data, one layer per mapped parameter
    color_array = np.zeros((nx,ny),dtype=np.uint8)
    
    x,y = 0,0
    for y in range(ny):    
        for x in range(nx):
            # Leave out datapoints where the fit failed
            if masked_data.mask[x,y] == False:
                value = masked_data[x,y] # Fitresults contain: 0-b, 1-c, 2-x0, 3-lw. Require i = 1,2,3.
                color_array[x,y] = color(value, minvalue, maxvalue)
    
    # Mask data according to thresholds and errors, if desired
    masked_color_array = np.ma.masked_where(data > ub, color_array)
    masked_color_array = np.ma.masked_where(masked_data < lb, masked_color_array)
    if mask == True:
        masked_color_array = np.ma.masked_where(fiterrors != 0, masked_color_array)
    
    # Rotate the arrays to get them horizontal
    masked_color_array = np.rot90(masked_color_array)
    
    # Define position and value of ticks and labels
    numbers = np.linspace(minvalue, maxvalue, n_ticks)
    print(numbers)
    numbers = np.around(numbers, n_decimals)
    print(numbers)
    ticklocs = [255*(n-numbers[0])/(numbers[-1]-numbers[0]) for n in numbers]
    ticklabels = ["{0:.0f}".format(n) for n in numbers] 
    print(ticklabels)
         
    fig = plt.figure(figsize=plotsize)
    ax = fig.add_subplot(111)
    ax.set_xticks([])
    ax.set_yticks([])
    cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
    im = ax.imshow(masked_color_array, cmap=cmap, origin="lower")
    
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax, orientation="horizontal", cmap=cmap, ticks=ticklocs)
    pos = cbar.ax.get_position()
    cax.set_xticklabels(ticklabels, size=fontsize)
    
    if scale2==True:
        if mask == True:
            # The second dataset is thetas, already rotated, so rotate errors as well
            fiterrors = np.rot90(fiterrors,3)
            masked_data2 = np.ma.masked_where(fiterrors != 0, data2)
        
        # Get min/maxvalues
        minvalue2, maxvalue2 = np.amin(masked_data2), np.amax(masked_data2)
        
        # Define position and value of ticks and labels
        numbers2 = np.linspace(minvalue2, maxvalue2, n_ticks2)
        print(numbers2)
        numbers2 = np.around(numbers2, n_decimals2)
        print(numbers2)
        ticklocs2 = [255*(n-numbers2[0])/(numbers2[-1]-numbers2[0]) for n in numbers2]
        ticklabels2 = ["{0:.0f}".format(n) for n in numbers2] 
        print(ticklabels2)
        # Define new colorbar with both axes
        cax2 = cax.twiny()
        cax2.set_xticks(ticklocs2)
        cax2.set_xticklabels(ticklabels2, size=fontsize)
        cax2.set_position([0.166, 0.025, 0.667, 0.05]) # for 150*90 measurement
        #cax2.set_position([0.1, 0.025, 0.8, 0.05]) # for 100*100 measurement

    
    plt.tight_layout()

    if save == True:
        plt.savefig(str(folder)+"/" +"Papermap " + str(title)+".png", format="png", dpi=900, bbox_inches='tight')



###############################################################################       
"""Get fit and plot for certain peak"""
###############################################################################

def check_peak(xlist, data, fitresults, fiterrors, sample_size, plot_size, pdict, maps_data, mapindex, peak, hx, hy):   
    nx, ny = sample_size
    fitrange = pdict["fitrange"]
    p0 = pdict["startparams"]
    
    #Map with highlighted datapoint
    create_maps(pdict, maps_data, fiterrors, plot_size, False, mapindex, hx=hx, hy=hy, highlight=True)  

    ylist = data[hx,hy]

    # This expression searches the xlist for the values closest to the chosen fit_range and retrieves their indexes    
    start = min(range(len(xlist)), key=lambda i: abs(xlist[i]-fitrange[0]))
    stop = min(range(len(xlist)), key=lambda i: abs(xlist[i]-fitrange[1]))

    #Slice data to keep only the relevant peak
    ylist_fit = ylist[start:stop]
    xlist_fit = xlist[start:stop]

    #plot scan
    popt, popt_std = fit1L(xlist_fit, ylist_fit, p0, True, True)

    #Print whether or not this fit was accepted
    print("Result of the fit:")
    print(fiterrors[hx,hy])
        
    

###############################################################################       
"""Substract given fit from data set"""
###############################################################################

def substract_fit(sample_size, xlist, data, fitresults):
    nx,ny = sample_size
    data_substracted = np.zeros((nx,ny,len(data[0,0])))
    clist = []
    x0list = []
    lwlist = []
    blist = []

    for x in range(nx):
        for y in range(ny):
            # Read resulting parameters of the fit in this datapoint
            b = fitresults[x,y,0]
            c = fitresults[x,y,1]
            x0 = fitresults[x,y,2]
            lw = fitresults[x,y,3]
            blist.append(b)
            clist.append(c)
            x0list.append(x0)
            lwlist.append(lw)

            #create a list with y-values of the fit-function
            fitlist = [lorentzian(x, c, x0, lw)+b for x in xlist]

            # Substract those values from the data
            hlist = [y-yg for y,yg in zip(data[x,y], fitlist)]
            data_substracted[x,y] = hlist

    return(data_substracted)



###############################################################################    
"""Fancy Scatterplot of two peak properties, e.g. positions of G and R'"""
###############################################################################

def fancy_scatter_plot(fitdata1, fitdata2, index1, index2, dict_1, dict_2, error_array, title, save=False, folder=0):
    nx, ny = len(fitdata1[:,0,0]), len(fitdata1[0,:,0])

    #Read out relevant parameters from fitdata
    data1 = fitdata1[:,:,index1]
    data2 = fitdata2[:,:,index2]    

    #mask array where fits failed
    masked1 = np.ma.masked_where(data1 == 0, data1)
    masked2 = np.ma.masked_where(data2 == 0, data2)

    #Lists of markers and colors
    markerlist = [".", "^", "s", "p", "x", "P", "*", "D", "X", "o"]
    colorlist = ["deeppink", "violet", "darkviolet", "b", "deepskyblue", "darkturquoise", "mediumspringgreen", "greenyellow", "yellow", "lightsalmon", "tomato", "red", "crimson", "brown"]

    #Find out which rows and columns in the data contain nonzero elements
    nonzeroindices_y = [] #List of y-values that acutally show a peak
    nonzeroindices_x = [] #List of x-values that acutally show a peak
    
    for y in range(ny):
        if any(data2[:,y]):
            nonzeroindices_y.append(y)
    for x in range(nx):
        if any(data2[x,:]):
            nonzeroindices_x.append(x)
    
    #Scatter plot
    fig = plt.figure(figsize=(11,16))
    gridspec.GridSpec(15,10)
    ax1 = plt.subplot2grid((15,10), (0,0), colspan=7, rowspan=5)
    ax2 = plt.subplot2grid((15,10), (6,0), colspan=7, rowspan=5)
    ax1.set_xlabel("G position (1/cm)")
    ax1.set_ylabel("R' position (1/cm)")
    
    #step through the sample and add one point after another to the scatter plot
    for y in range(ny):
        for x in range(nx):
            
            if error_array[x,y] == 0:
                #c and m are used to assign a color and marker to each datapoint
                m = floor(len(markerlist)* (y-nonzeroindices_y[0]) / (nonzeroindices_y[-1] - nonzeroindices_y[0]) - 0.01)
                c = floor(len(colorlist)* (x-nonzeroindices_x[0]) / (nonzeroindices_x[-1] - nonzeroindices_x[0]) - 0.01)
                ax1.scatter([masked1[x,y]], [masked2[x,y]], marker=markerlist[m], color=colorlist[c])
                if x%2==0:
                    #only plot every second datapoint in the sample map, otherwise it gets too overloaded
                    ax2.scatter([x], [ny-y], marker=markerlist[m], color=colorlist[c])

    fig.suptitle(title)             
    #plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    if save == True:
        if folder == 0:
            plt.savefig("scatter_"+str(title)+str(index1)+str(index2)+".png", format="png", dpi=900)
        else:
            plt.savefig(str(folder)+"/"+"Scatter_"+str(title)+".png", format="png", dpi=900)

     

###############################################################################    
"""Scatterplot of two peak properties, optionally color coded"""
###############################################################################
def scatter_plot(data1, data2, figsize, xlim, ylim, title, xlabel, ylabel, color_code=False, color_data=0, color_threshold=0, errorbars=False, std_data=0, widthgradprediction=False, save=False, folder=0):
    # Set up scatter plot
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(111) 
   
    ax1.set_xlabel(str(xlabel))
    ax1.set_ylabel(str(ylabel))
    fig.suptitle(title)     
    
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    
    if errorbars == True:
        fraction = 1 # Determines how many data points will be in the errorbar plot, e.g. 3 --> every third data point    
        data1 = np.ndarray.flatten(data1[0::fraction, 0::fraction])
        data2 = np.ndarray.flatten(data2[0::fraction, 0::fraction])
        std_data = np.ndarray.flatten(std_data[0::fraction, 0::fraction])
    
    if color_code == False: 
        if errorbars == False:
            ax1.scatter(data1, data2, color="b", marker=".", s=1)
            
        if errorbars == True:
            ax1.errorbar(data1, data2, std_data, fmt=".", color="cornflowerblue", markersize=5, capsize=5, zorder=0)
            ax1.scatter(data1, data2, color="b", s=5, zorder=1)
            
    if color_code == True:
        
        # Rotate color array if necessary
        if data1.shape != color_data.shape:
            print("rotacion")
            color_data = np.rot90(color_data)
        
        # Mask data arrays according to color threshold
        data1 = ma.masked_where(color_data < color_threshold, data1)
        data2 = ma.masked_where(color_data < color_threshold, data2)
        
        # Define color map
        cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
        
        if errorbars == False:
            ax1.scatter(data1, data2, s=2, c=color_data, cmap=cmap)
        
        if errorbars == True:
            std_data = np.ma.masked_where(np.ma.getmask(data1), std_data)
            ax1.errorbar(data1, data2, std_data, fmt=".", markersize=5, capsize=5, c=color_data, cmap=cmap)    
    
    
    if widthgradprediction == True:
        # Load calculation of Grad(Pos(TA)) from Jens
        # Transform to Grad(Theta) by slope factor m obtained by linear fit of Theta - Pos(TA)
        # Transformation has maximum deviation of 0.014deg in range of 220-340 1/cm / 5.8-9.5 deg
        m = 0.02732685075118969
        lines = tuple(open("width_vs_grad_2FWHM.txt", 'r'))
        data_str1 = [x.split('\n') for x in lines]
        data_str2 = [[x.split() for x in stringlist] for stringlist in data_str1]
        xlist = [m*float(x)*1000 for x in data_str2[0][0]]
        # Set y axis interception to 1.08333
        ylist = [float(y) for y in data_str2[1][0]] #Factor 1/2 converts FWHM to linewidth 
        ax1.plot(xlist, ylist, label="Theoretical value", color="red")
        #plt.legend()
    
    if save == True:
        if folder == 0:
            if errorbars == False:
                plt.savefig("Scatter_"+str(title)+".png", format="png", dpi=900, bbox_inches='tight')
            if errorbars == True:
                plt.savefig("Errorbar_"+str(title)+".png", format="png", dpi=900, bbox_inches='tight')
                
        if folder != 0:
            if errorbars == False:
                plt.savefig(str(folder)+"/"+"Scatter_"+str(title)+".png", format="png", dpi=900, bbox_inches='tight')
            if errorbars == True:
                plt.savefig(str(folder)+"/"+"Errorbar_"+str(title)+".png", format="png", dpi=900, bbox_inches='tight')





###############################################################################    
"""Map of distance between two peaks"""
###############################################################################

def difference_map(pdict_1, pdict_2, fitdata1, fitdata2, error1, error2, variable_index, plot_size, save=False, folder=0):
    nx, ny = pdict_1["sample_size"]
    peak_name_1 = pdict_1["peak_name"]
    peak_name_2 = pdict_2["peak_name"]
    variable = pdict_1["variables"][variable_index]
    
    # Read out data
    position1 = fitdata1[:,:,variable_index]
    position2 = fitdata2[:,:,variable_index]
    
    #create array of differences
    difference = np.subtract(position2, position1)
    masked_difference = np.ma.masked_where(error1 != 0, difference)
    masked_difference = np.ma.masked_where(error2 != 0, masked_difference)
    dmax = np.amax(masked_difference)
    dmin = np.amin(masked_difference)
    #create array of differences fit for plotting, mask where there is no R'-peak
    image = np.zeros((nx,ny))

    x,y = 0,0
    for y in range(ny):    
        for x in range(nx):
            d = difference[x,y]
            image[x,y] = color(d, dmin, dmax)
    masked_image = np.ma.masked_where(error1 != 0, image)
    masked_image = np.ma.masked_where(error2 != 0, masked_image)
    masked_image = np.rot90(masked_image)    

    #Define ticks on colorbar
    ticklocs = np.linspace(0, np.amax(masked_image), 5)
    numbers = np.linspace(np.amin(masked_difference), np.amax(masked_difference), 5)
    ticklabels = ["{0:.2f}".format(n) for n in numbers]

    # Do the plots
    fig = plt.figure(figsize=plot_size)
    fig.suptitle(str(peak_name_1)+" - "+str(peak_name_2)+" "+str(variable))
    ax = fig.add_subplot(111)
    ax.set_xticks([])
    ax.set_yticks([])
    cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
    im = ax.imshow(masked_image, cmap=cmap, origin="lower")
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax, orientation="horizontal", cmap=cmap, ticks=ticklocs)
    cax.set_xticklabels(ticklabels, size=12)

    if save == True:
        if folder == 0:
            plt.savefig("Difference "+str(peak_name_1)+" - "+str(peak_name_2)+" "+str(variable)+".png", format="png", dpi=900)        
        else:
            plt.savefig(str(folder)+"/"+"Difference "+str(peak_name_1)+" - "+str(peak_name_2)+" "+str(variable)+".png", format="png", dpi=900)        

    return(masked_difference)


###############################################################################    
"""Map of ratio between two peak quantities"""
###############################################################################

def ratio_map(pdict_1, pdict_2, fitdata1, fitdata2, error1, error2, variable_index, plot_size, save=False, folder=0):
    nx, ny = pdict_1["sample_size"]
    peak_name_1 = pdict_1["peak_name"]
    peak_name_2 = pdict_2["peak_name"]
    variable = pdict_1["variables"][variable_index]
    
    # Read out data
    value1 = fitdata1[:,:,variable_index]
    value2 = fitdata2[:,:,variable_index]
    
    #create array of differences
    ratio = np.divide(value2, value1)
    masked_ratio = np.ma.masked_where(error1 != 0, ratio)
    masked_ratio = np.ma.masked_where(error2 != 0, masked_ratio)
    dmax = np.amax(masked_ratio)
    dmin = np.amin(masked_ratio)
    #create array of differences fit for plotting, mask where there is no R'-peak
    image = np.zeros((nx,ny))

    x,y = 0,0
    for y in range(ny):    
        for x in range(nx):
            d = ratio[x,y]
            image[x,y] = color(d, dmin, dmax)
    masked_image = np.ma.masked_where(error1 != 0, image)
    masked_image = np.ma.masked_where(error2 != 0, masked_image)
    masked_image = np.rot90(masked_image)    

    #Define ticks on colorbar
    ticklocs = np.linspace(0, np.amax(masked_image), 5)
    numbers = np.linspace(np.amin(masked_ratio), np.amax(masked_ratio), 5)
    ticklabels = ["{0:.2f}".format(n) for n in numbers]

    # Do the plots
    fig = plt.figure(figsize=plot_size)
    fig.suptitle(str(peak_name_1)+" - "+str(peak_name_2)+" "+str(variable))
    ax = fig.add_subplot(111)
    ax.set_xticks([])
    ax.set_yticks([])
    cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
    im = ax.imshow(masked_image, cmap=cmap, origin="lower")
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax, orientation="horizontal", cmap=cmap, ticks=ticklocs)
    cax.set_xticklabels(ticklabels, size=12)

    if save == True:
        if folder == 0:
            plt.savefig("Difference "+str(peak_name_1)+" - "+str(peak_name_2)+" "+str(variable)+".png", format="png", dpi=900)        
        else:
            plt.savefig(str(folder)+"/"+"Difference "+str(peak_name_1)+" - "+str(peak_name_2)+" "+str(variable)+".png", format="png", dpi=900)



###############################################################################
"""Get the mean value of any fitting parameter in a certain area of the sample"""
###############################################################################

def get_mean(fitdata, errordata, mean_index, pdict, maps_images, x_bounds, y_bounds, cond_array, cond_var, lb, ub, save, folder=0):
    nx, ny = pdict["sample_size"]
    minvalues = pdict["minvalues"]
    maxvalues = pdict["maxvalues"]
    peak_name = pdict["peak_name"]
    parameter = pdict["variables"][mean_index]
    unit = pdict["units"][mean_index]
    
    # Roughly define area from which to extract the data
    xmin, xmax  = x_bounds
    ymin, ymax = y_bounds
    
    # Define array that will mark the used datapoints
    marker_array = np.zeros((nx,ny))
    
    # Rotate color array if necessary
# =============================================================================
#     if fitdata.shape != maps_images.shape:
#         print("rotacion")
#         maps_images = np.rot90(maps_images,4)
# =============================================================================
    
    # Mark used datapoints in marker_array
    for x in range(nx):
        for y in range(ny):
            
            # First condition: Fit must have worked
            if errordata[x,y] == 0:
                # Second condition: Datapoint must be inside the area that is to be considered
                if xmin <= x <= xmax and ymin <= ny-y <= ymax:
                    # Third condition: Defined by cond_array, cond_var
                    if lb < cond_array[x,y,cond_var] < ub:
                        marker_array[x,y] = 1

    
    # Take mean, std, min, max of considered datapoints
    mean = np.mean(fitdata[:,:,mean_index][marker_array > 0])
    std = np.std(fitdata[:,:,mean_index][marker_array > 0])
    minvalue_area = np.amin(fitdata[:,:,mean_index][marker_array > 0])
    maxvalue_area = np.amax(fitdata[:,:,mean_index][marker_array > 0])

    # Define the scale bar
    maxvalue = maxvalues[mean_index-1]
    minvalue = minvalues[mean_index-1]
    ticklocs = [0, 255]
    ticklabels = ["{0:.1f} ".format(minvalue), "{0:.1f} ".format(maxvalue)]
    
    # Extract the correct map from the map image data
    image = maps_images[:,:,mean_index-1]
    image = np.rot90(image)
    image = np.rot90(image)
    image = np.rot90(image)
    
    # Mask the array to exclude areas where fit failed
    masked_image = np.ma.masked_where(errordata != 0, image)

    # Rotate the arrays to get them horizontal
    masked_image = np.rot90(masked_image)
    marker_array = np.rot90(marker_array)     

    # Build plot
    fig_mean = plt.figure(figsize=(5,5))
    fig_mean.suptitle( str(peak_name) + " " + str(parameter) + ": " + "({0:.2f} +/- {1:.2f}) ".format(mean, std) + str(unit) + "\n Min/Max: {0:.2f}, {1:.2f}".format(minvalue_area, maxvalue_area), fontsize=12)
    ax = fig_mean.add_subplot(111)
    ax.set_xticks([])
    ax.set_yticks([])
        
    cmap = truncate_colormap(plt.get_cmap("gist_ncar"), 0.05, 0.8, 100)
    im = ax.imshow(masked_image, cmap=cmap, origin="lower")
    ax.imshow(marker_array, cmap = "Greys", alpha = 0.5, origin="lower")
    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax, orientation="horizontal", cmap=cmap, ticks=ticklocs)
    cax.set_xticklabels(ticklabels, size=12)

    if save == True:
        if folder == 0:
            plt.savefig(str(peak_name)+"_mean_"+str(parameter)+".png", format="png", dpi=900)
        else:
            plt.savefig(str(folder) +"/" + str(peak_name)+"_mean_"+str(parameter)+".png", format="png", dpi=900)


    return(mean, std)



###############################################################################
"""Get the ratio of any two fitting parameters in a certain area of the sample"""
###############################################################################

def get_ratio(sample_size, fitresults_1, fitresults_2, index_1, index_2):
    ratiolist = []
    nx, ny = sample_size

    for x in range(nx):
        for y in range(ny):            

            if fitresults_1[x,y,1] != 0:
                ratio = fitresults_2[x,y,index_1] / fitresults_1[x,y,index_2]
                ratiolist.append(ratio)

    ratio_mean = np.mean(ratiolist)
    ratio_std = np.std(ratiolist)

    print("Ratio = {0:.2E} +/- {1:.2E}".format(ratio_mean, ratio_std))

    return(ratio_mean, ratio_std)




###############################################################################
"""Normalize a given dataset, relevant for AI training"""
###############################################################################

def norm_data(df):
    normalized_df=(df-df.mean())/df.std()
    return(normalized_df)


###############################################################################
"""Compute averaged gradient of an array"""
###############################################################################
def averaged_gradient(data, mask_array, reach, stepsize):

    # Mask the array to exclude areas where the fit failed
    data = ma.masked_where(mask_array != 0, data)
    
    # Obtain first neighbour gradient by numpy
    gradient_x, gradient_y = np.gradient(data)

# =============================================================================
#     # Combine averaged gradient_x and gradient_y to an absolute value
#     gradient_abs_av = ma.sqrt(ma.power(gradient_x, 2) + ma.power(gradient_y, 2))/stepsize
#     gradient_abs_av = moving_average(gradient_abs_av, reach)
# =============================================================================
    
    
    # Average the x and y gradient components
    gradient_x_av = moving_average(gradient_x, reach)
    gradient_y_av = moving_average(gradient_y, reach)
    
    # Combine averaged gradient_x and gradient_y to an absolute value
    gradient_abs_av = ma.sqrt(ma.power(gradient_x_av, 2) + ma.power(gradient_y_av, 2))/(stepsize)

    # Create new error array that screens the edges and spots where moving_average yields no results
    gradient_abs_av.mask[0:reach,:] = 1
    gradient_abs_av.mask[-reach:,:] = 1
    gradient_abs_av.mask[:,0:reach] = 1
    gradient_abs_av.mask[:,-reach:] = 1
    
    # Extract extremal values of the gradient array
    exvalues = np.amin(gradient_abs_av), np.amax(gradient_abs_av)

    return(gradient_abs_av, exvalues)




###############################################################################
"""Average a given array in a given reach around the current data point"""
###############################################################################
def moving_average(data_array, reach, mask=False, error_data=0):
    nx,ny = data_array.shape
    averaged_array = ma.array(np.zeros((nx,ny)), mask = np.zeros((nx,ny)))
    
    # If an array with fiterrors is provided, use it to mask the real data
    if mask == True:
        data_array = ma.masked_where(error_data != 0, data_array)
    
    # Scan over the array, avoiding the edges so far as to always have #reach neighbours around each data point in which a gradient is computed
    for x in range(nx-2*reach):
        x = x+reach
        
        for y in range(ny-2*reach):
            y = y+reach
            fitfails = 0
            
            # Compute the average in the viscinity of the current datapoint, considering all other data points that lie in a square with
            # side length 2*reach+1 around the current data point
            dividend = 0
            divisor = 0
            
            # reach+1 makes sure that e.g. for reach=1, i = -1,0,1
            for i in range(-reach, reach+1):
            
                for j in range(-reach, reach+1):
                    # Do not include the local data point into the calculation of the local average
                    if i == 0  and j==0: continue
                    
                    # If there is a masked element in the environment of the current data point:
                    # Raise the counter of fit fails. If at the end there have been >=4 fitfails
                    # around the current data points, skip it.
                    if data_array.mask[x+i, y+j] == True:
                        fitfails += 1
                        continue
                    
                    else:
                        value = data_array[x+i, y+j]
                        weight = 1/sqrt(i**2 + j**2)
                        
                        dividend += value*weight
                        divisor += weight
                    
            # If at the end there have been <= 4 fitfails, gradient can be calculated
            if fitfails <= 3:
                averaged_array[x,y] = dividend/divisor
            
            # Otherwise skip the current data point and mask it 
            else:
                averaged_array.mask[x,y] = True
                
            
    # Apply mask around the edge where no average value was calculated
    averaged_array.mask[0:reach,:] = 1
    averaged_array.mask[-reach:,:] = 1
    averaged_array.mask[:,0:reach] = 1
    averaged_array.mask[:,-reach:] = 1
                
    return(averaged_array)        

