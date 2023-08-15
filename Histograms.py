#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 14:20:56 2023

@author: Aaron
"""
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from numpy import pi, arcsin, sqrt
from scipy.interpolate import interp1d


np.set_printoptions(threshold=np.inf)

a = 0.246 # Graphene lattice constant in nm


###############################################################################
""" 1. User input """
###############################################################################

# Directory of the Raman data
folder = "/Users/Aaron/Desktop/Code/Test_data" 

# When loading the results from an old fit: Which peaks shall be loaded?
b_load_TA = True
b_load_G = False
b_load_LO = False
b_load_2D = False

# Which peaks shall be mapped?
b_hist_TA = True
b_hist_G = False
b_hist_LO = False
b_hist_2D = False

# Number of bins in the histograms
N = 50

# Ranges covered on the x axes of the histograms
histrange_TA_int = (10, 500)
histrange_TA_pos = (250, 265)
histrange_TA_lw = (0.1, 5)
histrange_G_int = ()
histrange_G_pos = ()
histrange_G_lw = ()
histrange_LO_int = ()
histrange_LO_pos = ()
histrange_LO_lw = ()
histrange_2D_int = ()
histrange_2D_pos = ()
histrange_2D_lw = ()
histrange_theta = (5, 8)
histrange_grad_theta = ()


###############################################################################
""" 2. Don't touch section A (functions) """
###############################################################################

""" Lorentzian as used in spectroscopy """
def lorentzian(x, b, c, x0, lw):
   f = b + c / pi * (lw/2) / ((x-x0)**2 + (lw/2)**2)
   return(f)


""" Function to build interactive plot window """
def make_figure(fitresults, fiterrors, pdict):
    
    # Retrieve required variables from the dictionnary
    nx, ny = pdict["size_px"] # Scan size in pixels
    sx, sy = pdict["size_um"] # Scan size in microns
    peakname = pdict["peakname"] # Name of the peak
    (thresh_c, thresh_x0, thresh_lw) = pdict["params_thresh"] # Threshold values
    crange_int = pdict["crange_int"] # Colorbar range of intensity
    crange_pos = pdict["crange_pos"] # Colorbar range of position
    crange_lw = pdict["crange_lw"] # Colorbar range of linewidth
    
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
    fig = plt.figure(figsize = (12, 12))
    grid = mpl.gridspec.GridSpec(6, 6, figure = fig, hspace=0.6, wspace=0.8)
    ax0 = fig.add_subplot(grid[:3, 0:3])
    ax1 = fig.add_subplot(grid[:3, 3:6])
    ax2 = fig.add_subplot(grid[3:6, 0:3])
    ax3 = fig.add_subplot(grid[3:6, 3:6])
    
    # Plot the maps
    # Check if colorbar ranges were specified, if so: Apply them.
    if crange_int != ():
        im0 = ax0.imshow(fitresults_ma[:,:,1], extent = [0, sx, 0, sy], cmap="gist_rainbow", vmin=crange_int[0], vmax=crange_int[1]) # Intensity
    else:
        im0 = ax0.imshow(fitresults_ma[:,:,1], extent = [0, sx, 0, sy], cmap="gist_rainbow") # Intensity
    
    if crange_pos != ():
        im1 = ax1.imshow(fitresults_ma[:,:,2], extent = [0, sx, 0, sy], cmap="gist_rainbow", vmin=crange_pos[0], vmax=crange_pos[1]) # Position
    else:
        im1 = ax1.imshow(fitresults_ma[:,:,2], extent = [0, sx, 0, sy], cmap="gist_rainbow") # Position
    
    if crange_lw != ():
        im2 = ax2.imshow(fitresults_ma[:,:,3], extent = [0, sx, 0, sy], cmap="gist_rainbow", vmin=crange_lw[0], vmax=crange_lw[1]) # Linewidth
    else:
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


""" Handle click events """
# This function is activated when double clicking in tne interactive map.
# It creates a matplotlib path object from the clicked points in the map.
def onclick(event, fitresults, fiterrors, pdict):
    if event.dblclick:
        if event.button == 1:
            if isinstance(event.ydata, float): # Check if the doubleclick event happened inside the plot and returns correct values
                # Exact location of click in axis units
                x_value, y_value = event.xdata, event.ydata
                
                # Find index coordinates of click position. The +1 at y_map en-
                # sures that the shown spectrum really belongs to the clicked
                # pixel, otherwise there is a mismatch of exactly one pixel.
                x_index = np.abs(xaxis - x_value).argmin()
                y_index = np.abs(yaxis - y_value).argmin() + 1
                                
                # Add clicked point to the list of path points
                pathpoints_value.append((x_value, y_value))
                pathpoints_index.append((x_index, y_index))
                                
                # Mark the clicked point in the maps
                ax0.scatter([x_value], [y_value], s=15, color="black")
                ax1.scatter([x_value], [y_value], s=15, color="black")
                ax2.scatter([x_value], [y_value], s=15, color="black")
                
                # Two things happen only if any point has been clicked before:
                # a) Draw a line from the last point to the current one
                # b) Check if the current point closes the path
                if len(pathpoints_index) > 1:
                    
                    # Get coordinates of last point
                    x_value_old, y_value_old = pathpoints_value[-2]
                    x_index_old, y_index_old = pathpoints_index[-2]
                    
                    # Plot the connecting lines in all maps
                    ax0.plot([x_value, x_value_old], [y_value, y_value_old], color="black")
                    ax1.plot([x_value, x_value_old], [y_value, y_value_old], color="black")
                    ax2.plot([x_value, x_value_old], [y_value, y_value_old], color="black")
                
                    # Get coordinates of first point
                    x_index_1st, y_index_1st = pathpoints_index[0]
                
                    # Check if the clicked point is close to the first point that 
                    # was clicked, thereby closing the path. In that case, proceed
                    # to creating the histogram.
                    if np.abs(x_index - x_index_1st) < 2 and np.abs(y_index - y_index_1st) < 2:
                        make_histograms(fitresults, fiterrors, pdict)
                            
                # Update the plot in the window
                ax3.set_xlabel("Raman shift (rel. 1/cm)")
                ax3.set_ylabel("CCD counts")
                ax3.figure.canvas.draw()
                
            else:
                print('Please double click inside a map!')
                return None


""" Plot and save histograms """
# This function is called when the user closes the path in the interactive plot.
# It creates histograms of the Lorentz parameters considering only the area inside
# the created path. If the TA peak was selected, it also makes a histogram of the
# twist angle.
def make_histograms(fitresults, fiterrors, pdict, theta=0):
    # Retrieve required variables from the dictionnary
    nx, ny = pdict["size_px"] # Scan size in pixels
    sx, sy = pdict["size_um"] # Scan size in microns
    peakname = pdict["peakname"] # Name of the peak
    (thresh_c, thresh_x0, thresh_lw) = pdict["params_thresh"] # Threshold values
    histrange_int = pdict["histrange_int"]
    histrange_pos = pdict["histrange_pos"]
    histrange_lw = pdict["histrange_lw"]

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
    
    # Create the fiterror mask in 2D
    # MASK = TRUE = 1 MEANS DATA POINT IS MASKED
    # MASK = FALSE = 0 MEANS DATA POINT IS NOT MASKED
    # Convert to int so it can be combined with the path mask which is also int
    mask_fiterrors = np.logical_or.reduce(conditions).astype(int)
    
    # Convert the list of clicked coordinates into a path
    userpath = Path(pathpoints_index)
    
    # Create a grid of coordinates
    x_grid, y_grid = np.meshgrid(np.arange(0, nx), np.arange(0, ny))
    
    # Convert grid coordinates to 1D arrays
    x_coords = x_grid.flatten()
    y_coords = y_grid.flatten()
    
    # Check if each coordinate is inside the path, resulting in a 1D mask
    # The 1 - makes sure that points inside the path get assigned a 0 (not masked)
    # and points outside the path get assigned a 1 (masked)
    mask_path_1D = 1 - userpath.contains_points(np.column_stack((x_coords, y_coords)))
    
    # Reshape the boolean mask to match the grid shape
    mask_path_2D = mask_path_1D.reshape(x_grid.shape)
    
    # Turn the mask upside down so it corresponds to the user input (otherwise
    # it stands upside down)
    mask_path_2D = np.flipud(mask_path_2D)
        
    # Combine the fiterror and the path mask
    mask = np.ma.mask_or(mask_fiterrors, mask_path_2D)

    # The mapped quantities
    quantities = ["intensity", "position", "linewidth"]
    # Their units
    x_labels = ["Intensity (arb. u.)", r"$\omega$ (1/cm)", r"$\Gamma$ (1/cm)"]
    # Their ranges
    histranges = [histrange_int, histrange_pos, histrange_lw]
    
    # In a loop, create the histograms of the Lorentz parameters
    for i in range(3):
        hist_data = (fitresults[:,:,i+1]*(1-mask)).flatten()
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        ax.hist(hist_data, bins=N, range=histranges[i])
        ax.set_xlabel(x_labels[i])
        ax.set_ylabel("N")
        fig.suptitle(peakname + " " + quantities[i])
        plt.savefig(folder+"/"+peakname + "_hist_" + quantities[i]+".pdf", format="pdf", dpi=300)
     
    # If an array of twist angles was provided, make a histogram of that too
    if b_hist_TA:
        
        ############################### Part 1 ####################################
        # Convert peak position to twist angle #
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
        
        # Create an array of TA position values where all values above the 
        # interpolation range maximum (748 1/cm) are replaced by a dummy value of 0
        posTA_cutoff = np.where(fitresults[:,:,2] > 748, 0, fitresults[:,:,2])
        
        # Calculate the twist angle array
        theta = TA_position_to_theta(posTA_cutoff)
        
        # Make the histogram
        hist_data = (theta*(1-mask)).flatten()
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        ax.hist(hist_data, bins=N, range=histrange_theta)
        ax.set_xlabel("Twist angle (°)")
        ax.set_ylabel("N")
        plt.savefig(folder+"/twist_angle_hist.pdf", format="pdf", dpi=300)
     
    return()


###############################################################################
""" 4. Don't touch section B (executing the code) """
###############################################################################
    
###############################################################################
# 4.1 Load results of a previous fit along with xdata and the peak dictionnary
###############################################################################
    
if b_load_TA == True:
    with open(folder+"/TA", "rb") as file:
        xdata, data, pdict_TA, fitresults_TA, fitresults_std_TA, fiterrors_TA = pickle.load(file)
      
if b_load_G == True:
    with open(folder+"/G", "rb") as file:
        xdata, data, pdict_G, fitresults_G, fitresults_std_G, fiterrors_G = pickle.load(file)   
    
if b_load_LO == True:
    with open(folder+"/LO", "rb") as file:
        xdata, data, pdict_LO, fitresults_LO, fitresults_std_LO, fiterrors_LO = pickle.load(file)

if b_load_2D == True:
    with open(folder+"/2D", "rb") as file:
        xdata, data, pdict_2D, fitresults_2D, fitresults_std_2D, fiterrors_2D = pickle.load(file)


###############################################################################
# 4.3 Open the interactive figures
###############################################################################
if b_hist_TA == True:
    # Add the range of the histogram to the pdict
    pdict_TA["histrange_int"] = histrange_TA_int
    pdict_TA["histrange_pos"] = histrange_TA_pos
    pdict_TA["histrange_lw"] = histrange_TA_lw
    
    # Extract map size from pdict, needed at the end of the code!
    nx, ny = pdict_TA["size_px"] # Scan size in pixels
    sx, sy = pdict_TA["size_um"] # Scan size in microns
    
    # Build the interactive plot
    ax0, ax1, ax2, ax3 = make_figure(fitresults_TA, fiterrors_TA, pdict_TA)
    
if b_hist_G == True:
    # Add the range of the histogram to the pdict
    pdict_G["histrange_int"] = histrange_G_int
    pdict_G["histrange_pos"] = histrange_G_pos
    pdict_G["histrange_lw"] = histrange_G_lw
    
    # Extract map size from pdict, needed at the end of the code!
    nx, ny = pdict_G["size_px"] # Scan size in pixels
    sx, sy = pdict_G["size_um"] # Scan size in microns
    
    # Build the interactive plot
    ax0, ax1, ax2, ax3 = make_figure(fitresults_G, fiterrors_G, pdict_G)
    
if b_hist_LO == True:
    # Add the range of the histogram to the pdict
    pdict_LO["histrange_int"] = histrange_LO_int
    pdict_LO["histrange_pos"] = histrange_LO_pos
    pdict_LO["histrange_lw"] = histrange_LO_lw
    
    # Extract map size from pdict, needed at the end of the code!
    nx, ny = pdict_LO["size_px"] # Scan size in pixels
    sx, sy = pdict_LO["size_um"] # Scan size in microns
    
    # Build the interactive plot
    ax0, ax1, ax2, ax3 = make_figure(fitresults_LO, fiterrors_LO, pdict_LO)
    
if b_hist_2D == True:
    # Add the range of the histogram to the pdict
    pdict_2D["histrange_int"] = histrange_2D_int
    pdict_2D["histrange_pos"] = histrange_2D_pos
    pdict_2D["histrange_lw"] = histrange_2D_lw
    
    # Extract map size from pdict, needed at the end of the code!
    nx, ny = pdict_2D["size_px"] # Scan size in pixels
    sx, sy = pdict_2D["size_um"] # Scan size in microns
    
    # Build the interactive plot
    ax0, ax1, ax2, ax3 = make_figure(fitresults_2D, fiterrors_2D, pdict_2D)
    
# Create lists that are needed to locate the index coordinates of a click event
xaxis = np.linspace(0, sx, nx)
yaxis = np.linspace(0, sy, ny)

# Create list of clicked points that will together form the edge of the region
# we need a histogram of
pathpoints_index = []
pathpoints_value = []



    
    