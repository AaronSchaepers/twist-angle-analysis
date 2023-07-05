#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:07:14 2023

@author: Aaron
@date: 2023/06/29
@github: https://github.com/AaronSchaepers/twist-angle-analysis

Required packages:
    pyqtgraph
    
Next step:
    - Make the plots show up
    
At all times, there is one peak instance for TA, G, LO and 2D.
The peak selected in the dropdown menu needs to act as a kind of condition that 
defines on which peak instance all the other actions (set default values, fit, map)
are performed.  
"""

import sys
import pickle
import default_values

import numpy as np
import pyqtgraph as pg
import module as mod

from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QLineEdit, QPushButton, QComboBox, QMainWindow, QGraphicsView, QGraphicsItem


###############################################################################
""" The main window """
###############################################################################
# Create a base class that will load the .ui file in the constructor
# It inherits the QtWidgets.QMainWindow class because I created a new 
# "Main Window" when selecting the form type when first creating the .ui file 
# in PyQt Designer
class GUI(QMainWindow):
    def __init__(self):
        
        # Call the inherited classes __init__ method
        super(GUI, self).__init__() 
        
        # Load the .ui file
        uic.loadUi('gui.ui', self) 
        
        # Connect the plot frames to attributes of this class
        self.view_int = self.findChild(QGraphicsView, 'view_int')
        self.view_pos = self.findChild(QGraphicsView, 'view_pos')
        self.view_lw = self.findChild(QGraphicsView, 'view_lw')
        self.view_theta = self.findChild(QGraphicsView, 'view_theta')
        self.view_grad_theta = self.findChild(QGraphicsView, 'view_grad_theta')
        self.view_spectrum = self.findChild(QGraphicsView, 'view_spectrum')
        
        # Connect the GUI input elements to attributes of this class
        self.input_min_mean_range = self.findChild(QLineEdit, 'in_min_mean_range')
        self.input_max_mean_range = self.findChild(QLineEdit, 'in_max_mean_range')
        self.input_sx_px = self.findChild(QLineEdit, 'in_sx_px')
        self.input_sy_px = self.findChild(QLineEdit, 'in_sy_px')
        self.input_sx_mu = self.findChild(QLineEdit, 'in_sx_mu')
        self.input_sy_mu = self.findChild(QLineEdit, 'in_sy_mu')
        self.input_init_int = self.findChild(QLineEdit, 'in_init_int')
        self.input_init_lw = self.findChild(QLineEdit, 'in_init_lw')
        self.input_min_int = self.findChild(QLineEdit, 'in_min_int')
        self.input_max_int = self.findChild(QLineEdit, 'in_max_int')
        self.input_min_pos = self.findChild(QLineEdit, 'in_min_pos')
        self.input_max_pos = self.findChild(QLineEdit, 'in_max_pos')
        self.input_min_lw = self.findChild(QLineEdit, 'in_min_lw')
        self.input_max_lw = self.findChild(QLineEdit, 'in_max_lw')
        
        # Connect the "Import Raman scan" button to a method of this class
        self.import_scan_button = self.findChild(QPushButton, 'b_import_scan')
        self.import_scan_button.clicked.connect(self.import_button_clicked)
        
        # Connect peak dropdown menu to a method of this class
        self.which_peak = self.findChild(QComboBox, 'in_peak')
        
        # Connect the "Use default values" button to a method of this class
        self.default_values_button = self.findChild(QPushButton, 'b_default_values')
        self.default_values_button.clicked.connect(self.default_values_button_clicked)
        
        # Connect the "Start fit" button to a method of this class
        self.start_fit_button = self.findChild(QPushButton, 'b_start_fit')
        self.start_fit_button.clicked.connect(self.start_fit_button_clicked)
        
        # Connect the "Save fit" button to a method of this class
        self.save_fit_button = self.findChild(QPushButton, 'b_save_fit')
        self.save_fit_button.clicked.connect(self.save_fit_button_clicked)
        
        # Connect the "Load fit" button to a method of this class
        self.import_fit_button = self.findChild(QPushButton, 'b_import_fit')
        self.import_fit_button.clicked.connect(self.import_fit_button_clicked)
                
        # Show the GUI
        self.show()
        
    # Upon button click, create an instance of the importer class, make it an
    # attribute of the GUI instance and load Raman scan data into it
    def import_button_clicked(self):
        self.data = importer()
    
    # Upon button click, set default values in the text boxes
    def default_values_button_clicked(self):
        # Get the default pdict for the current peak
        default_pdict = default_pdicts[self.which_peak.currentIndex()]
        # Display default values in the text boxes, convert float to string
        self.input_init_int.setText("{0:.0f}".format(default_pdict["init_int"]))
        self.input_init_lw.setText("{0:.1f}".format(default_pdict["init_lw"]))
        self.input_min_int.setText("{0:.0f}".format(default_pdict["min_int"]))
        self.input_max_int.setText("{0:.0f}".format(default_pdict["max_int"]))
        self.input_min_pos.setText("{0:.0f}".format(default_pdict["min_pos"]))
        self.input_max_pos.setText("{0:.0f}".format(default_pdict["max_pos"]))
        self.input_min_lw.setText("{0:.1f}".format(default_pdict["min_lw"]))
        self.input_max_lw.setText("{0:.0f}".format(default_pdict["max_lw"]))
    
    # Upon button click, start the fit method of the currently selected peak
    # instance
    def start_fit_button_clicked(self):
        selected_peak = peaks[self.which_peak.currentIndex()]
        # The self argument passes the current GUI instance of which the data
        # is an attribute
        selected_peak.fit(self) 
                                
    
    # Upon button click, open a dialog to save fitresults of the currently
    # selected peak
    def save_fit_button_clicked(self):
        # Open a dialog where the user can choose the file path and name
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getSaveFileName(self, "Save File", "", "Numpy Files (*.npy)")
        if file_path:
            # Check which peak is currently selected in the dropdown menu
            index = self.which_peak.currentIndex()
            selected_peak = peaks[index]
            # Save fiterrors and index of the currently selected peak
            with open(file_path, "wb") as file:
                pickle.dump([selected_peak.fitresults, selected_peak.pdict, index], file)
                
    # Upon button click, open a dialog to open fitresults of a previously
    # fitted peak.
    def import_fit_button_clicked(self):
        # Open a dialog where the user can choose the file path and name
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(window, 'Import File')
        
        if file_path:
            # Load the fitresults and peak index
            with open(file_path, "rb") as file:
                fitresults, pdict, index = pickle.load(file)
                
            # Set loaded peak to currently selected and give it its fitresults
            selected_peak = peaks[index]
            selected_peak.fitresults = fitresults
            selected_peak.pdict = pdict
            self.which_peak.setCurrentIndex(index)
            
            # Display initial and threshold values used for the fit in the
            # text boxes, convert float to string
            self.input_init_int.setText("{0:.0f}".format(pdict["init_int"]))
            self.input_init_lw.setText("{0:.1f}".format(pdict["init_lw"]))
            self.input_min_int.setText("{0:.0f}".format(pdict["min_int"]))
            self.input_max_int.setText("{0:.0f}".format(pdict["max_int"]))
            self.input_min_pos.setText("{0:.0f}".format(pdict["min_pos"]))
            self.input_max_pos.setText("{0:.0f}".format(pdict["max_pos"]))
            self.input_min_lw.setText("{0:.1f}".format(pdict["min_lw"]))
            self.input_max_lw.setText("{0:.0f}".format(pdict["max_lw"]))
            
            # Display scan size in the text boxes
            self.input_sx_px.setText("{0:.0f}".format(pdict["size_x_px"]))
            self.input_sy_px.setText("{0:.0f}".format(pdict["size_y_px"]))
            self.input_sx_mu.setText("{0:.0f}".format(pdict["size_x_mu"]))
            self.input_sy_mu.setText("{0:.0f}".format(pdict["size_y_mu"]))
            self.input_min_mean_range.setText("{0:.0f}".format(pdict["min_mean_range"]))
            self.input_max_mean_range.setText("{0:.0f}".format(pdict["max_mean_range"]))
            
            # Display the maps in the canvasses
            selected_peak.plot_maps(self)

###############################################################################
""" A class for file imports """
###############################################################################
class importer:
    def __init__(self):
        
        # Only for debugging the fitting part of this code:
        self.xdata = np.load("xdata.npy")
        self.ydata = np.load("ydata.npy")
        
# =============================================================================
#         file_dialog = QFileDialog()
#         file_path, _ = file_dialog.getOpenFileName(window, 'Import File')
#         
#         if file_path:
#             # Get size of the scan in pixels from the text boxes
#             size_x_px = int(window.input_sx_px.text())
#             size_y_px = int(window.input_sy_px.text())
#             
#             # Get mean range from the text boxes
#             min_mean_range = int(window.input_min_mean_range.text())
#             max_mean_range = int(window.input_max_mean_range.text())
#             
#             # Call import function from module
#             self.xdata, self.ydata = mod.import_raman_scan(file_path, size_x_px, size_y_px, min_mean_range, max_mean_range)
# =============================================================================
    
        

###############################################################################
""" A class to be called for each new peak that is fitted """
###############################################################################
class peak:
    
    # The peak instance is passed the peak's pdict when creating it to extract
    # initial and threshold values
    def __init__(self, pdict):
        # Set pdict as attribute of the peak instance
        self.pdict = pdict
        # Initialize the fitresults attribute
        self.fitresults = 0
    
    # Fitting method: Retrieve all relevant user input from the GUI and call 
    # the fit function from the module. Results (fitresults and fiterrors) are
    # stored as attributes of the currently selected peak instance.
    def fit(self, window):
        
        # Get size of the scan and spectral range from the text boxes and add
        # them to the pdict 
        self.pdict["size_x_px"] = int(window.input_sx_px.text())
        self.pdict["size_y_px"] = int(window.input_sy_px.text())
        self.pdict["size_x_mu"] = float(window.input_sx_mu.text())
        self.pdict["size_y_mu"] = float(window.input_sy_mu.text())
        self.pdict["min_mean_range"] = float(window.input_min_mean_range.text())
        self.pdict["max_mean_range"] = float(window.input_max_mean_range.text())
        
        # Get initial and threshold values from the text boxes, write into pdict
        self.pdict["init_int"] = float(window.input_init_int.text())
        self.pdict["init_lw"] = float(window.input_init_lw.text())
        self.pdict["min_int"] = float(window.input_min_int.text())
        self.pdict["max_int"] = float(window.input_max_int.text())
        self.pdict["min_pos"] = float(window.input_min_pos.text())
        self.pdict["max_pos"] = float(window.input_max_pos.text())
        self.pdict["min_lw"] = float(window.input_min_lw.text())
        self.pdict["max_lw"] = float(window.input_max_lw.text())
        
        # Call the fit function
        self.fitresults = mod.fit_to_map(window.data.xdata, window.data.ydata, self.pdict)
        
        # When done, call the plot function
        self.plot_maps(window)
        
    # Plot the results    
    def plot_maps(self, window):
        
        # Create one instance of the ImageView widget for each plot
        map_int = pg.ImageView()
        map_pos = pg.ImageView()
        map_lw = pg.ImageView()
        
        # Connect them to the corresponding canvas in the GUI instance
        window.view_int.setCentralItem(map_int)
        window.view_pos.setCentralItem(map_pos)
        window.view_lw.setCentralItem(map_lw)
        
        # Display the image data
        window.view_int.setCentralItem(map_int)
        map_int.setImage(self.fitresults[:,:,1])
        map_pos.setImage(self.fitresults[:,:,2])
        map_lw.setImage(self.fitresults[:,:,3])
        
    
        
    
        
###############################################################################
""" Run the code """
###############################################################################
if __name__ == '__main__':
    
    # Create an instance of QtWidgets.QApplication
    app = QtWidgets.QApplication(sys.argv) 
    
    # Create an instance of the GUI
    window = GUI() 
    
    # A list of peak instances with one entry per peak (TA, G, LO, 2D)
    # Storing them in a list makes it easier to switch between peaks as the 
    # user selects them via the dropdown menu        
    peaks = [peak(default_values.TA()), peak(default_values.G()), peak(default_values.LO()), peak(default_values.TwoD())]
    
    # A list of dictionaries containing the default values of each peak
    default_pdicts = [default_values.TA(), default_values.G(), default_values.LO(), default_values.TwoD()]
    
    # Create an instance of the file import class
    #data = importer()
    
    # Start the application
    sys.exit(app.exec_())
