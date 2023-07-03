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
    Connect the peak drop down menu
    Make default values appear in the GUI based on which peak was chosen
    
At all times, there is one peak instance for TA, G, LO and 2D.
The peak selected in the dropdown menu needs to act as a kind of condition that 
defines on which peak instance all the other actions (set default values, fit, map)
are performed.  
"""

import sys
import numpy as np

from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QLineEdit, QPushButton, QComboBox, QMainWindow
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

import default_values


###############################################################################
""" The main window """
###############################################################################
# Create a base class that will load the .ui file in the constructor
# It inherits the QtWidgets.QMainWindow class because I created a new 
# "Main Window" when selecting the form type when first creating the .ui file 
# in PyQt Designer
class GUI(QMainWindow):
    
    # Define instances of pyqtSignal for all buttons that emit a signal
    import_scan_button_clicked = pyqtSignal()

    def __init__(self):
        # Call the inherited classes __init__ method
        super(GUI, self).__init__() 
        
        # Load the .ui file
        uic.loadUi('gui.ui', self) 
        
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

        # Make the "Import Raman scan" button emit a signal 
        self.import_scan_button = self.findChild(QPushButton, 'b_import_scan')
        self.import_scan_button.clicked.connect(self.import_scan_button_clicked.emit)
        
        # Connect peak dropdown menu to an attribute of this class
        self.which_peak = self.findChild(QComboBox, 'in_peak')
        
        # Make the "Use default values" button emit a signal
        self.default_values_button = self.findChild(QPushButton, 'b_default_values')
        self.default_values_button.clicked.connect(self.set_default_values)
                
        # Show the GUI
        self.show()
    
    # Upon button click, set default values in the lineEdits
    #@pyqtSlot()
    def set_default_values(self):
        selected_peak = peaks[self.which_peak.currentIndex()]
        # Convert float to string
        self.input_init_int.setText(str(selected_peak.init_int))
        self.input_init_lw.setText(str(selected_peak.init_lw))
        
        self.input_min_int.setText(str(selected_peak.min_int))
        self.input_max_int.setText(str(selected_peak.max_int))
        
        self.input_min_pos.setText(str(selected_peak.min_pos))
        self.input_max_pos.setText(str(selected_peak.max_pos))
        
        self.input_min_lw.setText(str(selected_peak.min_lw))
        self.input_max_lw.setText(str(selected_peak.max_lw))
        
    
        


###############################################################################
""" A class for file imports """
###############################################################################
class file_import(QObject):
    def __init__(self):
        super().__init__()
        
    # Open a file dialog to choose a .txt file containing a Raman scan
    # Then read the Raman scan data and set the baseline to zero
    @pyqtSlot()
    def scan(self):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(window, 'Import File')
        
        if file_path:
            # Get size of the scan in pixels
            nx = int(window.input_sx_px.text())
            ny = int(window.input_sy_px.text())
            
            # Get mean range
            min_mean_range = int(window.input_min_mean_range.text())
            max_mean_range = int(window.input_max_mean_range.text())
            
            # This is a tuple of strings with 1600 entries, i.e., one entry per 
            # CCD pixel, i.e., one entry per data point in the Raman spectrum.
            # Each of the 1600 strings is a succession of numbers giving the CCD counts
            # on the corresponding pixel for all the spatial data points measured in the scan.
            lines = tuple(open(file_path, 'r')) 
            
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
            self.data = np.zeros((ny,nx,len(data_str)))
            
            # Get indices of the averaging interval for baseline shifting
            i_xmin_mean = np.abs(xlist - min_mean_range).argmin()
            i_xmax_mean = np.abs(xlist - max_mean_range).argmin()
            
            # Scan over all spatial data points
            x,y = 0,0
            for i in range(nx*ny):
                # Calculate the average CCD counts in a range without peaks
                ymean = np.mean(data_float[i, i_xmin_mean:i_xmax_mean])
                # Subtract that average from the entire spectrum and stack the spectrum
                # into the data array
                self.data[y,x] = [y-ymean for y in data_float[i]]          

                x = x+1
                if x == nx:
                    y = y+1
                    x = 0
                    
            # Convert xlist to array for more efficient function evaluation
            self.xdata = np.asarray(xlist)


###############################################################################
""" A class to be called for each new peak that is fitted """
###############################################################################
class peak:
    
    # The peak instance is passed the peak's pdict when creating it
    def __init__(self, pdict):
        # Set initial and threshold values as attributes of the peak instance
        self.init_int = pdict["init_values"][1]
        self.init_lw = pdict["init_values"][3]
        self.min_int = pdict["min_values"][0]
        self.max_int = pdict["max_values"][0]
        self.min_pos = pdict["min_values"][1]
        self.max_pos = pdict["max_values"][1]
        self.min_lw = pdict["min_values"][2]
        self.max_lw = pdict["max_values"][2]
    
    #def use_default_values(self, index):
        
        
        
###############################################################################
""" Run the code """
###############################################################################
if __name__ == '__main__':
    
    # Create an instance of QtWidgets.QApplication
    app = QtWidgets.QApplication(sys.argv) 
    
    # Create an instance of the GUI
    window = GUI() 
    
    # Introduce a list of peak class instances,
    # one entry per peak (TA, G, LO, 2D)
    # Storing them in a list makes it easier to switch between peaks as the 
    # user selects them via the dropdown menu        
    peaks = [peak(default_values.TA()), peak(default_values.G()), peak(default_values.LO()), peak(default_values.TwoD())]
    
    # Create an instance of the file import class
    file_importer = file_import()
    
    # Connect signals emitted by buttons to their respective methods in other
    # classes
    window.import_scan_button_clicked.connect(file_importer.scan)
    
    # Start the application
    sys.exit(app.exec_())
