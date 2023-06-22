# twist-angle-analysis
This tool analyses spatially resolved Raman measurements recorded on twisted graphene heterostructures. It maps the properties of the TA, G, LO and 2D peak. From the TA peak position, it also creates a map of the local twist angle.

PREPARATION

    Export Raman maps from ProjectFive as .txt by choosing "Table" in the
    export menu. Make sure to set the spectral unit to "rel. 1/cm".


HOW TO USE THIS CODE

    As you scroll down the code, you will see that it has four sections:
        1. User input section
        2. Advanced user input section
        3. Don't touch section A (functions)
        4. Don't touch section B (code execution)
        
    Obviously, the last two sections are not meant to be edited as they contain
    the main body of this code (kind of the back end).
    As a user, here is what you have to do to analyse your Raman Scan:
    
    1. User input section
        - Provide the required information about the scan (see comments)
        - Choose which peaks to fit and to map
        - Run this code to start the fitting and mapping, or proceed to the

    2. Advanced user input section
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
