# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:39:42 2020

@author: Schäpers
"""

import Library as lib 
import matplotlib as mpl
import matplotlib.pyplot as plt

###############################################################################
"""Read spectrum data and customize the plot layout"""
###############################################################################
folder = "//serveri2a/Transfer/Aaron/Work on tBLG/TA Paper Workshop"
filename = "Examplary spectrum 7deg.txt"
xlist, ylist = lib.read_single_spectrum(filename, folder)

plotsize = 10,7
save = True
title = filename


###############################################################################
"""Plot a single spectrum"""
###############################################################################
#lib.plot_spectrum(xlist, ylist, plotsize, save, filename, folder)




###############################################################################
"""Plot several spectra in one graph"""
###############################################################################
folder = "//serveri2a/Transfer/Aaron/Work on tBLG/TA Paper Workshop"
filename1 = "Examplary spectrum 7deg.txt"
filename2 = "Examplary spectrum 8deg.txt"
title = "Spectra 7+8"

xlist1, ylist1 = np.array(lib.read_single_spectrum(filename1, folder))
xlist2, ylist2 = np.array(lib.read_single_spectrum(filename2, folder))

#Normalize the spectra
ylist1 = ylist1/max(ylist1)
ylist2 = ylist2/max(ylist2)



fig = plt.figure(figsize=plotsize)
ax1 = fig.add_subplot(211)
ax1.set_yticks([])
ax1.set_ylim(0.3,0.6)
ax1.set_xlabel(r"Raman shift (cm$^{-1})$")
ax1.set_ylabel("Intensity (a.u.)")
ax1.plot(xlist1, ylist1, color="blue", label="7deg")

ax2 = ax1.twinx()
ax2.set_yticks([])
ax2.set_ylim(0,0.6)
ax2.plot(xlist2, ylist2-0.1, color="red", label="8deg")


plt.legend()
plt.tight_layout()

#plt.savefig(str(folder)+"/"+str(title)+".png", format="png", dpi=900)
