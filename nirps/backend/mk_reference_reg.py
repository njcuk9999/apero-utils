import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt


# From the Optical model, generate a DS9 .reg file to flag
# orders and wavelength solution.
#
# The file is assumed to be with the "nominal" orientation, having
# redder orders at low Y values and longer wavelengths within each
# order on the low X part of the order. Cross-dispersion is along the
# Y axis.

# only input that may eventually be updated
optical_model_file = 'NIRPS-4310-ULA-007_BE_F8_v01_all_orders.txt'

# read as an astropy table
optical_model=Table.read(optical_model_file,format = 'ascii')

# define orientation of orders on the array
pixscale = 0.015  # pix in mm on the science array. We will express everything in pixels, not microns
# pixel considered to be the optical center of the array
x_detector_center = 2047
y_detector_center = 2047


order = optical_model['ORDER']
wave = optical_model['WAVE']
xpix = -optical_model['Y']/pixscale+x_detector_center
ypix = -optical_model['X']/pixscale+y_detector_center


# header of the .reg file
reg = """# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
physical"""


# loop through orders
for uorder in np.unique(order):
    # find optical model points for that order
    g = uorder == order

    # add some nice colors
    color = 'red'
    if wave[g][0] < 1.4:
        color = 'green'
    if wave[g][0] < 1.14:
        color = 'blue'

    # adding the line along orders
    for i in range(len(xpix[g])-1):
        reg+='\n'
        reg+='line({0},{1},{2},{3}) # color={4} line=0 0 edit=0 move=0 rotate=0 delete=0'.format(int(xpix[g][i]),int(ypix[g][i]),int(xpix[g][i+1]),int(ypix[g][i+1]),color)

    # adding some text
    for i in range(0,len(g[g]),len(g[g])//5):
        reg+='\n# text('+str(xpix[g][i])+','+str(ypix[g][i])+') text={'+str(int(wave[g][i]*1000))+' nm, order '+str(int(order[g][i]))+'} edit=0 move=0 rotate=0 delete=0 background'

# save to a file
f = open('nominal_order_disposition.reg', "w")
f.write(reg)
f.close()