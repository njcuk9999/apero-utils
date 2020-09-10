import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline


def bisector(rv, ccf, doplot = False, low_high_cut = 0.1):
    # use the props from the CCF determination code
    # Could be per-order or with the mean
    #rv = props['RV_CCF']
    #ccf = props['MEAN_CCF']

    # get minima
    imin = np.argmin(ccf)

    # get point where the derivative changes sign at the edge of the line
    # the bisector is ambiguous passed this poind


    width_blue =  imin - np.max(np.where(np.gradient(ccf[0:imin])>0))
    width_red = np.min(np.where(np.gradient(ccf[imin:])<0))

    # get the width from the side of the center that reaches
    # that point first
    width = np.min([width_blue, width_red])

    # set depth to zero
    ccf -= np.min(ccf)

    # set continuum to one
    ccf /= np.min( ccf[ [imin - width,imin + width] ])

    # interpolate each side of the ccf slope at a range of depths
    depth = np.arange(0.1,0.91,0.01)

    # blue and red side of line
    g1 = (ccf[imin:imin - width:-1]>low_high_cut) & (ccf[imin:imin - width:-1]<(1-low_high_cut))
    spline1 = InterpolatedUnivariateSpline(ccf[imin:imin - width:-1][g1],rv[imin:imin - width:-1 ][g1], k=2)

    g2 = (ccf[imin : imin + width]>low_high_cut) & (ccf[imin : imin + width]<(1-low_high_cut))
    spline2 = InterpolatedUnivariateSpline(ccf[imin : imin + width][g2],rv[imin : imin + width][g2], k=2)

    # get midpoint
    bisector = (spline2(depth)+spline1(depth))/2

    # get bisector widht
    width = (spline2(depth)-spline1(depth))

    if doplot:
        # some nice plots
        plt.plot(rv[imin - width : imin+ width],ccf[imin - width : imin+ width],label = 'ccf')
        plt.plot(bisector,depth,label = 'bisector')
        plt.plot((bisector-np.mean(bisector))*100+np.mean(bisector),depth, label = 'bisector * 100')
        plt.legend()
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('Depth')
        plt.savefig('ccf_bisector.png')
        plt.show()

    return depth, bisector, width
