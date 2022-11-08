#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tester code for selecting points on a graph
For use with default reference wave solution

or maybe for localisation (selecting limits and number of orders)

Created on 2022-11-08 at 12:05

@author: cook
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import warnings


# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
# Simple mouse click function to store coordinates
class Click():
    def __init__(self, ax, func, button=1):
        self.ax = ax
        self.func = func
        self.button = button
        self.press = False
        self.move = False
        self.c1 = self.ax.figure.canvas.mpl_connect('button_press_event',
                                                    self.onpress)
        self.c2 = self.ax.figure.canvas.mpl_connect('button_release_event',
                                                    self.onrelease)
        self.c3 = self.ax.figure.canvas.mpl_connect('motion_notify_event',
                                                    self.onmove)

    def onclick(self, event):
        if event.inaxes == self.ax:
            if event.button == self.button:
                self.func(event, self.ax)

    def onpress(self, event):
        self.press = True

    def onmove(self, event):
        if self.press:
            self.move = True

    def onrelease(self, event):
        if self.press and not self.move:
            self.onclick(event)
        self.press = False
        self.move = False



class ZoomWindow():
    def __init__(self, ax, func, plot_args, scale=0.1):
        self.ax = ax
        self.func = func
        self.plot_args = plot_args
        self.scale = scale
        self.c1 = self.ax.figure.canvas.mpl_connect('motion_notify_event',
                                                    self.onmove)
        self.limits = self.ax.axis()

    def onmove(self, event):
        self.func(self.ax, None, self.plot_args)
        self.zoomax = self.ax.inset_axes([0.75, 0.75, 0.22, 0.22])
        self.zoomax.set_xticklabels([])
        self.zoomax.set_yticklabels([])
        self.func(self.zoomax, event, self.plot_args, self.scale,
                  self.limits)

        # self.ax.indicate_inset_zoom(self.zoomax, edgecolor='black')
        self.ax.figure.canvas.draw()


class KeyAction():
    def __init__(self, ax):
        self.c1 = ax.figure.canvas.mpl_connect('key_press_event',
                                                    self.onpress)

    def onpress(self, event):

        global newcoords, oldcoords, count

        if event.key == 'r':
            newcoords = []
            oldcoords = []
            count = 0
        if event.key == 'z':
            if count % 2 != 0:
                newcoords = newcoords[:-1]
            else:
                oldcoords = oldcoords[:-1]
            count -= 1

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------

    import matplotlib.pyplot as plt





    newcoords = []
    oldcoords = []
    count = 0


    # plot code here
    def click_action(event, ax):
        global newcoords, oldcoords, count

        if count % 2 == 0:
            newcoords.append((event.xdata, event.ydata))
            print('NEW', event.xdata, event.ydata)
        else:
            oldcoords.append((event.xdata, event.ydata))
            print('OLD', event.xdata, event.ydata)
        # update count
        count += 1

    fig, frame = plt.subplots(ncols=1, nrows=1)

    x = np.arange(1000.0)
    y = np.random.random_integers(0, 100, size=1000).astype(float)
    y[y < 80] = 0
    y = y * np.sin(x/100)
    y[y < 0] = 0.01 * abs(y[y < 0])

    def our_plot(ax, event, plot_args, scale=1, limits=None):
        _x = plot_args['x']
        _y = plot_args['y']
        ax.cla()
        ax.plot(_x, _y, zorder=0, color='k', alpha=0.5)
        ax.plot(_x + 200, _y, zorder=0, color='r', alpha=0.5)
        if hasattr(event, 'xdata') and event is not None:
            xmin, xmax, ymin, ymax = limits
            deltax, deltay = xmax - xmin, ymax - ymin

            if event.xdata is not None and event.ydata is not None:
                print(event.xdata, event.ydata,  xmin, xmax, ymin, ymax, scale)
                ax.set(xlim=[event.xdata-scale*deltax,
                             event.xdata+scale*deltax],
                       ylim=[event.ydata-scale*deltay,
                             event.ydata+scale*deltay])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.plot([event.xdata], [event.ydata], marker='+', color='k',
                        ms=20, alpha=0.25)
        else:
            ax.set(title='Press "r" to reset, "z" to undo last point')

        # plot coords
        for it in range(len(newcoords)):
            ax.scatter(*newcoords[it],
                       marker=r"$ {} $".format(it),
                       color='red', zorder=10)
        for it in range(len(oldcoords)):
            ax.scatter(*oldcoords[it],
                       marker=r"$ {} $".format(it),
                       color='black', zorder=10)


    plot_args = dict(x=x, y=y)
    our_plot(frame, None, plot_args)

    # deal with clicking
    click = Click(frame, click_action, button=1)
    zoom = ZoomWindow(frame, our_plot, plot_args, scale=0.05)
    action = KeyAction(frame)
    # then show
    plt.show()

# =============================================================================
# End of code
# =============================================================================
