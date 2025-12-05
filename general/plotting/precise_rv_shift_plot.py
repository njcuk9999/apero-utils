#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-06-12 at 08:52

@author: cook
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.ndimage import shift

# =============================================================================
# Define variables
# =============================================================================
N_PIXELS = 100
IMG_HEIGHT = 5
TRACE_WID = 10
PIXEL_SHIFT = 0.1  # Try 1, 0.1, 0.01 etc.
SHIFT_NAME = '1/10'
VELO_NAME = '100 m/s'
# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
class Shifter():
    def __init__(self, n_pixels=N_PIXELS, img_height=IMG_HEIGHT,
                 trace_wid = TRACE_WID, pixel_shift=PIXEL_SHIFT,
                 shift_name=SHIFT_NAME, velo_name=VELO_NAME):
        # Create x and base trace
        self.x = np.arange(n_pixels)
        self.center = n_pixels // 2
        self.trace = np.exp(-0.5 * ((self.x - self.center) / (trace_wid / 2)) ** 2)
        self.trace_shifted = shift(self.trace, pixel_shift, mode='nearest')

        # Tile vertically
        self.img_original = np.tile(self.trace, (img_height, 1))
        self.img_shifted = np.tile(self.trace_shifted, (img_height, 1))

        self.fig = None
        self.frames = None
        self.im = None
        self.line0 = None
        self.line1 = None
        self.text1 = None
        self.text2 = None

        self.shift_name = shift_name
        self.pixel_shift = pixel_shift
        self.velo_name = velo_name

    def plot(self):
        # Set up figure and axes
        self.fig, self.frames = plt.subplots(2, 1, figsize=(8, 8),
                                             gridspec_kw={'height_ratios': [3, 1]})

        # imshow plot
        shift_text = (f'Shift = {self.shift_name} px \n'
                      f'        = {self.velo_name}')

        self.im = self.frames[0].imshow(self.img_original,
                                   cmap='viridis', aspect='auto')
        self.text1 = self.frames[0].text(0.05, 0.95, shift_text,
                                         color='white', fontsize=24,
                                         ha='left', va='top',
                                         transform=self.frames[0].transAxes)
        self.text2 = self.frames[1].text(0.05, 0.9, 'Frame {0}'.format(1),
                                         color='black', fontsize=18,
                                         ha='left', va='top',
                                         transform=self.frames[1].transAxes)

        self.frames[0].set_xticks([])
        self.frames[0].set_yticks([])
        # line plot
        self.line0, = self.frames[1].plot(self.x, self.trace, color='orange',
                                     label='Unshifted')
        self.line1, = self.frames[1].plot(self.x, self.trace, color='blue',
                                     label=f'Shifted')
        self.frames[1].set_ylim(0, 1.1)
        self.frames[1].set_xticks([])
        self.frames[1].set_yticks([])
        self.frames[1].legend(loc=0, ncol=1, fontsize=18)

        return 0

    # Animation update function
    def update(self, frame):
        if frame % 2 == 0:
            self.im.set_data(self.img_original)
            self.line1.set_ydata(self.trace)
        else:
            self.im.set_data(self.img_shifted)
            self.line1.set_ydata(self.trace_shifted)

        self.text2.set_text(f'Frame = {1 + frame % 2}')

        return [self.im, self.line1, self.text1, self.text2]

    def run(self, interval: int = 1000):
        # Animate: blink once per second
        ani = FuncAnimation(self.fig, self.update, frames=100,
                            interval=interval, blit=True)
        plt.tight_layout()
        ani.save(f'Shift_{self.pixel_shift}.gif', writer='pillow', fps=1)
        plt.show()
        plt.close()


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    shifts = [10.0, 1.0, 1/10, 1/100, 1/1000]
    shift_names = ['10', '1', '1/10', '1/100', '1/1000']
    velo_names = ['10 km/s', '1 km/s', '100 m/s', '10 m/s', '1 m/s']

    for it in range(len(shifts)):

        print('Shifting {}'.format(shift_names[it]))

        plot = Shifter(pixel_shift=shifts[it], shift_name=shift_names[it],
                       velo_name=velo_names[it])
        plot.plot()
        plot.run()

        uinput = input('Next? [Y]es or [N]o: ')
        # get out
        if 'N' in uinput.upper():
            break




# =============================================================================
# End of code
# =============================================================================
