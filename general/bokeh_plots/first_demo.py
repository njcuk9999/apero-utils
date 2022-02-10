#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-02-08

@author: cook
"""
import numpy as np
import os

from bokeh.themes import built_in_themes
from bokeh.io import curdoc, show, output_file
from bokeh.plotting import figure
from bokeh.models import HoverTool
from bokeh.layouts import column, row, grid
from bokeh.models import ColumnDataSource, Slider, TextInput


# =============================================================================
# Define variables
# =============================================================================
workspace = '/data/spirou/misc/bokeh/'


testtooltip = """
<div>
    <div>
        <span style="font-size: 10px; color: #696;">x: @x</span><br>
        <span style="font-size: 10px; color: #696;">y: @y</span>
    </div>
</div>
"""

TEST = False


# =============================================================================
# Define functions
# =============================================================================
class App:
    def __init__(self, fig, func, create=True):

        self.fig = fig
        self.func = func
        # set starting values of options
        self.title_value = 'Current title'
        self.offset_value = 0.0
        self.amplitude_value = 1.0
        self.phase_value = 0.0
        self.freq_value = 1.0
        self.npoints_value = 100
        # set up the data source
        self.source = ColumnDataSource()

        # set up widgets
        if create:
            self.title = TextInput(title='title', value=self.title_value)
            self.offset = Slider(title='offset', value=self.offset_value,
                                 start=-5.0, end=5.0, step=0.1)
            self.amplitude = Slider(title='amplitude', value=self.amplitude_value,
                                    start=-5.0, end=5.0, step=0.1)
            self.phase = Slider(title='phase', value=self.phase_value,
                                start=0.0, end=2*np.pi)
            self.freq = Slider(title='frequency', value=self.freq_value,
                               start=0.1, end=5.1, step=0.1)
            self.npoints = Slider(title='Npoints', value=self.npoints_value,
                                  start=0, end=500, step=10)
            # update data source using starting values of options
            self.calculate()
            # define those inputs that use update_data function
            self.update_data_inputs = [self.offset, self.amplitude, self.phase,
                                       self.freq, self.npoints]
            # link title to change
            self.title.on_change('value', self.update_title)
            # link them to the updater
            for w in self.update_data_inputs:
                w.on_change('value', self.update_data)
            # create list of inputs
            self.inputs = [self.title] + self.update_data_inputs

    def calculate(self):
        # Get the current slider values
        a = self.amplitude.value
        b = self.offset.value
        w = self.phase.value
        k = self.freq.value
        n = self.npoints_value

        # deal with low values of n
        if n < 1:
            n = 1
        # Generate the new curve
        namespace = dict(a=a, b=b, w=w, k=k, n=n)
        x, y = self.func(**namespace)
        self.source.data = dict(x=x, y=y)

    def update_data(self, attrname, old, new):

        self.amplitude_value = self.amplitude.value
        self.offset_value = self.offset.value
        self.phase_value = self.phase.value
        self.freq_value = self.freq.value
        self.npoints_value = self.npoints.value
        self.calculate()

    def update_title(self, attrname, old, new):
        self.fig.title.text = self.title.value

    def newcopy(self, fig, func) -> 'App':
        new = App(fig, func=func, create=False)
        # set up the data source
        new.source = ColumnDataSource()
        # current values
        new.title_value = self.title_value
        new.offset_value = self.offset_value
        new.amplitude_value = self.amplitude_value
        new.phase_value = self.phase_value
        new.freq_value = self.freq_value
        new.npoints_value = self.npoints_value
        # set up widgets
        new.title = self.title
        new.offset = self.offset
        new.amplitude = self.amplitude
        new.phase = self.phase
        new.freq = self.freq
        new.npoints = self.npoints
        # define those inputs that use update_data function
        new.update_data_inputs = [new.offset, new.amplitude, new.phase,
                                   new.freq, new.npoints]
        # link them to the updater
        for w in new.update_data_inputs:
            w.on_change('value', new.update_data)
        # update data source using starting values of options
        new.calculate()
        return new

# =============================================================================
# Start of code
# =============================================================================


# set up our plotting window
plot1 = figure(height=400, width=800,
               tools="crosshair,pan,reset,save,wheel_zoom",
               x_range=[0, 4*np.pi], y_range=[-2.5, 2.5])
plot2 = figure(height=400, width=800, title='dy/dx',
               tools="crosshair,pan,reset,save,wheel_zoom",
               x_range=[0, 4*np.pi], y_range=[-2.5, 2.5])

def func1(a, k, w, b, n):
    x = np.linspace(0, 4 * np.pi, n)
    return x, a * np.sin(k * x + w) + b

def func2(a, k, w, b, n):
    x = np.linspace(0, 4 * np.pi, n)
    return x, a *k * np.cos(k * x + w)

# create the data source
app = App(plot1, func=func1)
app2 = app.newcopy(fig=plot2, func=func2)

# plot scatter plot
plot1.circle(source=app.source, x='x', y='y', color='red',
             legend_label='a*sin(k*x + w) + b')
plot2.circle(source=app2.source, x='x', y='y', color='blue',
             legend_label='a*k*cos(k*x + w)')

plot1.title.text = app.title_value

# add x label
plot1.xaxis.axis_label = 'X value'
plot1.yaxis.axis_label = 'Y value'
plot2.xaxis.axis_label = 'X value'
plot2.yaxis.axis_label = 'Y value'

# construct hover tool
hover = HoverTool(tooltips=testtooltip)
plot1.add_tools(hover)

inputs = column(*app.inputs)

if TEST:
    output_file(os.path.join(workspace, 'first_demo.html'))
    show(row(inputs, plot1, width=1200))
else:
    curdoc().add_root(grid([inputs, [plot1, plot2]]))
    curdoc().title = 'first_test'
    curdoc().theme = 'night_sky'

# =============================================================================
# End of code
# =============================================================================
