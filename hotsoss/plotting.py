# -*- coding: utf-8 -*-
"""
A module to plot SOSS data

Authors: Joe Filippazzo
"""
import copy
from pkg_resources import resource_filename
import os

from astropy.io import fits
from bokeh.plotting import figure, show
from bokeh.models import HoverTool, LogColorMapper, LogTicker, LinearColorMapper, ColorBar, Span
from bokeh.layouts import column
import numpy as np

from . import utils
from . import locate_trace as lt


def plot_frame(frame, scale='linear', trace_coeffs=None, saturation=0.8, title=None):
    """
    Plot a SOSS frame

    Parameters
    ----------
    frame: sequence
        The 2D frame to plot
    scale: str
        Plot scale, ['linear', 'log']
    trace_coeffs: sequence
        Plot the traces for the given coefficients
    draw: bool
        Render the figure instead of returning it
    """
    # Determine subarray
    nrows, ncols = frame.shape

    # Get data, snr, and saturation for plotting
    vmax = int(np.nanmax(frame[frame < np.inf]))
    dat = frame
    snr = np.sqrt(frame)
    fullWell = 65536.0
    sat = dat > saturation * fullWell
    sat = sat.astype(int)

    # Make the figure
    tooltips = [("(x,y)", "($x{int}, $y{int})"), ("ADU/s", "@data"), ("SNR", "@snr"), ('Saturation', '@saturation')]
    fig = figure(x_range=(0, dat.shape[1]), y_range=(0, dat.shape[0]),
                 tooltips=tooltips, width=1024, height=int(nrows/2.)+50,
                 title=title, toolbar_location='above', toolbar_sticky=True)

    # Plot the frame
    if scale == 'log':
        dat[dat < 1.] = 1.
        source = dict(data=[dat], snr=[snr], saturation=[sat])
        color_mapper = LogColorMapper(palette="Viridis256", low=dat.min(), high=dat.max())
        fig.image(source=source, image='data', x=0, y=0, dw=dat.shape[1],
                  dh=dat.shape[0], color_mapper=color_mapper)
        color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(),
                             orientation="horizontal", label_standoff=12,
                             border_line_color=None, location=(0, 0))

    else:
        source = dict(data=[dat], snr=[snr], saturation=[sat])
        color_mapper = LinearColorMapper(palette="Viridis256", low=dat.min(), high=dat.max())
        fig.image(source=source, image='data', x=0, y=0, dw=dat.shape[1],
                  dh=dat.shape[0], palette='Viridis256')
        color_bar = ColorBar(color_mapper=color_mapper,
                             orientation="horizontal", label_standoff=12,
                             border_line_color=None, location=(0, 0))

    # Plot the polynomial too
    if trace_coeffs is not None:
        X = np.linspace(0, 2048, 2048)

        # Order 1
        Y = np.polyval(trace_coeffs[0], X)
        fig.line(X, Y, color='red')

        # Order 2
        Y = np.polyval(trace_coeffs[1], X)
        fig.line(X, Y, color='red')

    return fig

def plot_slice(frame, col, idx=0, **kwargs):
    """
    Plot a column of a frame to see the PSF in the cross dispersion direction

    Parameters
    ----------
    data: np.ndarray
        The datacube
    col: int, sequence
        The column index(es) to plot
    idx: int
        The frame index to plot
    """
    # Transpose data
    slc = frame.T

    # Turn one column into a list
    if isinstance(col, int):
        col = [col]

    # Plot the frame
    dfig = plot_frame(frame)

    # Make the figure
    fig = figure(width=1024, height=500)
    fig.xaxis.axis_label = 'Row'
    fig.yaxis.axis_label = 'Count Rate [ADU/s]'
    fig.legend.click_policy = 'mute'
    for c in col:
        color = next(utils.COLORS)
        fig.line(np.arange(slc[c, :].size), slc[c, :], color=color, legend='Column {}'.format(c))
        vline = Span(location=c, dimension='height', line_color=color, line_width=3)
        dfig.add_layout(vline)

    return column(fig, dfig)

def plot_ramp(data):
    """
    Plot the total flux on each frame to display the ramp
    """
    fig = figure()
    dims3 = data.shape[0]*data.shape[1], data.shape[2], data.shape[3]
    x = range(dims3[0])
    y = np.sum(data.reshape(dims3), axis=(-1, -2))
    fig.circle(x, y, size=8)
    fig.xaxis.axis_label = 'Group'
    fig.yaxis.axis_label = 'Count Rate [ADU/s]'

    return fig

def plot_lightcurve(self, column, time_unit='s', resolution_mult=20, draw=True):
    """
    Plot a lightcurve for each column index given

    Parameters
    ----------
    column: int, float, sequence
        The integer column index(es) or float wavelength(s) in microns
        to plot as a light curve
    time_unit: string
        The string indicator for the units that the self.time array is in
        ['s', 'min', 'h', 'd' (default)]
    resolution_mult: int
        The number of theoretical points to plot for each data point
    draw: bool
        Render the figure instead of returning it
    """
    # Check time_units
    if time_unit not in ['s', 'min', 'h', 'd']:
        raise ValueError("time_unit must be 's', 'min', 'h' or 'd'")

    # Get the scaled flux in each column for the last group in
    # each integration
    flux_cols = np.nansum(self.datacube.reshape(self.dims3)[self.ngrps-1::self.ngrps], axis=1)
    flux_cols = flux_cols/np.nanmax(flux_cols, axis=1)[:, None]

    # Make it into an array
    if isinstance(column, (int, float)):
        column = [column]

    # Make the figure
    lc = figure()

    for kcol, col in enumerate(column):

        color = next(utils.COLORS)

        # If it is an index
        if isinstance(col, int):
            lightcurve = flux_cols[:, col]
            label = 'Column {}'.format(col)

        # Or assumed to be a wavelength in microns
        elif isinstance(col, float):
            waves = np.mean(self.wave[0], axis=0)
            lightcurve = [np.interp(col, waves, flux_col) for flux_col in flux_cols]
            label = '{} um'.format(col)

        else:
            print('Please enter an index, astropy quantity, or array thereof.')
            return

        # Plot the theoretical light curve
        if str(type(self.tmodel)) == "<class 'batman.transitmodel.TransitModel'>":

            # Make time axis and convert to desired units
            time = np.linspace(min(self.time), max(self.time), self.ngrps*self.nints*resolution_mult)
            time = time*q.d.to(time_unit)

            tmodel = batman.TransitModel(self.tmodel, time)
            tmodel.rp = self.rp[col]
            theory = tmodel.light_curve(tmodel)
            theory *= max(lightcurve)/max(theory)

            lc.line(time, theory, legend=label+' model', color=color, alpha=0.1)

        # Convert datetime
        data_time = self.time[self.ngrps-1::self.ngrps].copy()
        data_time*q.d.to(time_unit)

        # Plot the lightcurve
        lc.circle(data_time, lightcurve, legend=label, color=color)

    lc.xaxis.axis_label = 'Time [{}]'.format(time_unit)
    lc.yaxis.axis_label = 'Transit Depth'

    if draw:
        show(lc)
    else:
        return lc

def plot_spectrum(self, frame=0, order=None, noise=False, scale='log', draw=True):
    """
    Parameters
    ----------
    frame: int
        The frame number to plot
    order: sequence
        The order to isolate
    noise: bool
        Plot with the noise model
    scale: str
        Plot scale, ['linear', 'log']
    draw: bool
        Render the figure instead of returning it
    """
    if order in [1, 2]:
        tso = getattr(self, 'tso_order{}_ideal'.format(order))
    else:
        if noise:
            tso = self.tso
        else:
            tso = self.tso_ideal

    # Get extracted spectrum (Column sum for now)
    wave = np.mean(self.wave[0], axis=0)
    flux_out = np.sum(tso.reshape(self.dims3)[frame].data, axis=0)
    response = 1./self.order1_response

    # Convert response in [mJy/ADU/s] to [Flam/ADU/s] then invert so
    # that we can convert the flux at each wavelegth into [ADU/s]
    flux_out *= response/self.time[np.mod(self.ngrps, frame)]

    # Trim wacky extracted edges
    flux_out[0] = flux_out[-1] = np.nan

    # Plot it along with input spectrum
    flux_in = np.interp(wave, self.star[0], self.star[1])

    # Make the spectrum plot
    spec = figure(x_axis_type=scale, y_axis_type=scale, width=1024, height=500)
    spec.step(wave, flux_out, mode='center', legend='Extracted', color='red')
    spec.step(wave, flux_in, mode='center', legend='Injected', alpha=0.5)
    spec.yaxis.axis_label = 'Flux Density [{}]'.format(self.star[1].unit)

    # Get the residuals
    res = figure(x_axis_type=scale, x_range=spec.x_range, width=1024, height=150)
    res.step(wave, flux_out-flux_in, mode='center')
    res.xaxis.axis_label = 'Wavelength [{}]'.format(self.star[0].unit)
    res.yaxis.axis_label = 'Residuals'

    if draw:
        show(column(spec, res))
    else:
        return column(spec, res)
