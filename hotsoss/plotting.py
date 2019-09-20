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
from bokeh.models import HoverTool, LogColorMapper, LogTicker, LinearColorMapper, ColorBar, Span, CustomJS, Slider
from bokeh.layouts import column
import numpy as np

from . import utils
from . import locate_trace as lt


# def plot_frames(data, idx=0, scale='linear', trace_coeffs=None, saturation=0.8, title=None, wavecal=None):
#     """
#     Plot a SOSS frame
#
#     Parameters
#     ----------
#     data: sequence
#         The 4D data to plot
#     scale: str
#         Plot scale, ['linear', 'log']
#     trace_coeffs: sequence
#         Plot the traces for the given coefficients
#     saturation: float
#         The full-well fraction to be considered saturated, (0, 1)
#     title: str
#         A title for the plot
#     wavecal: np.ndarray
#         A wavelength calibration map for each pixel
#     """
#     # Reshape into (frames, nrows, ncols)
#     dims = data.shape
#
#     # Test
#     data[1, 50:100, 50:100] *= 100.
#
#     # Get data, snr, and saturation for plotting
#     vmin = int(np.nanmin(data[data >= 0]))
#     vmax = int(np.nanmax(data[data < np.inf]))
#     dat = data
#     snr = np.sqrt(data)
#     fullWell = 65536.0
#     sat = dat > saturation * fullWell
#     sat = sat.astype(int)
#
#     # Wrap the data in two ColumnDataSources
#     source_visible = dict(data=[dat[idx]], snr=[snr[idx]], saturation=[sat[idx]])
#     source_available = dict(data=[dat], snr=[snr], saturation=[sat])
#
#     # Set the tooltips
#     tooltips = [("(x,y)", "($x{int}, $y{int})"), ("ADU/s", "@data"), ("SNR", "@snr"), ('Saturation', '@saturation')]
#
#     # Add wavelength calibration if possible
#     if isinstance(wavecal, np.ndarray):
#         if wavecal.shape == dat[idx].shape:
#             source_visible['wave1'] = [wavecal]
#             tooltips.append(("Wavelength", "@wave1"))
#         if wavecal.ndim == 3 and wavecal.shape[0] == 3:
#             source_visible['wave1'] = [wavecal[0]]
#             source_visible['wave2'] = [wavecal[1]]
#             source_visible['wave3'] = [wavecal[2]]
#             tooltips.append(("Wave 1", "@wave1"))
#             tooltips.append(("Wave 2", "@wave2"))
#             tooltips.append(("Wave 3", "@wave3"))
#
#     # Make the figure
#     fig = figure(x_range=(0, dims[2]), y_range=(0, dims[1]),
#                  tooltips=tooltips, width=1024, height=int(dims[1]/2.)+50,
#                  title=title, toolbar_location='above', toolbar_sticky=True)
#
#     # Plot the frame
#     if scale == 'log':
#         source['data'][source['data'] < 1.] = 1.
#         color_mapper = LogColorMapper(palette="Viridis256", low=vmin, high=vmax)
#         fig.image(source=source_visible, image='data', x=0, y=0, dw=dims[2], dh=dims[1], color_mapper=color_mapper)
#         color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(), orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))
#
#     else:
#         color_mapper = LinearColorMapper(palette="Viridis256", low=vmin, high=vmax)
#         fig.image(source=source_visible, image='data', x=0, y=0, dw=dims[2], dh=dims[1], palette='Viridis256')
#         color_bar = ColorBar(color_mapper=color_mapper, orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))
#
#     # Plot the trace polynomials
#     if trace_coeffs is not None:
#         X = np.linspace(0, 2048, 2048)
#
#         for coeffs in trace_coeffs:
#             Y = np.polyval(coeffs, X)
#             fig.line(X, Y, color='red')
#
#     slider = Slider(title='Frame', value=idx, start=0, end=dims[0], step=1)
#
#     # Define CustomJS callback, which updates the plot based on selected function
#     # by updating the source_visible ColumnDataSource
#     slider.callback = CustomJS(
#         args=dict(source_visible=source_visible,
#                   source_available=source_available),
#                   code="""
#                        var selected_function = cb_obj.value;
#
#                        // Get the data from the data sources
#                        var data_visible = source_visible.data;
#                        var data_available = source_available.data;
#
#                        // Change y-axis data according to the selected value
#                        data_visible = data_available[selected_function];
#
#                        // Update the plot
#                        source_visible.change.emit();
#                        """)
#
#     return column(fig, slider)


def plot_frame(frame, scale='linear', trace_coeffs=None, saturation=0.8, title=None, wavecal=None):
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
    saturation: float
        The full-well fraction to be considered saturated, (0, 1)
    title: str
        A title for the plot
    wavecal: np.ndarray
        A wavelength calibration map for each pixel
    """
    # Determine subarray
    nrows, ncols = frame.shape

    # Get data, snr, and saturation for plotting
    vmin = int(np.nanmin(frame[frame >= 0]))
    vmax = int(np.nanmax(frame[frame < np.inf]))
    dat = frame
    snr = np.sqrt(frame)
    fullWell = 65536.0
    sat = dat > saturation * fullWell
    sat = sat.astype(int)

    # Set the source data
    source = dict(data=[dat], snr=[snr], saturation=[sat])

    # Set the tooltips
    tooltips = [("(x,y)", "($x{int}, $y{int})"), ("ADU/s", "@data"), ("SNR", "@snr"), ('Saturation', '@saturation')]

    # Add wavelength calibration if possible
    if isinstance(wavecal, np.ndarray):
        if wavecal.shape == frame.shape:
            source['wave1'] = [wavecal]
            tooltips.append(("Wavelength", "@wave1"))
        if wavecal.ndim == 3 and wavecal.shape[0] == 3:
            source['wave1'] = [wavecal[0]]
            source['wave2'] = [wavecal[1]]
            source['wave3'] = [wavecal[2]]
            tooltips.append(("Wave 1", "@wave1"))
            tooltips.append(("Wave 2", "@wave2"))
            tooltips.append(("Wave 3", "@wave3"))

    # Make the figure
    fig = figure(x_range=(0, dat.shape[1]), y_range=(0, dat.shape[0]),
                 tooltips=tooltips, width=1024, height=int(nrows/2.)+50,
                 title=title, toolbar_location='above', toolbar_sticky=True)

    # Plot the frame
    if scale == 'log':
        source['data'][0][source['data'][0] < 1.] = 1.
        color_mapper = LogColorMapper(palette="Viridis256", low=vmin, high=vmax)
        fig.image(source=source, image='data', x=0, y=0, dw=dat.shape[1], dh=dat.shape[0], color_mapper=color_mapper)
        color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(), orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))

    else:
        color_mapper = LinearColorMapper(palette="Viridis256", low=vmin, high=vmax)
        fig.image(source=source, image='data', x=0, y=0, dw=dat.shape[1], dh=dat.shape[0], palette='Viridis256')
        color_bar = ColorBar(color_mapper=color_mapper, orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))

    # Plot the trace polynomials
    if trace_coeffs is not None:
        X = np.linspace(0, 2048, 2048)

        for coeffs in trace_coeffs:
            Y = np.polyval(coeffs, X)
            fig.line(X, Y, color='red')

    return fig


def plot_slice(frame, col, **kwargs):
    """
    Plot a column of a frame to see the PSF in the cross dispersion direction

    Parameters
    ----------
    data: np.ndarray
        The datacube
    col: int, sequence
        The column index(es) to plot
    """
    # Transpose data
    slc = frame.T

    # Turn one column into a list
    if isinstance(col, int):
        col = [col]

    # Plot the frame
    dfig = plot_frame(frame, **kwargs)

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


# def plot_lightcurve(self, column, time_unit='s', resolution_mult=20, draw=True):
#     """
#     Plot a lightcurve for each column index given
#
#     Parameters
#     ----------
#     column: int, float, sequence
#         The integer column index(es) or float wavelength(s) in microns
#         to plot as a light curve
#     time_unit: string
#         The string indicator for the units that the self.time array is in
#         ['s', 'min', 'h', 'd' (default)]
#     resolution_mult: int
#         The number of theoretical points to plot for each data point
#     draw: bool
#         Render the figure instead of returning it
#     """
#     # Check time_units
#     if time_unit not in ['s', 'min', 'h', 'd']:
#         raise ValueError("time_unit must be 's', 'min', 'h' or 'd'")
#
#     # Get the scaled flux in each column for the last group in
#     # each integration
#     flux_cols = np.nansum(self.datacube.reshape(self.dims3)[self.ngrps-1::self.ngrps], axis=1)
#     flux_cols = flux_cols/np.nanmax(flux_cols, axis=1)[:, None]
#
#     # Make it into an array
#     if isinstance(column, (int, float)):
#         column = [column]
#
#     # Make the figure
#     lc = figure()
#
#     for kcol, col in enumerate(column):
#
#         color = next(utils.COLORS)
#
#         # If it is an index
#         if isinstance(col, int):
#             lightcurve = flux_cols[:, col]
#             label = 'Column {}'.format(col)
#
#         # Or assumed to be a wavelength in microns
#         elif isinstance(col, float):
#             waves = np.mean(self.wave[0], axis=0)
#             lightcurve = [np.interp(col, waves, flux_col) for flux_col in flux_cols]
#             label = '{} um'.format(col)
#
#         else:
#             print('Please enter an index, astropy quantity, or array thereof.')
#             return
#
#         # Plot the theoretical light curve
#         if str(type(self.tmodel)) == "<class 'batman.transitmodel.TransitModel'>":
#
#             # Make time axis and convert to desired units
#             time = np.linspace(min(self.time), max(self.time), self.ngrps*self.nints*resolution_mult)
#             time = time*q.d.to(time_unit)
#
#             tmodel = batman.TransitModel(self.tmodel, time)
#             tmodel.rp = self.rp[col]
#             theory = tmodel.light_curve(tmodel)
#             theory *= max(lightcurve)/max(theory)
#
#             lc.line(time, theory, legend=label+' model', color=color, alpha=0.1)
#
#         # Convert datetime
#         data_time = self.time[self.ngrps-1::self.ngrps].copy()
#         data_time*q.d.to(time_unit)
#
#         # Plot the lightcurve
#         lc.circle(data_time, lightcurve, legend=label, color=color)
#
#     lc.xaxis.axis_label = 'Time [{}]'.format(time_unit)
#     lc.yaxis.axis_label = 'Transit Depth'
#
#     if draw:
#         show(lc)
#     else:
#         return lc


def plot_spectrum(wavelength, flux, fig=None, scale='log', legend=None, ylabel='Flux Density', xlabel='Wavelength [um]', width=1024, height=500, **kwargs):
    """
    Plot a generic spectrum

    Parameters
    ----------
    wavelength: array-like
        The 1D wavelength array
    flux: array-like
        The 1D counts or flux
    fig: bokeh.plotting.figure.Figure (optional)
        The plot to add the spectrum to 
    scale: str
        Plot scale, ['linear', 'log']
    legend: str
        The text for the legend
    ylabel: str
        The text for the y-axis
    xlabel: str
        The text for the x-axis
    width: int
        The width of the plot
    height: int
        The height of the plot

    Returns
    -------
    bokeh.plotting.figure.Figure
        The figure
    """
    # Make the figure
    if fig is None:
        fig = figure(x_axis_type=scale, y_axis_type=scale, width=width, height=height)

    # Add the spectrum plot
    fig.step(wavelength, flux, mode='center', legend=legend, **kwargs)
    fig.yaxis.axis_label = ylabel
    fig.xaxis.axis_label = xlabel

    return fig
