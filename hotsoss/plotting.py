# -*- coding: utf-8 -*-
"""
A module to plot SOSS data

Authors: Joe Filippazzo
"""
import copy
from pkg_resources import resource_filename
import os

from astropy.io import fits
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool, LogColorMapper, FixedTicker, BasicTickFormatter, FuncTickFormatter, BasicTicker, LogTicker, LinearColorMapper, ColorBar, Span, CustomJS, Slider
from bokeh.models.widgets import Panel, Tabs
from bokeh.layouts import gridplot, column
import numpy as np

from . import utils
from . import locate_trace as lt


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
    dat = frame
    snr = np.sqrt(frame)
    fullWell = 65536.0
    sat = dat > saturation * fullWell
    sat = sat.astype(int)
    dh, dw = dat.shape

    # Fix if log scale
    if scale == 'log':
        dat[dat < 1.] = 1.

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

    # Set shared plot params
    x_range = (0, dat.shape[1])
    y_range = (0, dat.shape[0])
    height = int(nrows/2.)+160
    toolbar = 'above'

    # Draw the figures
    tabs = []
    for pname, ptype in zip(['Counts', 'SNR', 'Saturation ({}% Full Well)'.format(saturation*100)], ['data', 'snr', 'saturation']):

        # Make the figure
        fig_title = '{} - {}'.format(title, pname)
        fig = figure(x_range=x_range, y_range=y_range, tooltips=tooltips, width=1024, height=height, title=fig_title, toolbar_location=toolbar, toolbar_sticky=True)

        # Get the data
        vals = source[ptype][0]

        # Saturation plot is different
        if ptype == 'saturation':
            vmin = 0
            vmax = 1
            formatter = FuncTickFormatter(code="""return {0: 'Unsaturated', 1: 'Saturated'}[tick]""")
            color_map = ['#404387', '#FDE724']
            ticker = FixedTicker(ticks=[vmin, vmax])

        # Counts and SNR are similar plots
        else:
            vmin = int(np.nanmin(vals[vals >= 0]))
            vmax = int(np.nanmax(vals[vals < np.inf]))
            formatter = BasicTickFormatter()
            color_map = 'Viridis256'
            ticker = BasicTicker()

        # Set the plot scale
        if scale == 'log':
            mapper = LogColorMapper(palette=color_map, low=vmin, high=vmax)
        else:
            mapper = LinearColorMapper(palette=color_map, low=vmin, high=vmax)

        # Plot the frame
        fig.image(source=source, image=ptype, x=0, y=0, dw=dw, dh=dh, color_mapper=mapper)
        color_bar = ColorBar(color_mapper=mapper, ticker=ticker, formatter=formatter, orientation="horizontal", location=(0, 0))

        # Plot the trace polynomials
        if trace_coeffs is not None:
            X = np.linspace(0, 2048, 2048)

            for coeffs in trace_coeffs:
                Y = np.polyval(coeffs, X)
                fig.line(X, Y, color='red')

        # Add the colorbar
        fig.add_layout(color_bar, 'below')

        # Add the figure to the tab list
        tabs.append(Panel(child=fig, title=pname))

    # Make the final tabbed figure
    final = Tabs(tabs=tabs)

    return final


def plot_frames(data, idx=0, scale='linear', trace_coeffs=None, saturation=0.8, width=1024, height=300, title=None, wavecal=None):
    """
    Plot a SOSS frame

    Parameters
    ----------
    data: sequence
        The 3D data to plot
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
    output_file('soss_frames.html')

    # Fix log scale plot values
    if scale == 'log':
        data[data < 1.] = 1.

    # Get data, snr, and saturation for plotting
    dims = data.shape
    vmin = int(np.nanmin(data[data >= 0]))
    vmax = int(np.nanmax(data[data < np.inf]))
    dat = data
    snr = np.sqrt(data)
    fullWell = 65536.0
    sat = dat > saturation * fullWell
    sat = sat.astype(int)

    # Wrap the data in two ColumnDataSources
    frames = range(len(dat))
    source_available = ColumnDataSource(data=dict(**{'counts{}'.format(n): dat[n] for n in frames}, **{'snr{}'.format(n): snr[n] for n in frames}, **{'saturation{}'.format(n): sat[n] for n in frames}))
    source_visible = ColumnDataSource(data=dict(counts=[dat[idx]], snr=[snr[idx]], saturation=[sat[idx]]))

    # Set the tooltips
    tooltips = [("(x,y)", "($x{int}, $y{int})"), ("ADU/s", "@counts"), ("SNR", "@snr"), ('Saturation', '@saturation')]

    # Add wavelength calibration if possible
    if isinstance(wavecal, np.ndarray):
        if wavecal.shape == dat[idx].shape:
            source_visible.data['wave1'] = [wavecal]
            tooltips.append(("Wavelength", "@wave1"))
        if wavecal.ndim == 3 and wavecal.shape[0] == 3:
            source_visible.data['wave1'] = [wavecal[0]]
            source_visible.data['wave2'] = [wavecal[1]]
            source_visible.data['wave3'] = [wavecal[2]]
            tooltips.append(("Wave 1", "@wave1"))
            tooltips.append(("Wave 2", "@wave2"))
            tooltips.append(("Wave 3", "@wave3"))

    # Make the figure
    fig = figure(x_range=(0, dims[2]), y_range=(0, dims[1]),
                 tooltips=tooltips, width=width, height=height,
                 title=title, toolbar_location='above', toolbar_sticky=True)

    # Plot the frame
    if scale == 'log':
        color_mapper = LogColorMapper(palette="Viridis256", low=vmin, high=vmax)
        fig.image(source=source_visible, image='counts', x=0, y=0, dw=dims[2], dh=dims[1], color_mapper=color_mapper)
        color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(), orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))

    else:
        color_mapper = LinearColorMapper(palette="Viridis256", low=vmin, high=vmax)
        fig.image(source=source_visible, image='counts', x=0, y=0, dw=dims[2], dh=dims[1], palette='Viridis256')
        color_bar = ColorBar(color_mapper=color_mapper, orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))

    # Plot the trace polynomials
    if trace_coeffs is not None:
        X = np.linspace(0, 2048, 2048)

        for coeffs in trace_coeffs:
            Y = np.polyval(coeffs, X)
            fig.line(X, Y, color='red')

    # Make the frame slider
    slider = Slider(title='Frame', value=idx, start=0, end=dims[0]-1, step=1)

    # CustomJS callback to update the three plots on slider changes
    callback = CustomJS(args=dict(visible=source_visible, available=source_available, slide=slider), code="""
        var vis = visible.data;
        var avail = available.data;
        var frame = slide.value;
        vis['counts'] = [avail['counts'.concat(frame.toString(10))]];
        vis['snr'] = [avail['snr'.concat(frame.toString(10))]];
        vis['saturation'] = [avail['saturation'.concat(frame.toString(10))]];
        visible.change.emit();
    """)

    # Add callback to spectrum slider
    slider.js_on_change('value', callback)

    return column(fig, slider)


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
    for c in col:
        color = next(utils.COLORS)
        fig.line(np.arange(slc[c, :].size), slc[c, :], color=color, legend='Column {}'.format(c))
        vline = Span(location=c, dimension='height', line_color=color, line_width=3)
        dfig.add_layout(vline)
    fig.xaxis.axis_label = 'Row'
    fig.yaxis.axis_label = 'Count Rate [ADU/s]'
    fig.legend.click_policy = 'mute'

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


def plot_time_series_spectra(fluxes, wavelength, width=1024, height=300, title=None, **kwargs):
    """
    Plot time series 1D spectra as an image

    Parameters
    ----------
    fluxes: array-like
        The 2D counts or flux
    wavelength: array-like
        The 1D wavelength array
    scale: str
        Plot scale, ['linear', 'log']
    width: int
        The width of the plot
    height: int
        The height of the plot
    title: str
        A title for the plot

    Returns
    -------
    bokeh.plotting.figure.Figure
        The figure
    """
    output_file('time_series_spectra.html')

    # Get plot params
    dh, dw = fluxes.shape
    fmin = np.nanmin(fluxes)
    fmax = np.nanmax(fluxes)
    wmin = np.nanmin(wavelength)
    wmax = np.nanmax(wavelength)
    lightcurves = fluxes.T

    # Set the source data
    sourceX = ColumnDataSource(data=dict(y=np.zeros(dw), wavelength=wavelength, flux=fluxes[0], **{'flux{}'.format(n): flux for n, flux in enumerate(fluxes)}))
    sourceY = ColumnDataSource(data=dict(x=np.zeros(dh), frames=np.arange(dh), lightcurve=lightcurves[0], **{'lightcurve{}'.format(round(w,3)): lc for w, lc in zip(wavelength, lightcurves)}))

    # ====================================================================

    # Make the 2D spectra figure
    spec_fig = figure(x_range=(wmin, wmax), y_range=(0, dh), x_axis_label='Wavelength', y_axis_label='Frame', width=width, height=height, title=title, toolbar_location='above', toolbar_sticky=True)

    # Plot the image
    fluxes[fluxes < 1.] = 1.
    color_mapper = LogColorMapper(palette="Viridis256", low=fmin, high=fmax)
    spec_fig.image(image=[fluxes], x=wmin, y=0, dw=wmax, dh=dh, color_mapper=color_mapper, alpha=0.8)
    color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(), orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))

    # Add current lightcurve line to plot
    spec_fig.line(x='x', y='frames', source=sourceY, color='red', line_width=3)

    # Add current spectrum line to plot
    spec_fig.line(x='wavelength', y='y', source=sourceX, color='blue', line_width=3)

    # ====================================================================

    # Make the 1D spectrum figure
    sp_fig = figure(x_range=(wmin, wmax), y_range=(fmin, fmax), width=width, height=height, x_axis_label='Wavelength', y_axis_label='Flux Density', title='Spectrum')

    # Draw the spectrum
    sp_fig.line('wavelength', 'flux', source=sourceX, color='blue', line_width=3, line_alpha=0.6)

    # Make the spectrum slider
    sp_slider = Slider(value=0, start=0, end=dh-1, step=1, width=30, title="Frame", orientation='vertical', direction='rtl', bar_color='blue')

    # ====================================================================

    # Make the 1D lightcurve figure
    lc_fig = figure(x_range=(0, dh), y_range=(0, dw), width=width, height=height, x_axis_label='Frame', y_axis_label='Flux Density', title='Lightcurve')

    # Draw the lightcurve
    lc_fig.line('frames', 'lightcurve', source=sourceY, color='red', line_width=3, line_alpha=0.6)

    # Make the lightcurve slider
    lc_slider = Slider(value=0, start=wmin, end=wmax, step=wavelength[1]-wavelength[0], width=width+40, title="Wavelength [um]", bar_color='red')

    # ====================================================================

    # CustomJS callback to update the three plots on slider changes
    callback = CustomJS(args=dict(sourcex=sourceX, sourcey=sourceY, sp_slide=sp_slider, lc_slide=lc_slider), code="""
        var datax = sourcex.data;
        var datay = sourcey.data;
        var sp = sp_slide.value;
        var lc = lc_slide.value;
        var wavelength = datax['wavelength'];
        var frames = datay['frames'];
        var x = datay['x'];
        var y = datax['y'];
        datax['flux'] = datax['flux'.concat(sp.toString(10))];
        datay['lightcurve'] = datay['lightcurve'.concat(lc.toFixed(3).toString(10))];
        for (var xi = 0; xi < x.length; xi++) {
            x[xi] = lc.toFixed(3);
        };
        for (var yi = 0; yi < y.length; yi++) {
            y[yi] = sp;
        };
        sourcex.change.emit();
        sourcey.change.emit();
    """)

    # Add callback to spectrum slider
    sp_slider.js_on_change('value', callback)

    # Add callback to lightcurve slider
    lc_slider.js_on_change('value', callback)

    return gridplot([[sp_fig, None],[spec_fig, sp_slider], [lc_slider, None], [lc_fig, None]])
