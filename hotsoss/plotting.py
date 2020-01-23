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
from bokeh.models import ColumnDataSource, HoverTool, LogColorMapper, FixedTicker, BasicTickFormatter, FuncTickFormatter, BasicTicker, LogTicker, LinearColorMapper, ColorBar, Span, CustomJS, Slider, Range1d
from bokeh.models.widgets import Panel, Tabs
from bokeh.layouts import gridplot, column
import numpy as np

from . import utils
from . import locate_trace as lt


def plot_frame(frame, cols=0, scale='linear', trace_coeffs=None, saturation=0.8, title=None, wavecal=None):
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


def plot_frames(data, idx=0, col=0, scale='linear', trace_coeffs=None, saturation=0.8, width=1000, title=None, wavecal=None):
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
    # Determine subarray
    nframes, nrows, ncols = data.shape

    # Remove infs
    data[data == np.inf] = np.nan

    # Get data, snr, and saturation for plotting
    dat = data
    snr = np.sqrt(data)
    fullWell = 65536.0
    sat = dat > saturation * fullWell
    sat = sat.astype(int)

    # Fix log scale plot values
    if scale == 'log':
        dat[dat < 1.] = 1.
        snr[snr < 1.] = 1.

    # Broadcast the data
    frames = np.arange(nframes)
    columns = np.arange(ncols)
    rows = np.arange(nrows)
    verticals = np.tile(np.arange(ncols), (nrows, 1)).T

    # Wrap the data in ColumnDataSources
    source_available = ColumnDataSource(data=dict(**{'counts{}'.format(n): dat[n] for n in frames}, **{'snr{}'.format(n): snr[n] for n in frames}, **{'saturation{}'.format(n): sat[n] for n in frames}))
    source_visible = ColumnDataSource(data=dict(counts=[dat[idx]], snr=[snr[idx]], saturation=[sat[idx]]))
    vertical_available = ColumnDataSource(data=dict(**{'vertical{}'.format(n): vert for n, vert in enumerate(verticals)}))
    vertical_visible = ColumnDataSource(data=dict(column=rows, vertical=verticals[col]))
    col_visible = ColumnDataSource(data=dict(columns=rows, counts=dat[0, :, col], snr=snr[0, :, col], saturation=sat[0, :, col]))
    col_dict = {}
    for fnum in frames:
        for cnum in columns:
            for datacube, pname in zip([dat, snr, sat], ['counts', 'snr', 'saturation']):
                col_dict['{}{}_{}'.format(pname, cnum, fnum)] = datacube[fnum, :, cnum]
    col_available = ColumnDataSource(data=col_dict)

    # Set the tooltips
    tooltips = [("(x,y)", "($x{int}, $y{int})"), ("ADU/s", "@counts"), ("SNR", "@snr"), ('Saturation', '@saturation')]
    col_color = 'blue'

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

    # Set shared plot params
    x_range = (0, ncols)
    y_range = (0, nrows)
    height = int(nrows/2.)+160
    toolbar = 'right'

    # Draw the figures
    tabs = []
    for pdata, pname, ptype, ylabel in zip([dat, snr, sat], ['Counts', 'SNR', 'Saturation ({}% Full Well)'.format(saturation*100)], ['counts', 'snr', 'saturation'], ['Count Rate [ADU/s]', 'SNR', 'Saturated']):

        # Make the figure
        fig_title = '{} - {}'.format(title, pname)

        # Get min and max
        vmin = np.nanmin(pdata)
        vmax = np.nanmax(pdata)

        # Saturation plot is different
        if ptype == 'saturation':
            formatter = FuncTickFormatter(code="""return {0: 'Unsaturated', 1: 'Saturated'}[tick]""")
            color_map = ['#404387', '#FDE724']
            ticker = FixedTicker(ticks=[vmin, vmax])

        # Counts and SNR are similar plots
        else:
            formatter = BasicTickFormatter()
            color_map = 'Viridis256'
            ticker = BasicTicker()

        # Set the plot scale
        if scale == 'log':
            mapper = LogColorMapper(palette=color_map, low=vmin, high=vmax)
        else:
            mapper = LinearColorMapper(palette=color_map, low=vmin, high=vmax)

        # Plot the image data
        fig = figure(x_range=x_range, y_range=y_range, tooltips=tooltips, width=width, height=height, title=fig_title, toolbar_location=toolbar, toolbar_sticky=True)
        fig.image(source=source_visible, image=ptype, x=0, y=0, dw=ncols, dh=nrows, color_mapper=mapper)

        # Plot the line indicating the column
        fig.line('vertical', 'column', source=vertical_visible, color=col_color, line_width=3)

        # Add the colorbar
        color_bar = ColorBar(color_mapper=mapper, ticker=ticker, formatter=formatter, orientation="horizontal", location=(0, 0))
        fig.add_layout(color_bar, 'below')

        # Plot the column data
        col_fig = figure(width=width, height=300, toolbar_location=toolbar, toolbar_sticky=True)
        col_fig.step('columns', ptype, source=col_visible, color=col_color)
        col_fig.xaxis.axis_label = 'Row'
        col_fig.yaxis.axis_label = ylabel
        col_fig.y_range = Range1d(vmin*0.9, vmax*1.1)
        col_fig.x_range = Range1d(*y_range)

        # Plot the trace polynomials
        if trace_coeffs is not None:
            for coeffs in trace_coeffs:
                Y = np.polyval(coeffs, columns)
                fig.line(columns, Y, color='red')

        # Add the figure to the tab list
        tabs.append(Panel(child=column(fig, col_fig), title=pname))

    # Make the final tabbed figure
    final = Tabs(tabs=tabs)

    # Write JS code
    code ="""
        var vis = visible.data;
        var avail = available.data;
        var frame = fr_slide.value.toString(10);

        var viscol = col_vis.data;
        var availcol = col_avail.data;
        var col = col_slide.value.toString(10);

        var visvert = vert_vis.data;
        var availvert = vert_avail.data;

        vis['counts'] = [avail['counts'.concat(frame)]];
        vis['snr'] = [avail['snr'.concat(frame)]];
        vis['saturation'] = [avail['saturation'.concat(frame)]];

        viscol['counts'] = availcol['counts'.concat(col, '_', frame)];
        viscol['snr'] = availcol['snr'.concat(col, '_', frame)];
        viscol['saturation'] = availcol['saturation'.concat(col, '_', frame)];

        visvert['vertical'] = availvert['vertical'.concat(col)];

        visible.change.emit();
        col_vis.change.emit();
        vert_vis.change.emit();
    """

    # Make the column slider
    column_slider = Slider(title='Column', value=col, start=0, end=ncols-1, step=1)

    # Make the frame slider
    if nframes-1 > 0:
        frame_slider = Slider(title='Frame', value=idx, start=0, end=nframes-1, step=1)
    else:
        frame_slider = None
        code = code.replace('fr_slide.value.toString(10);', '0')

    # CustomJS callback to update the three plots on slider changes
    callback = CustomJS(args=dict(visible=source_visible, available=source_available, col_vis=col_visible, col_avail=col_available, vert_vis=vertical_visible, vert_avail=vertical_available, fr_slide=frame_slider, col_slide=column_slider), code=code)

    # Add callback to column slider
    column_slider.js_on_change('value', callback)

    # Add callback to frame slider
    if frame_slider is not None:
        frame_slider.js_on_change('value', callback)
        return column(final, frame_slider, column_slider)

    else:
        return column(final, column_slider)


def plot_ramp(data):
    """
    Plot the total flux on each frame to display the ramp
    """
    fig = figure()
    x = range(len(data))
    y = np.sum(data, axis=(-1, -2))
    fig.circle(x, y, size=8)
    fig.xaxis.axis_label = 'Group'
    fig.yaxis.axis_label = 'Count Rate [ADU/s]'

    return fig


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


def plot_time_series_spectra(flux, wavelength=None, time=None, xlabel='Column', ylabel='Frame', width=1024, height=300, title=None, **kwargs):
    """
    Plot time series 1D spectra as an image

    Parameters
    ----------
    flux: array-like
        The 2D counts or flux
    wavelength: sequence (optional)
        The 1D wavelength array
    time: sequence (optional)
        The 1D time array
    xlabel: str
        The label for the data x-axis
    ylabel: str
        The label for the data y-axis
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
    # Check that flux is 2D
    if not flux.ndim == 2:
        raise ValueError("{}: 'flux' must be a 2D array.".format(flux.shape))

    # Copy flux array
    flx = copy.copy(flux)

    # Get plot params
    dh, dw = flx.shape
    fmin = np.nanmin(flx)
    fmax = np.nanmax(flx)
    wmin, wmax = 0, dw
    tmin, tmax = 0, dh
    lightcurves = flx.T

    # Make sure time array is correct length
    if time is not None:
        if len(time) != dh:
            raise ValueError("{} != {}: 'time' array must be the same length as 'flux' cube.".format(len(time), flux.shape[0]))

    # Make sure wavelength array is correct length
    if wavelength is not None:
        if len(wavelength) != dw:
            raise ValueError("{} != {}: 'wavelength' array must be the same depth as 'flux' cube.".format(len(wavelength), flux.shape[1]))

    # Major tick labels
    waxis = np.arange(wmin, wmax)
    wstart, wskip = 128, 256
    taxis = np.arange(tmin, tmax)
    tstart, tskip = 1, max(4, min(tmax//2, 10))

    # Set the source data
    sourceX = ColumnDataSource(data=dict(wavelength=waxis, flux=flx[0], **{'flux{}'.format(n): fx for n, fx in enumerate(flx)}))
    sourceY = ColumnDataSource(data=dict(time=taxis, lightcurve=lightcurves[0], **{'lightcurve{}'.format(n): lc for n, lc in enumerate(lightcurves)}))
    sourceZ = ColumnDataSource(data=dict(x=[0.5], y=[0.5]))

    # Set the tools
    tools =  "pan,hover,wheel_zoom,box_zoom,reset"

    # ====================================================================

    # Make the 2D spectra figure
    spec_fig = figure(x_range=(wmin, wmax), y_range=(tmin, tmax), x_axis_label=xlabel, y_axis_label=ylabel, plot_width=width, plot_height=height, title=title, tools=tools, toolbar_location='above', toolbar_sticky=True, match_aspect=True)

    # Plot the image
    flx[flx < 1.] = 1.
    color_mapper = LogColorMapper(palette="Viridis256", low=fmin, high=fmax)
    spec_fig.image(image=[flx], x=wmin, y=tmin, dw=wmax, dh=tmax, color_mapper=color_mapper, alpha=0.8)
    color_bar = ColorBar(color_mapper=color_mapper, ticker=LogTicker(), orientation="horizontal", label_standoff=12, border_line_color=None, location=(0, 0))

    # Add current lightcurve line to plot
    spec_fig.vbar(x='x', width=1, top=tmax, source=sourceZ, color='red', alpha=0.3)

    # Add current spectrum line to plot
    spec_fig.hbar(y='y', height=1, right=wmax, source=sourceZ, color='blue', alpha=0.3)

    # Change y tick labels
    if time is not None:
        spec_fig.yaxis.ticker = taxis[tstart::tskip]
        spec_fig.yaxis.major_label_overrides = {int(n): '{:.2f}'.format(t) for n, t in zip(taxis[tstart::tskip], time[tstart::tskip])}

    # Change x tick labels
    if wavelength is not None:
        spec_fig.xaxis.ticker = waxis[wstart::wskip]
        spec_fig.xaxis.major_label_overrides = {int(n): '{:.3f}'.format(w) for n, w in zip(waxis[wstart::wskip], wavelength[wstart::wskip])}

    # ====================================================================

    # Make the 1D spectrum figure
    sp_fig = figure(x_range=(wmin, wmax), y_range=(fmin, fmax), width=width, height=height, x_axis_label=xlabel, y_axis_label='Flux Density', title='Spectrum')

    # Draw the spectrum
    sp_fig.step('wavelength', 'flux', source=sourceX, color='blue', line_width=3, line_alpha=0.6, mode='center')

    # Change x tick labels
    if wavelength is not None:
        sp_fig.xaxis.ticker = waxis[wstart::wskip]
        sp_fig.xaxis.major_label_overrides = {int(n): '{:.3f}'.format(w) for n, w in zip(waxis[wstart::wskip], wavelength[wstart::wskip])}

    # Make the spectrum slider
    sp_slider = Slider(value=0, start=tmin, end=tmax-1, step=1, width=30, title=ylabel, orientation='vertical', direction='rtl', bar_color='blue')

    # ====================================================================

    # Make the 1D lightcurve figure
    lc_fig = figure(x_range=(tmin, tmax), y_range=(fmin, fmax), width=width, height=height, x_axis_label=ylabel, y_axis_label='Flux Density', title='Lightcurve')

    # Draw the lightcurve
    lc_fig.step('time', 'lightcurve', source=sourceY, color='red', line_width=3, line_alpha=0.6, mode='center')

    # Change x tick labels
    if time is not None:
        lc_fig.xaxis.ticker = taxis[tstart::tskip]
        lc_fig.xaxis.major_label_overrides = {int(n): '{:.2f}'.format(t) for n, t in zip(taxis[tstart::tskip], time[tstart::tskip])}

    # Make the lightcurve slider
    lc_slider = Slider(value=0, start=wmin, end=wmax-1, step=1, width=width, title=xlabel, bar_color='red')

    # ====================================================================

    # CustomJS callback to update the three plots on slider changes
    callback = CustomJS(args=dict(sourcex=sourceX, sourcey=sourceY, sourcez=sourceZ, sp_slide=sp_slider, lc_slide=lc_slider), code="""
        var datax = sourcex.data;
        var datay = sourcey.data;
        var dataz = sourcez.data;
        var sp = sp_slide.value;
        var lc = lc_slide.value;
        var wavelength = datax['wavelength'];
        var time = datay['time'];
        datax['flux'] = datax['flux'.concat(sp.toString(10))];
        datay['lightcurve'] = datay['lightcurve'.concat(lc.toString(10))];
        dataz['x'] = [lc+0.5];
        dataz['y'] = [sp+0.5];
        sourcex.change.emit();
        sourcey.change.emit();
        sourcez.change.emit();
    """)

    # Add callback to spectrum slider
    sp_slider.js_on_change('value', callback)

    # Add callback to lightcurve slider
    lc_slider.js_on_change('value', callback)

    return gridplot([[sp_fig, None],[spec_fig, sp_slider], [lc_slider, None], [lc_fig, None]])
