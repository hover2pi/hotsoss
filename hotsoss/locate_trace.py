"""Module to locate the signal traces in SOSS 2D frames"""
import os
from functools import partial
from multiprocessing.dummy import Pool as ThreadPool
from pkg_resources import resource_filename
import time
import warnings

from astropy.io import fits
from bokeh.plotting import figure, show
from bokeh.layouts import column
from bokeh.models import Span
from bokeh.models.glyphs import Step
import numpy as np
from scipy.optimize import curve_fit

from . import crossdispersion as xdsp
from . import utils

warnings.simplefilter('ignore')


def isolate_signal(idx, frame, bounds=None, sigma=3, err=None, radius=None, filt='CLEAR', plot=False):
    """
    Fit a mixed gaussian function to the signal in a column of data. 
    Identify all pixels within n-sigma as signal.
    
    Parameters
    ----------
    idx: int
        The index of the column in the 2048 pixel wide subarray
    frame: array-like
         The 2D frame to pull the column from
    err: array-like (optional)
        The errors in the 1D data
    bounds: tuple
        A sequence of length-n (lower,upper) bounds on the n-parameters of func
    sigma: float
        The number of standard deviations to use when defining the signal
    err: int (optional)
        The uncertainty of the fit
    radius: int (optional)
        A constant radius for the trace at all wavelengths
    filt: str
        The filter used in the observations, ['CLEAR', 'F277W']
    plot: bool
        Plot the signal with the fit function
    
    Returns
    -------
    np.ndarray
        The values of signal pixels and the upper and lower bounds on the fit function
    """
    # Get the column of data
    col = frame[:, idx]

    # Make a mask for each order
    ord1 = np.ones_like(col)
    ord2 = np.ones_like(col)

    # Use the trace centers as the position of the center peak
    x1, x2 = trace_polynomial(evaluate=True)[:, idx]

    # Set the column at which order 2 ends
    order2end = 1900 if col.size == 256 else 1050

    # No order 2 if the trace is off the detector or the filter is F277W
    if idx >= order2end or filt.lower() == 'f277w':

        # Use the batman function to find only the first order
        func = xdsp.batman

        # Same as above, just one signal
        bounds = ([x1-3, 2, 100, 2, 300, 5], [x1+3, 4, 1e6, 8, 1e6, 10])

    # Otherwise there are two orders
    else:

        # Use the batmen function to find both orders
        func = xdsp.batmen

        # Set (lower, upper) bounds for each parameter
        # --------------------------------------------
        # x-position of the center peak, second psf
        # stanfard deviation of the center peak, second psf
        # amplitude of the center peak, second psf
        # stanfard deviation of the outer peaks, second psf
        # amplitude of the outer peaks, second psf
        # separation of the outer peaks from the center, second psf
        # x-position of the center peak, first psf
        # stanfard deviation of the center peak, first psf
        # amplitude of the center peak, first psf
        # stanfard deviation of the outer peaks, first psf
        # amplitude of the outer peaks, first psf
        # separation of the outer peaks from the center, first psf
        bounds = ([x2-3, 2, 3, 2, 3, 5, x1-3, 2, 100, 2, 300, 5],
                  [x2+3, 4, 1e3, 8, 1e3, 10, x1+3, 4, 1e6, 8, 1e6, 10])

    # Fit function to signal
    x = np.arange(col.size)
    params, cov = curve_fit(func, x, col, bounds=bounds, sigma=err)

    # Switch order 1 and 2 in list
    if len(params) == 12:
        params = np.concatenate([params[6:], params[:6]])

    # -----------------------------------------------------------------------
    # Order 1
    # -----------------------------------------------------------------------
    # Reduce to mixed gaussians with arguments (mu, sigma, A)
    p1 = params[:3]
    p2 = [params[0]-params[5], params[3], params[4]]
    p3 = [params[0]+params[5], params[3], params[4]]

    # Get the mu, sigma, and amplitude of each gaussian in each order
    params1 = np.array(sorted(np.array([p1,p2,p3]), key=lambda x: x[0]))

    # If radius is given use a set number of pixels as the radius
    if isinstance(radius, int):
        llim1 = params1[1][0]-radius
        ulim1 = params1[1][0]+radius

    # Otherwise, use the sigma value
    else:
        llim1 = params1[0][0]-params1[0][1]*sigma
        ulim1 = params1[-1][0]+params1[-1][1]*sigma

    # Unmask order 1
    ord1[(x > llim1) & (x < ulim1)] = 0

    # -----------------------------------------------------------------------
    # Order 2
    # -----------------------------------------------------------------------
    if func == xdsp.batmen:

        # Reduce to mixed gaussians with arguments (mu, sigma, A)
        p4 = params[6:9]
        p5 = [params[6]-params[11], params[9], params[10]]
        p6 = [params[6]+params[11], params[9], params[10]]

        # Get the mu, sigma, and amplitude of each gaussian in each order
        params2 = np.array(sorted(np.array([p4,p5,p6]), key=lambda x: x[0]))

        # If radius is given use a set number of pixels as the radius
        if isinstance(radius, int):
            llim2 = params2[1][0]-radius
            ulim2 = params2[1][0]+radius

        # Otherwise, use the sigma
        else:
            llim2 = params2[0][0]-params2[0][1]*sigma
            ulim2 = params2[-1][0]+params2[-1][1]*sigma

        # Unmask order 2
        ord2[(x > llim2) & (x < ulim2)] = 0

    # Mask the whole column for order 2
    else:
        p4 = p5 = p6 = llim2 = ulim2 = None

    # -----------------------------------------------------------------------

    if plot:

        # Make the figure
        fig = figure(x_range=(0, col.size), width=1000, height=600,
                     tooltips=[("x", "$x"), ("y", "$y")],
                     title='Column {}'.format(idx))
        fig.xaxis.axis_label = 'Row'
        fig.yaxis.axis_label = 'Count Rate [ADU/s]'

        # The data
        fig.step(x, col, color='black', legend='Data', mode='center')

        # Plot order 1 fit with limits
        color1 = '#2171b5'
        bm1 = xdsp.batman(x, *params[:6])
        fig.line(x, bm1, legend='Order 1', color=color1, line_width=2, alpha=0.8)
        low1 = Span(location=llim1, dimension='height', line_color=color1, line_dash='dashed')
        high1 = Span(location=ulim1, dimension='height', line_color=color1, line_dash='dashed')
        fig.renderers.extend([low1, high1])
        for g in [p1, p2, p3]:
            fig.line(x, xdsp.gaussian(x, *g), alpha=0.3, color=color1)

        # Plot order 2 fit with limits
        if llim2 is not None:
            color2 = '#DD4968'
            bm2 = xdsp.batman(x, *params[6:])
            fig.line(x, bm2, legend='Order 2', color=color2, line_width=2, alpha=0.8)
            low2 = Span(location=llim2, dimension='height', line_color=color2, line_dash='dashed')
            high2 = Span(location=ulim2, dimension='height', line_color=color2, line_dash='dashed')
            fig.renderers.extend([low2, high2])
            for g in [p4, p5, p6]:
                fig.line(x, xdsp.gaussian(x, *g), alpha=0.3, color=color2)

        fig.legend.click_policy = 'hide'
        show(fig)

    return ord1, ord2


def order_masks(frame, subarray='SUBSTRIP256', n_jobs=4, plot=False, save=False, **kwargs):
    """
    Generate a mask of the SOSS frame based on the column fits from isolate_signal

    Parameters
    ----------
    frame: array-like
        The 2D frame of SOSS data
    n_jobs: int
        The number of jobs for multiprocessing
    plot: bool
        Plot the masks
    save: bool
        Save the masks to file

    Returns
    -------
    sequence
        A masked frame for each trace order
    """
    # Get the file
    file = resource_filename('hotsoss', 'files/order_masks.npy')

    # Generate the trace masks
    if save:

        # Make 2 unmasked frames
        order1 = np.zeros_like(frame)
        order2 = np.zeros_like(frame)

        # Multiprocess spectral extraction for frames
        print("Coffee time! This takes about 3 minutes...")
        pool = ThreadPool(n_jobs)
        func = partial(isolate_signal, frame=frame, **kwargs)
        masks = pool.map(func, range(2048))
        pool.close()
        pool.join()

        # Set the mask in each column
        for n, (ord1, ord2) in enumerate(masks):
            order1[:, n] = ord1
            order2[:, n] = ord2

        # Save to file
        np.save(file, np.array([order1, order2]))

    # Or grab them from file
    else:

        # Open the file
        try:
            order1, order2 = np.load(file)
        except FileNotFoundError:
            print("No order mask file. Generating one now...")
            order1, order2 = order_masks(frame, save=True, plot=plot)
            return

    # Trim if SUBSTRIP96
    if subarray == 'SUBSTRIP96':
        order1 = order1[:96]
        order2 = order2[:96]

    if plot:

        # Make the figure
        height, width = order1.shape
        fig1 = figure(x_range=(0, width), y_range=(0, height),
                      tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                      width=int(width/2), height=height, title='Order 1 Mask')

        # Make the figure
        fig2 = figure(x_range=(0, width), y_range=(0, height),
                      tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")],
                      width=int(width/2), height=height, title='Order 2 Mask')

        # Plot the order mask
        fig1.image(image=[order1], x=0, y=0, dw=width, dh=height,
                   palette='Viridis256')

        # Plot the order mask
        fig2.image(image=[order2], x=0, y=0, dw=width, dh=height,
                   palette='Viridis256')

        show(column([fig1, fig2]))

    return order1, order2


def simulate_frame(filt='CLEAR', order_amps=(100, 10), plot=False):
    """
    Create a fake frame with traces

    Parameters
    ----------
    filt: str
        The filter to use, ['CLEAR', 'F277W']
    order_amps: sequence
        The amplitudes for each order
    """
    # Make empty frame
    nrows = 256
    frame = np.zeros((2048, nrows))

    # Get the traces
    traces = trace_polynomial(evaluate=True)
    waves = trace_wavelengths()

    # Determine cutoff for filter
    cutoff = 2.37 if filt == 'F277W' else 0

    # Get the batman PSF
    x = np.arange(40)
    psf = xdsp.batman(x, 20, 5, 10, 6, 50, 8)
    pix = int(len(psf)/2)

    # Iterate over the frame columns placing the psf, trimming where necessary
    for wave, amp, trace in zip(waves, order_amps, traces):
        for col, (w, y) in enumerate(zip(wave, trace)):
            if w > cutoff:
                y = int(y)
                start = max(0, y-pix)
                end = min(y+pix, nrows)
                frame[col, start:end] += psf[start-(y-pix):max(0, 2*pix-((y+pix)-end))]*amp

    # Transpose frame
    frame = frame.T

    if plot:

        # Make the figure
        fig = figure(x_range=(0, 2048), y_range=(0, nrows), width=1024, height=140)
        fig.image(image=[frame], x=0, y=0, dw=2048, dh=nrows, palette='Viridis256')
        show(fig)

    return frame


def trace_polynomial(subarray='SUBSTRIP256', order=None, evaluate=False, generate=False, poly_order=4):
    """The polynomial that describes the order trace

    Parameters
    ----------
    subarray:str
        The subarray to use, ['FULL', 'SUBSTRIP256', 'SUBSTRIP96']
    order: int, NoneType
        The order polynomial
    evaluate: bool
        Evaluate the polynomial on [0, 2048]
    poly_order: int
        The order polynomail to fit
    generate: bool
        Generate new coefficients

    Returns
    -------
    sequence
        The y values of the given order across the 2048 pixels
    """
    if generate:

        # Get the data
        file = resource_filename('hotsoss', 'files/soss_wavelength_trace_table1.txt')
        x1, y1, w1, x2, y2, w2 = np.genfromtxt(file, unpack=True)

        # Fit the polynomails
        fit1 = np.polyfit(x1, y1, poly_order)
        fit2 = np.polyfit(x2, y2, poly_order)

        coeffs = [fit1, fit2]

    else:

        # Coefficients to use
        coeffs = [[1.71164994e-11, -4.72119272e-08, 5.10276801e-05, -5.91535309e-02, 8.30680347e+01],
                  [2.35792131e-13, 2.42999478e-08, 1.03641247e-05, -3.63088657e-02, 9.96766537e+01]]

    # Specify the orders
    if order is None:
        slc = slice(None)
    elif isinstance(order, int) and order <= len(coeffs):
        slc = slice(order-1, order)
    else:
        raise ValueError("{}: order not understood".format(order))

    # Slice and dice
    result = coeffs[slc]

    # Evaluate
    if evaluate:
        return np.array([np.polyval(coeff, np.arange(2048)) for coeff in result]).squeeze()

    else:
        return result if order is None else result[0]


def trace_wavelengths(order=None, wavecal_file=None, npix=10, subarray='SUBSTRIP256'):
    """
    Return the average wavelength of N pixels in each column centered on trace
    """
    results = []

    # Get the pixel coordinates for the center of the trace
    trace_pixels = trace_polynomial(evaluate=True)

    # Load the wavelength calibration data
    wavecal = utils.wave_solutions(subarray=subarray, file=wavecal_file)[:2]

    # Get average wavelength of N pixels
    for pix, wave_map in zip(trace_pixels, wavecal):

        # Get average wavelength
        wave_cols = wave_map.swapaxes(0, 1)
        wavelengths = np.array([np.mean(col[i-npix:i+npix]) for i, col in zip(pix.astype(int), wave_cols)])

        # Add to the results
        results.append(wavelengths)

    # Select the right order
    if order in [1, 2]:
        results = results[order-1]

    return results


def wavelength_bins(save=False, subarray='SUBSTRIP256', wavecal_file=None):
    """Determine all the pixels in each wavelength bin for orders 1 and 2

    Parameters
    ----------
    save: bool
        Save the binned pixels to file
    subarray: str
        The subarray to use
    wavecal_file: str (optional)
        The path to the wavelength calibration file

    Returns
    -------
    list
        The (x, y) coordinates of all the pixels in each wavelength bin
    """
    file = resource_filename('hotsoss', 'files/wavelength_bins.npy')

    if save:

        # Get the pixel coordinates for the center of the trace
        trace_pixels = trace_polynomial(evaluate=True)

        # Load the wavelength calibration data
        wavecal = utils.wave_solutions(subarray=subarray, file=wavecal_file)[:2]
        wavelengths = trace_wavelengths(order=None, wavecal_file=wavecal_file, npix=10, subarray=subarray)
        signal_pixels = [[], []]

        # Store the pixel coordinates of each wavelength bin for each order
        for order, (pix, wave_map, wave_center) in enumerate(zip(trace_pixels, wavecal, wavelengths)):

            # Make a mask for each wavelength bin
            for n, w in enumerate(wave_center):

                # Edge cases
                try:
                    w0 = wave_center[n-1]
                except IndexError:
                    w0 = 0.1

                try:
                    w1 = wave_center[n+1]
                except IndexError:
                    w1 = 10

                # Define the width of the wavelength bin as half-way between neighboring points
                # and sort in case wavelength decreases to the right
                dw0 = np.mean([w0, w])
                dw1 = np.mean([w1, w])
                dw0, dw1 = sorted([dw0, dw1])

                # Isolate the signal pixels
                signal = np.where(np.logical_and(wave_map >= dw0, wave_map < dw1))

                # Add them to the list
                signal_pixels[order].append(signal)

        # Save to file
        np.save(file, signal_pixels)

    else:

        # Load from file
        try:
            signal_pixels = np.load(file)
        except FileNotFoundError:
            print("No wavelength bin file. Generating one now...")
            signal_pixels = wavelength_bins(save=True)

    return signal_pixels
