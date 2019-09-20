# -*- coding: utf-8 -*-
"""
A module to of utilities

Authors: Joe Filippazzo
"""

import itertools
from pkg_resources import resource_filename

import astropy.constants as ac
from astropy.io import fits
import astropy.units as q
from bokeh.palettes import Category20
import numpy as np


def spectral_response(wavelength, filt='CLEAR', subarray='SUBSTRIP256', order=1, calfile=None):
    """
    Get the relative spectral response for the given filter, order, and subarray

    Parameters
    ----------
    wavelength: array-like
        The wavelength array
    filt: str
        The filter name, ['CLEAR', 'F277W']
    subarray: str
        The subarray to use, ['FULL', 'SUBSTRIP256', 'SUBSTRIP96']
    order: int
        The dispersion order, [1, 2, 3]
    calfile: str (optional)
        A custom file to use for absolute calibration

    Returns
    -------
    np.ndarray
        The relative spectral response
    """
    # Valid filters
    filters = ['CLEAR', 'F277W']

    # Check the filter
    if not isinstance(filt, str) or filt not in filters:
        raise ValueError("'{}' not a supported filter. Try {}".format(filt, filters))

    # Check the order
    if not order in [1, 2, 3]:
        raise ValueError("{}: Not a valid order. Please select from [1, 2, 3].".format(order))

    # Get absolute calibration reference file
    calfile = calfile or resource_filename('hotsoss', 'files/niriss_ref_photom.fits')
    caldata = fits.getdata(calfile)

    # Get relative spectral response for the order (from
    # /grp/crds/jwst/references/jwst/jwst_niriss_photom_0028.fits)
    throughput = caldata[(caldata['pupil'] == 'GR700XD') & (caldata['filter'] == filt) & (caldata['order'] == order)]
    ph_wave = throughput.wavelength[throughput.wavelength > 0][1:-2]
    ph_resp = throughput.relresponse[throughput.wavelength > 0][1:-2]
    response = np.interp(wavelength, ph_wave, ph_resp)

    return response


def color_gen():
    """
    Generator for a Bokeh color palette
    """
    yield from itertools.cycle(Category20[20])


def counts_to_flux(wavelength, counts, filt='CLEAR', subarray='SUBSTRIP256', order=1, units=q.erg/q.s/q.cm**2/q.AA, **kwargs):
    """
    Convert the given counts to flux density in the given units

    Parameters
    ----------
    wavelength: array-like
        The wavelength values for each count
    counts: array-like
        The counts in [ADU/s]
    filt: str
        The filter name, ['CLEAR', 'F277W']
    subarray: str
        The subarray to use, ['FULL', 'SUBSTRIP256', 'SUBSTRIP96']
    order: int
        The dispersion order, [1, 2, 3]
    units: astropy.units.quantity.Quantity
        The flux density units for the result

    Returns
    -------
    np.ndarray
        The flux density
    """
    # Valid input
    filters = ['CLEAR', 'F277W']
    orders = [1, 2, 3]
    subarrays = ['SUBSTRIP256', 'SUBSTRIP96', 'FULL']

    # Check the filter
    if filt not in filters:
        raise ValueError("'{}' not a valid filter. Try {}".format(filt, filters))

    # Check the subarray
    if subarray not in subarrays:
        raise ValueError("'{}' not a valid subarray. Try {}".format(subarray, subarrays))

    # Check the order
    if order not in orders:
        raise ValueError("'{}' not a valid order. Try {}".format(order, orders))

    # Check the units
    if not (units.is_equivalent(q.erg/q.s/q.cm**2/q.AA) or units.is_equivalent(q.Jy)):
        raise ValueError('{}: Flux density must be in units of F_nu or F_lambda'.format(units))

    # Get the frame time for the given subarray
    frame_time = subarray_specs(subarray)['tfrm']

    # Get the relative spectral response
    response = spectral_response(wavelength, filt=filt, subarray=subarray, order=order, **kwargs)

    # Convert response in [mJy/ADU/s] to [Flam/ADU/s] or [Fnu/ADU/s]
    if units.is_equivalent(q.Jy):
        response = (response*q.mJy).to(units).value
    else:
        response = (response*q.mJy*ac.c/(wavelength*q.um)**2).to(units).value

    # Multiply response in [Flam/ADU/s] by counts in [ADU/s] for flux (is frame_time needed?)
    flux = counts*response/frame_time
    flux[flux == np.inf] = 0

    return flux


def planet_data():
    """
    Dummy data for time series simulations

    Returns
    -------
    sequence
        The wavelength and atmospheric transmission of the planet
    """
    planet_file = resource_filename('hotsoss', '/files/WASP107b_pandexo_input_spectrum.dat')
    planet = np.genfromtxt(planet_file, unpack=True)
    planet1D = [planet[0]*q.um, planet[1]]

    return planet1D


def star_data():
    """
    Dummy data for time series simulations

    Returns
    -------
    sequence
        The wavelength and flux of the star
    """
    star_file = resource_filename('hotsoss', 'files/scaled_spectrum.txt')
    star = np.genfromtxt(star_file, unpack=True)
    star1D = [star[0]*q.um, (star[1]*q.W/q.m**2/q.um).to(q.erg/q.s/q.cm**2/q.AA)]

    return star1D


def subarray_specs(subarr):
    """
    Get the pixel information for a NIRISS subarray

    The returned dictionary defines the extent ('x' and 'y'),
    the starting pixel ('xloc' and 'yloc'), and the number
    of reference pixels at each subarray edge ('x1', 'x2',
    'y1', 'y2) as defined by SSB/DMS coordinates shown below:
        ___________________________________
       |               y2                  |
       |                                   |
       |                                   |
       | x1                             x2 |
       |                                   |
       |               y1                  |
       |___________________________________|
    (1,1)

    Parameters
    ----------
    subarr: str
        The subarray name

    Returns
    -------
    dict
        The dictionary of the specified subarray
        or a nested dictionary of all subarrays
    """
    pix = {'FULL': {'xloc': 1, 'x': 2048, 'x1': 4, 'x2': 4, 'yloc': 1, 'y': 2048, 'y1': 4, 'y2': 4, 'tfrm': 10.737, 'tgrp': 10.737},
           'SUBSTRIP96': {'xloc': 1, 'x': 2048, 'x1': 4, 'x2': 4, 'yloc': 1803, 'y': 96, 'y1': 0, 'y2': 0, 'tfrm': 2.213, 'tgrp': 2.213},
           'SUBSTRIP256': {'xloc': 1, 'x': 2048, 'x1': 4, 'x2': 4, 'yloc': 1793, 'y': 256, 'y1': 0, 'y2': 4, 'tfrm': 5.491, 'tgrp': 5.491}}

    return pix[subarr]


def transit_params(time, teff=3500, logg=5., feh=0.):
    """
    Dummy transit parameters for time series simulations

    Parameters
    ----------
    time: sequence
        The time axis of the transit observation
    teff: float
        The effective temperature of the host star
    logg: float
        The surface gravity of the host star
    feh: float
        The metallicity of the host star

    Returns
    -------
    batman.transitmodel.TransitModel
        The transit model
    """
    try:
        import batman
        params = batman.TransitParams()
        params.t0 = 0.                                # time of inferior conjunction
        params.per = 5.7214742                        # orbital period (days)
        params.a = 0.0558*q.AU.to(q.R_sun)*0.66       # semi-major axis (in units of stellar radii)
        params.inc = 89.8                             # orbital inclination (in degrees)
        params.ecc = 0.                               # eccentricity
        params.w = 90.                                # longitude of periastron (in degrees)
        params.limb_dark = 'quadratic'                # limb darkening profile to use
        params.u = [0.1, 0.1]                         # limb darkening coefficients
        params.rp = 0.                                # planet radius (placeholder)
        tmodel = batman.TransitModel(params, time)
        tmodel.teff = teff                            # effective temperature of the host star
        tmodel.logg = logg                            # log surface gravity of the host star
        tmodel.feh = feh                              # metallicity of the host star

        return tmodel

    except ImportError:
        return None

def wave_solutions(subarray=None, order=None, file=None):
    """
    Get the wavelength maps for SOSS orders 1, 2, and 3
    This will be obsolete once the apply_wcs step of the JWST pipeline
    is in place.

    Parameters
    ----------
    subarray: str
        The subarray to return, ['SUBSTRIP96', 'SUBSTRIP256', 'FULL']
    order: int (optional)
        The trace order, [1, 2, 3]
    file: str
        The file containing the wavelength calibration

    Returns
    -------
    np.ndarray
        An array of the wavelength solutions for orders 1, 2, and 3
    """
    # Get the directory
    if file is None:
        file = resource_filename('hotsoss', '/files/soss_wavelengths_fullframe.fits')

    # Trim to the correct subarray
    if subarray == 'SUBSTRIP256':
        idx = slice(0, 256)
    elif subarray == 'SUBSTRIP96':
        idx = slice(160, 256)
    else:
        idx = slice(0, 2048)

    # Select the right order
    if order in [1, 2]:
        order = int(order)-1
    else:
        order = slice(0, 3)

    # Get the data from file and trim
    wave = fits.getdata(file).swapaxes(-2, -1)[order, idx, ::-1]

    return wave


COLORS = color_gen()
STAR_DATA = star_data()
PLANET_DATA = planet_data()
