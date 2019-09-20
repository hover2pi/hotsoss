#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `utils` module."""

import unittest

import astropy.units as q
import numpy as np

from hotsoss import utils


class TestSpectralResponse(unittest.TestCase):
    """Test spectral_response function"""
    def setUp(self):
        """Test instance setup"""
        # Make wavelength for testing
        self.n = 100
        self.wave = np.linspace(1, 2, self.n)

    def test_response(self):
        """Test spectral_response function"""
        # CLEAR, SUBSTRIP256, order 1
        resp = utils.spectral_response(self.wave, filt='CLEAR', subarray='SUBSTRIP256', order=1, calfile=None)
        self.assertEqual(self.n, len(resp))

        # CLEAR, SUBSTRIP96, order 1
        resp = utils.spectral_response(self.wave, filt='CLEAR', subarray='SUBSTRIP96', order=1, calfile=None)
        self.assertEqual(self.n, len(resp))

        # CLEAR, FULL, order 1
        resp = utils.spectral_response(self.wave, filt='CLEAR', subarray='FULL', order=1, calfile=None)
        self.assertEqual(self.n, len(resp))

        # CLEAR, SUBSTRIP256, order 2
        resp = utils.spectral_response(self.wave, filt='CLEAR', subarray='SUBSTRIP256', order=2, calfile=None)
        self.assertEqual(self.n, len(resp))

        # CLEAR, SUBSTRIP96, order 2
        resp = utils.spectral_response(self.wave, filt='CLEAR', subarray='SUBSTRIP96', order=2, calfile=None)
        self.assertEqual(self.n, len(resp))

        # CLEAR, FULL, order 2
        resp = utils.spectral_response(self.wave, filt='CLEAR', subarray='FULL', order=2, calfile=None)
        self.assertEqual(self.n, len(resp))

        # F277W, SUBSTRIP256, order 1
        resp = utils.spectral_response(self.wave, filt='F277W', subarray='SUBSTRIP256', order=1, calfile=None)
        self.assertEqual(self.n, len(resp))

        # F277W, SUBSTRIP96, order 1
        resp = utils.spectral_response(self.wave, filt='F277W', subarray='SUBSTRIP96', order=1, calfile=None)
        self.assertEqual(self.n, len(resp))

        # F277W, FULL, order 1
        resp = utils.spectral_response(self.wave, filt='F277W', subarray='FULL', order=1, calfile=None)
        self.assertEqual(self.n, len(resp))

    def test_bad_input(self):
        """Check exceptions are raised"""
        # Test bad filt
        kwargs = {'wavelength': self.wave, 'filt': 'F430M'}
        self.assertRaises(ValueError, utils.spectral_response, **kwargs)

        # Test bad order
        kwargs = {'wavelength': self.wave, 'order': 4}
        self.assertRaises(ValueError, utils.spectral_response, **kwargs)


def test_colorgen():
    """Check color generator works"""
    colors = utils.color_gen()
    assert next(colors) != next(colors)


class TestCountsToFlux(unittest.TestCase):
    """Test counts_to_flux function"""
    def setUp(self):
        """Test instance setup"""
        # Make wavelength for testing
        self.n = 100
        self.wave = np.linspace(1, 2, self.n)
        self.counts = np.ones(self.n)*20000.

    def test_flux(self):
        """Test counts_to_flux function"""
        # CLEAR, SUBSTRIP256, order 1
        flux = utils.counts_to_flux(self.wave, self.counts, filt='CLEAR', subarray='SUBSTRIP256', order=1)
        self.assertEqual(self.n, len(flux))

        # CLEAR, SUBSTRIP96, order 1
        flux = utils.counts_to_flux(self.wave, self.counts, filt='CLEAR', subarray='SUBSTRIP96', order=1)
        self.assertEqual(self.n, len(flux))

        # CLEAR, FULL, order 1
        flux = utils.counts_to_flux(self.wave, self.counts, filt='CLEAR', subarray='FULL', order=1)
        self.assertEqual(self.n, len(flux))

        # CLEAR, SUBSTRIP256, order 2
        flux = utils.counts_to_flux(self.wave, self.counts, filt='CLEAR', subarray='SUBSTRIP256', order=2)
        self.assertEqual(self.n, len(flux))

        # CLEAR, SUBSTRIP96, order 2
        flux = utils.counts_to_flux(self.wave, self.counts, filt='CLEAR', subarray='SUBSTRIP96', order=2)
        self.assertEqual(self.n, len(flux))

        # CLEAR, FULL, order 2
        flux = utils.counts_to_flux(self.wave, self.counts, filt='CLEAR', subarray='FULL', order=2)
        self.assertEqual(self.n, len(flux))

        # F277W, SUBSTRIP256, order 1
        flux = utils.counts_to_flux(self.wave, self.counts, filt='F277W', subarray='SUBSTRIP256', order=1)
        self.assertEqual(self.n, len(flux))

        # F277W, SUBSTRIP96, order 1
        flux = utils.counts_to_flux(self.wave, self.counts, filt='F277W', subarray='SUBSTRIP96', order=1)
        self.assertEqual(self.n, len(flux))

        # F277W, FULL, order 1
        flux = utils.counts_to_flux(self.wave, self.counts, filt='F277W', subarray='FULL', order=1)
        self.assertEqual(self.n, len(flux))

    def test_bad_input(self):
        """Check exceptions are raised"""
        # Test bad filt
        kwargs = {'wavelength': self.wave, 'counts': self.counts, 'filt': 'F430M'}
        self.assertRaises(ValueError, utils.counts_to_flux, **kwargs)

        # Test bad order
        kwargs = {'wavelength': self.wave, 'counts': self.counts, 'order': 4}
        self.assertRaises(ValueError, utils.counts_to_flux, **kwargs)

        # Test bad subarray
        kwargs = {'wavelength': self.wave, 'counts': self.counts, 'subarray': None}
        self.assertRaises(ValueError, utils.counts_to_flux, **kwargs)

    def test_units(self):
        """Check the units are being processed"""
        # Fnu
        flux = utils.counts_to_flux(self.wave, self.counts, units=q.Jy)
        self.assertEqual(self.n, len(flux))

        # Bad units
        kwargs = {'wavelength': self.wave, 'counts': self.counts, 'units': q.cm**2}
        self.assertRaises(ValueError, utils.counts_to_flux, **kwargs)


def test_planet_data():
    """Check planet_data works"""
    data = utils.planet_data()
    assert len(data) == 2


def test_star_data():
    """Check star_data works"""
    data = utils.star_data()
    assert len(data) == 2


def test_transit_params():
    """Check transit_params model loads"""
    try:
        import batman

        time = np.linspace(0, 1, 100)
        tp = utils.transit_params(time)
        assert str(type(tp)) == "<class 'batman.transitmodel.TransitModel'>"

    except ImportError:
        pass


def test_subarray_specs():
    """Check subarray_specs works"""
    # SUBSTRIP256
    sub = utils.subarray_specs('SUBSTRIP256')
    assert sub['y'] == 256

    # SUBSTRIP96
    sub = utils.subarray_specs('SUBSTRIP96')
    assert sub['y'] == 96

    # FULL
    sub = utils.subarray_specs('FULL')
    assert sub['y'] == 2048


def test_wave_solutions():
    """Check wave_solutions works"""
    # SUBSTRIP256
    wav = utils.wave_solutions('SUBSTRIP256')
    assert wav.shape == (3, 256, 2048)

    # SUBSTRIP96
    wav = utils.wave_solutions('SUBSTRIP96')
    assert wav.shape == (3, 96, 2048)

    # FULL
    wav = utils.wave_solutions('FULL')
    assert wav.shape == (3, 2048, 2048)

    # Single order
    wav = utils.wave_solutions('FULL', order=1)
    assert wav.shape == (2048, 2048)


def test_globals():
    """Check the globals load"""
    colors = utils.COLORS
    star = utils.STAR_DATA
    planet = utils.PLANET_DATA
