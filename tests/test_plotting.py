#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `plotting` module."""

import unittest

import numpy as np

from hotsoss import plotting as plt
from hotsoss import locate_trace as lt


class TestPlotFrame(unittest.TestCase):
    """Test plot_frame function"""
    def setUp(self):
        """Test instance setup"""
        # Make frame for testing
        self.frame = np.ones((256, 2048))

    def test_linearscale(self):
        """Test that a plot of linear data can be creted"""
        fig = plt.plot_frame(self.frame, scale='linear')
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Tabs'>")

    def test_logscale(self):
        """Test that a plot of log scale can be created"""
        fig = plt.plot_frame(self.frame, scale='log')
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Tabs'>")

    def test_coeffs(self):
        """Check the traces are drawn"""
        coeffs = lt.trace_polynomial()
        fig = plt.plot_frame(self.frame, scale='log', trace_coeffs=coeffs)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Tabs'>")

    def test_wavecal(self):
        """Check a wavecal array can be supplied"""
        # Order 1 only
        wavecal = np.ones((256, 2048))
        fig = plt.plot_frame(self.frame, wavecal=wavecal)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Tabs'>")

        # All orders
        wavecal = np.ones((3, 256, 2048))
        fig = plt.plot_frame(self.frame, wavecal=wavecal)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Tabs'>")


class TestPlotFrames(unittest.TestCase):
    """Test plot_frames function"""
    def setUp(self):
        """Test instance setup"""
        # Make frame for testing
        self.frames = np.ones((4, 256, 2048))

    def test_single_frame(self):
        """Test edge case of one frame"""
        frames = self.frames[:1, :, :]
        fig = plt.plot_frames(self.frames, scale='linear')
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")

    def test_linearscale(self):
        """Test that a plot of linear data can be creted"""
        fig = plt.plot_frames(self.frames, scale='linear')
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")

    def test_logscale(self):
        """Test that a plot of log scale can be created"""
        fig = plt.plot_frames(self.frames, scale='log')
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")

    def test_coeffs(self):
        """Check the traces are drawn"""
        coeffs = lt.trace_polynomial()
        fig = plt.plot_frames(self.frames, scale='log', trace_coeffs=coeffs)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")

    def test_wavecal(self):
        """Check a wavecal array can be supplied"""
        # Order 1 only
        wavecal = np.ones((256, 2048))
        fig = plt.plot_frames(self.frames, wavecal=wavecal)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")

        # All orders
        wavecal = np.ones((3, 256, 2048))
        fig = plt.plot_frames(self.frames, wavecal=wavecal)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")


class TestPlotRamp(unittest.TestCase):
    """Test plot_ramp function"""
    def setUp(self):
        """Test instance setup"""
        # Make frames for testing
        self.frames = np.ones((4, 256, 2048))

    def test_ramp(self):
        """Test that a ramp plot can be created"""
        fig = plt.plot_ramp(self.frames)
        self.assertEqual(str(type(fig)), "<class 'bokeh.plotting.figure.Figure'>")


class TestPlotSpectrum(unittest.TestCase):
    """Test plot_spectrum function"""
    def setUp(self):
        """Test instance setup"""
        # Make spectrum for testing
        self.wave = np.linspace(1, 2, 100)
        self.flux = np.linspace(1, 1.1, 100)

    def test_spectrum(self):
        """Test that a spectrum plot can be created"""
        # No figure
        fig = plt.plot_spectrum(self.wave, self.flux)
        self.assertEqual(str(type(fig)), "<class 'bokeh.plotting.figure.Figure'>")

        # With figure
        fig2 = plt.plot_spectrum(self.wave, self.flux, fig=fig)
        self.assertEqual(str(type(fig2)), "<class 'bokeh.plotting.figure.Figure'>")


class TestPlotTimeSeriesSpectra(unittest.TestCase):
    """Test plot_time_series_spectra function"""
    def setUp(self):
        """Test instance setup"""
        # Make spectrum for testing
        self.wave = np.linspace(1, 2, 2048)
        self.flux = np.random.normal(loc=1000., size=(100, 2048), scale=10.)
        self.time = np.arange(100)

    def test_tso(self):
        """Test that a spectrum plot can be created"""
        # Working plot
        fig = plt.plot_time_series_spectra(self.flux, wavelength=self.wave, time=self.time)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")

    def test_bad_input(self):
        """Test the function fails properly"""
        # Bad flux
        self.assertRaises(ValueError, plt.plot_time_series_spectra, np.ones((2, 3, 4, 5)))

        # Bad wavelength
        kwargs = {'flux': np.ones((2, 2048)), 'wavelength': np.arange(35)}
        self.assertRaises(ValueError, plt.plot_time_series_spectra, **kwargs)

        # Bad time
        kwargs = {'flux': np.ones((2, 2048)), 'time': np.arange(35)}
        self.assertRaises(ValueError, plt.plot_time_series_spectra, **kwargs)
