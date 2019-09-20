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
        self.assertEqual(str(type(fig)), "<class 'bokeh.plotting.figure.Figure'>")

    def test_logscale(self):
        """Test that a plot of log scale can be created"""
        fig = plt.plot_frame(self.frame, scale='log')
        self.assertEqual(str(type(fig)), "<class 'bokeh.plotting.figure.Figure'>")

    def test_coeffs(self):
        """Check the traces are drawn"""
        coeffs = lt.trace_polynomial()
        fig = plt.plot_frame(self.frame, scale='log', trace_coeffs=coeffs)
        self.assertEqual(str(type(fig)), "<class 'bokeh.plotting.figure.Figure'>")

    def test_wavecal(self):
        """Check a wavecal array can be supplied"""
        # Order 1 only
        wavecal = np.ones((256, 2048))
        fig = plt.plot_frame(self.frame, wavecal=wavecal)
        self.assertEqual(str(type(fig)), "<class 'bokeh.plotting.figure.Figure'>")

        # All orders
        wavecal = np.ones((3, 256, 2048))
        fig = plt.plot_frame(self.frame, wavecal=wavecal)
        self.assertEqual(str(type(fig)), "<class 'bokeh.plotting.figure.Figure'>")


class TestPlotSlice(unittest.TestCase):
    """Test plot_slice function"""
    def setUp(self):
        """Test instance setup"""
        # Make frame for testing
        self.frame = np.ones((256, 2048))

    def test_col(self):
        """Test that a plot of a given column can be created"""
        # Single column
        fig = plt.plot_slice(self.frame, 500)
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")

        # Multiple columns
        fig = plt.plot_slice(self.frame, [500, 1000, 1500])
        self.assertEqual(str(type(fig)), "<class 'bokeh.models.layouts.Column'>")


class TestPlotRamp(unittest.TestCase):
    """Test plot_ramp function"""
    def setUp(self):
        """Test instance setup"""
        # Make frames for testing
        self.frames = np.ones((2, 2 256, 2048))

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