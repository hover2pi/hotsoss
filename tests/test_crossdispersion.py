#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `crossdispersion` module."""

import unittest

import numpy as np

from hotsoss import crossdispersion as xd


class TestCrossDiespersion(unittest.TestCase):
    """Test crossdispersion.py functions"""
    def setUp(self):
        """Test instance setup"""
        # Make x axis for testing
        self.n = 100
        self.x = np.linspace(0, 100, self.n)
        self.mu0 = 50.
        self.sigma0 = 1.
        self.A0 = 10.
        self.sigma1 = 0.5
        self.A1 = 50.
        self.sep = 5.

    def test_batman(self):
        """Check batman function works"""
        result = xd.batman(self.x, self.mu0, self.sigma0, self.A0, self.sigma1, self.A1, self.sep)
        self.assertEqual(self.n, len(result))

    def test_batmen(self):
        """Check batmen function works"""
        result = xd.batmen(self.x, self.mu0, self.sigma0, self.A0, self.sigma1, self.A1, self.sep, self.mu0+10., self.sigma0, self.A0/2., self.sigma1, self.A1/2., self.sep)
        self.assertEqual(self.n, len(result))

    def test_bimodal(self):
        """Check bimodal function works"""
        result = xd.bimodal(self.x, self.mu0, self.sigma0, self.A0, self.mu0+10., self.sigma1, self.A1)
        self.assertEqual(self.n, len(result))

    def test_gaussian(self):
        """Check gaussian function works"""
        result = xd.gaussian(self.x, self.mu0, self.sigma0, self.A0)
        self.assertEqual(self.n, len(result))

    def test_trimodal(self):
        """Check trimodal function works"""
        result = xd.trimodal(self.x, self.mu0, self.sigma0, self.A0, self.mu0+10., self.sigma0, self.A0, self.mu0-10., self.sigma0, self.A0)
        self.assertEqual(self.n, len(result))
