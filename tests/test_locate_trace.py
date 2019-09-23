#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `locate_trace` module."""

import unittest

import numpy as np

from hotsoss import locate_trace as lt


def test_simulate_frame():
    """Test the simulate_frame function"""
    # CLEAR and plot test
    assert lt.simulate_frame(plot=True).shape == (256, 2048)

    # F277W test
    assert lt.simulate_frame(filt='F277W').shape == (256, 2048)


def test_isolate_signal():
    """Test isolate_signal function"""
    # Make frame for testing
    frame = lt.simulate_frame()

    # CLEAR and plot
    assert len(lt.isolate_signal(500, frame, plot=True)) == 2

    # F277W and radius
    assert len(lt.isolate_signal(500, frame, filt='F277W', radius=20)) == 2


def test_order_masks():
    """Test order_masks function"""
    # Make frame for testing
    frame = lt.simulate_frame()

    # Test plot and defaults
    assert len(lt.order_masks(frame, plot=True)) == 2

    # Test save and subarray
    assert len(lt.order_masks(frame, subarray='SUBSTRIP96', save=True)) == 2


def test_trace_polynomial():
    """Test trace_polynomial function"""
    # No order specified
    assert len(lt.trace_polynomial(order=None, evaluate=False)) == 2

    # Single order
    assert len(lt.trace_polynomial(order=1, evaluate=False)) == 5

    # Test evaluate
    assert len(lt.trace_polynomial(order=1, evaluate=True)) == 2048


def test_trace_wavelengths():
    """Test trace_wavelengths function"""
    # No order specified
    assert len(lt.trace_wavelengths(order=None)) == 2

    # Single order
    assert len(lt.trace_wavelengths(order=1)) == 2048


def test_wavelength_bins():
    """Test wavelength_bins works"""
    # Default values for two orders
    assert len(lt.wavelength_bins()) == 2

    # Generate
    assert len(lt.wavelength_bins(save=True)) == 2
