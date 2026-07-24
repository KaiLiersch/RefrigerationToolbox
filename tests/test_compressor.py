"""Tests for the abstract ``Compressor`` base class and its helpers."""

import numpy as np
import pytest

from refrigerationtoolbox.cycle.Compressor import Compressor, celsius_to_kelvin_coeffs


def _poly(coeffs, Te, Tc):
    """Evaluate the AHRI540 / EN12900 10-term polynomial."""
    vars_poly = np.array([1, Te, Tc, Te**2, Te * Tc, Tc**2, Te**3, Te**2 * Tc, Te * Tc**2, Tc**3])
    return float(np.dot(coeffs, vars_poly))


def test_celsius_to_kelvin_returns_ten_coeffs():
    coeffs_c = np.arange(1.0, 11.0)
    coeffs_k = celsius_to_kelvin_coeffs(coeffs_c)
    assert isinstance(coeffs_k, np.ndarray)
    assert coeffs_k.shape == (10,)


def test_celsius_to_kelvin_preserves_polynomial_value():
    """A polynomial in Celsius must equal the converted one evaluated in Kelvin."""
    coeffs_c = np.array([433.4, 14.99, -0.33, 0.20, 0.0078, -0.0055, 0.0015, 3e-5, -1.8e-4, -7.6e-7])
    coeffs_k = celsius_to_kelvin_coeffs(coeffs_c)

    for Te_c, Tc_c in [(0.0, 30.0), (-10.0, 45.0), (5.0, 55.0)]:
        val_c = _poly(coeffs_c, Te_c, Tc_c)
        val_k = _poly(coeffs_k, Te_c + 273.15, Tc_c + 273.15)
        assert val_k == pytest.approx(val_c, rel=1e-9)

