"""Tests for the map-based ``PolynomialCompressor``."""

import numpy as np
import pytest

from refrigerationtoolbox.cycle.PolynomialCompressor import PolynomialCompressor

TE = 273.15
TC = 303.15
SH = 5.0
SC = 0.0


@pytest.fixture
def compressor(fluid_ref, poly_coeffs):
    return PolynomialCompressor(fluid_ref, **poly_coeffs)


def _poly(coeffs, Te, Tc):
    vars_poly = np.array([1, Te, Tc, Te**2, Te * Tc, Tc**2, Te**3, Te**2 * Tc, Te * Tc**2, Tc**3])
    return float(np.dot(coeffs, vars_poly))


def test_calc_matches_polynomial_evaluation(compressor, poly_coeffs):
    compressor.calc(TE, TC, SH, SC)
    assert compressor.m_flow == pytest.approx(_poly(poly_coeffs["m_flow_coeffs"], TE, TC))
    assert compressor.Q_ref == pytest.approx(_poly(poly_coeffs["Q_ref_coeffs"], TE, TC))
    assert compressor.P == pytest.approx(_poly(poly_coeffs["P_coeffs"], TE, TC))


def test_heat_is_sum_of_power_and_cooling(compressor):
    """Energy balance across the compressor+condenser map: Q_heat = P + Q_ref."""
    compressor.calc(TE, TC, SH, SC)
    assert compressor.Q_heat == pytest.approx(compressor.P + compressor.Q_ref)


def test_calc_populates_states(compressor):
    compressor.calc(TE, TC, SH, SC)
    for attr in ["T_in", "p_in", "h_in", "s_in", "d_in", "T_out", "p_out", "h_out", "s_out", "d_out"]:
        assert getattr(compressor, attr) is not None


def test_outlet_pressure_exceeds_inlet(compressor):
    compressor.calc(TE, TC, SH, SC)
    assert compressor.p_out > compressor.p_in


def test_str_includes_polynomial_duties(compressor):
    compressor.calc(TE, TC, SH, SC)
    text = str(compressor)
    assert "Q_heat" in text
    assert "Q_ref" in text
    assert "PolynomialCompressor" in text
