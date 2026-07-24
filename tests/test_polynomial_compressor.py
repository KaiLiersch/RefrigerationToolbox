"""Tests for the polaynomial-based ``PolynomialCompressor``."""


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


def test_verify_ploynomial_compressor(compressor):
    """Verification Testcase for Polynomial Compressor
    Polynomials and the reference solution are obtained using Bitzer Software for
    Bitzer ESH730Y scroll compressor at Te=0°C, Tc=30°C, sh=5K, sc=0K
    verification result is the compressor outlet temperature (48.2°C)"""
    compressor.calc(TE, TC, SH, SC)
    assert compressor.T_out == pytest.approx(48.2 + 273.15, rel=1e-2)


def test_calc_populates_states(compressor):
    for attr in ["T_in", "p_in", "h_in", "s_in", "d_in", "T_out", "p_out", "h_out", "s_out", "d_out"]:
        assert getattr(compressor, attr) is None
    compressor.calc(TE, TC, SH, SC)
    for attr in ["T_in", "p_in", "h_in", "s_in", "d_in", "T_out", "p_out", "h_out", "s_out", "d_out"]:
        assert getattr(compressor, attr) is not None


def test_subcooling_effects_the_discharge_enthalpy(compressor, fluid_ref, poly_coeffs):
    """The polynomial fixes the heating capacity, so subcooling shifts the whole high pressure side.

    Note, this test is non-physical because the polynomial for the cooling capacity is only ever
    valid for one subcooling value. The test is still valuable because it verifies that subcooling
    is taken into account. 
    """
    without = PolynomialCompressor(fluid_ref, **poly_coeffs)
    without.calc(TE, TC, SH, sc=0.0)
    h_out_without = without.h_out

    compressor.calc(TE, TC, SH, sc=5.0)

    assert compressor.h_out < h_out_without
    # Both share the same map, so the specific heating capacity is identical.
    assert compressor.Q_heat == pytest.approx(without.Q_heat)


