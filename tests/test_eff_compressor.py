"""Tests for the efficiency-based ``EffCompressor``."""

import CoolProp.CoolProp as CP
import pytest
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.EffCompressor import EffCompressor

# Operating point: evaporating at 0 degC, condensing at 30 degC, 5 K superheat.
TE = 273.15
TC = 303.15
SH = 5.0
SC = 0.0


@pytest.fixture
def compressor(fluid_ref):
    return EffCompressor(fluid_ref, is_eff=0.8, vol_eff=0.9, N=1500 / 60, Vd=5e-5)


def test_verify_efficiency_compressor_no_sh(compressor):
    """Verify solution against TLK Energy log(p)-h diagram without superheat
    for R134a, TE=0°C, TC=30°C, SH=0 K, SC=0 K https://tlk-energy.de/en/phase-diagrams/pressure-enthalpy"""
    compressor.calc(TE, TC, 0, SC)
    assert compressor.T_out == pytest.approx(38.5 + 273.15, abs=1e-1)


def test_verify_efficiency_compressor_sh(compressor):
    """Verify solution against TLK Energy log(p)-h diagram with superheat
    for R134a, TE=0°C, TC=30°C, SH=5 K, SC=0 K https://tlk-energy.de/en/phase-diagrams/pressure-enthalpy"""
    compressor.calc(TE, TC, SH, SC)
    assert compressor.T_out == pytest.approx(43.5 + 273.15, abs=1e-1)


def test_state_none_before_calc(compressor):
    assert compressor.m_flow is None
    assert compressor.P is None
    assert compressor.h_out is None


def test_calc_populates_states(compressor):
    for attr in ["T_in", "p_in", "h_in", "s_in", "d_in", "T_out", "p_out", "h_out", "s_out", "d_out"]:
        assert getattr(compressor, attr) is None
    compressor.calc(TE, TC, SH, SC)
    for attr in ["T_in", "p_in", "h_in", "s_in", "d_in", "T_out", "p_out", "h_out", "s_out", "d_out"]:
        assert getattr(compressor, attr) is not None


def test_compression_raises_enthalpy_and_pressure(compressor):
    compressor.calc(TE, TC, SH, SC)
    # Compression adds work: outlet enthalpy and pressure exceed inlet.
    assert compressor.h_out > compressor.h_in
    assert compressor.p_out > compressor.p_in


def test_lower_efficiency_needs_more_power():
    """A less efficient compressor must consume more power for the same duty."""
    good = EffCompressor(AbstractState("HEOS", "R134a"), is_eff=0.9, vol_eff=0.9, N=1500 / 60, Vd=5e-5)
    good.calc(TE, TC, SH, SC)

    bad = EffCompressor(AbstractState("HEOS", "R134a"), is_eff=0.5, vol_eff=0.9, N=1500 / 60, Vd=5e-5)
    bad.calc(TE, TC, SH, SC)

    assert bad.P > good.P


def test_inlet_pressure_is_saturation_at_te(compressor, fluid_ref):
    compressor.calc(TE, TC, SH, SC)
    fluid_ref.update(CP.QT_INPUTS, 1.0, TE)
    assert compressor.p_in == pytest.approx(fluid_ref.p())
