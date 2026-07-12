"""Tests for the abstract ``HeatExchanger`` base class."""

import pytest

from refrigerationtoolbox.cycle.HeatExchanger import HeatExchanger


def test_heat_exchanger_is_abstract(fluid_ref, fluid_sec):
    """``calc`` is abstract, so the base class cannot be instantiated."""
    with pytest.raises(TypeError):
        HeatExchanger(fluid_ref, fluid_sec)


def test_base_calc_sets_inputs_and_state(fluid_ref, fluid_sec):
    """A minimal subclass can delegate to the base ``calc`` to store inputs."""

    class _Passthrough(HeatExchanger):
        def calc(self, m_flow_ref, m_flow_sec, h_ref_in, h_ref_out, p_ref, p_sec, T_sec_in):
            super().calc(m_flow_ref, m_flow_sec, h_ref_in, h_ref_out, p_ref, p_sec, T_sec_in)

    hx = _Passthrough(fluid_ref, fluid_sec)
    assert hx.state == "Undefined"

    hx.calc(0.01, 0.1, h_ref_in=4e5, h_ref_out=2e5, p_ref=7e5, p_sec=1e5, T_sec_in=293.15)
    assert hx.state == "cond"  # h_ref_in > h_ref_out
    assert hx.m_flow_ref == 0.01
    assert hx.p_sec == 1e5

    hx.calc(0.01, 0.1, h_ref_in=2e5, h_ref_out=4e5, p_ref=7e5, p_sec=1e5, T_sec_in=293.15)
    assert hx.state == "evap"  # h_ref_in < h_ref_out


def test_str_before_calc_is_na(fluid_ref, fluid_sec):
    class _Passthrough(HeatExchanger):
        def calc(self, *args, **kwargs):  # pragma: no cover - not exercised
            pass

    text = str(_Passthrough(fluid_ref, fluid_sec))
    assert "R134a" in text
    assert "Water" in text
    assert "N/A" in text
    assert "Undefined" in text
