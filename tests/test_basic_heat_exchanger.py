"""Tests for ``BasicHeatExchanger`` (lumped energy-balance model)."""

import pytest

from refrigerationtoolbox.cycle.BasicHeatExchanger import BasicHeatExchanger
from tests.conftest import COND_STATE


@pytest.fixture
def condenser(fluid_ref, fluid_sec):
    return BasicHeatExchanger(fluid_ref, fluid_sec)


def test_state_none_before_calc(condenser):
    assert condenser.Q is None
    assert condenser.T_ref_in is None


def test_condenser_state_detected(condenser):
    condenser.calc(**COND_STATE)
    assert condenser.state == "cond"


def test_heat_duty_definition(condenser):
    """Q is defined directly from the refrigerant enthalpy drop and flow."""
    condenser.calc(**COND_STATE)
    expected = COND_STATE["m_flow_ref"] * (COND_STATE["h_ref_in"] - COND_STATE["h_ref_out"])
    assert condenser.Q == pytest.approx(expected)
    assert condenser.Q > 0  # condenser rejects heat


def test_secondary_energy_balance(condenser):
    """Heat rejected by the refrigerant is absorbed by the secondary fluid."""
    condenser.calc(**COND_STATE)
    q_secondary = condenser.m_flow_sec * (condenser.h_sec_out - condenser.h_sec_in)
    assert q_secondary == pytest.approx(condenser.Q, rel=1e-9)


def test_secondary_heats_up_in_condenser(condenser):
    condenser.calc(**COND_STATE)
    assert condenser.T_sec_out > condenser.T_sec_in


def test_evaporator_state_detected(fluid_ref, fluid_sec):
    """Swapping inlet/outlet enthalpies flips the exchanger into evaporator mode."""
    evap = BasicHeatExchanger(fluid_ref, fluid_sec)
    state = dict(COND_STATE)
    state["h_ref_in"], state["h_ref_out"] = state["h_ref_out"], state["h_ref_in"]
    state["T_sec_in"] = 10 + 273.15
    evap.calc(**state)
    assert evap.state == "evap"
    assert evap.Q < 0  # h_ref_in < h_ref_out
