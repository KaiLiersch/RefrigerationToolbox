"""Tests for the top-level ``Cycle`` orchestrator."""

import pytest
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.BasicHeatExchanger import BasicHeatExchanger
from refrigerationtoolbox.cycle.Cycle import Cycle
from refrigerationtoolbox.cycle.EffCompressor import EffCompressor
from refrigerationtoolbox.cycle.PolynomialCompressor import PolynomialCompressor

TE = 273.15
TC = 303.15


@pytest.fixture
def eff_cycle(fluid_ref):
    cycle = Cycle(fluid_ref, Te=TE, Tc=TC, sh=5.0, sc=0.0)
    cycle.set_compressor(EffCompressor(fluid_ref, is_eff=0.8, vol_eff=0.9, N=1500 / 60, Vd=5e-5))
    cycle.set_condenser(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.set_evaporator(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    return cycle


def test_states_none_before_calc(eff_cycle):
    assert eff_cycle.h1 is None
    assert eff_cycle.COP_heat is None


def test_setters_assign_components(fluid_ref):
    cycle = Cycle(fluid_ref, Te=TE, Tc=TC)
    comp = EffCompressor(fluid_ref, is_eff=0.8, vol_eff=0.9, N=25, Vd=5e-5)
    cycle.set_compressor(comp)
    assert cycle.compressor is comp


def test_calc_populates_all_states(eff_cycle):
    eff_cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    for i in (1, 2, 3, 4):
        assert getattr(eff_cycle, f"T{i}") is not None
        assert getattr(eff_cycle, f"h{i}") is not None


def test_pressure_levels_ordered(eff_cycle):
    eff_cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    # Condensing pressure is above evaporating pressure.
    assert eff_cycle.pc > eff_cycle.pe
    assert eff_cycle.p2 == pytest.approx(eff_cycle.pc)
    assert eff_cycle.p1 == pytest.approx(eff_cycle.pe)


def test_cop_definitions(eff_cycle):
    eff_cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    assert eff_cycle.COP_heat == pytest.approx(eff_cycle.Q_heat / eff_cycle.P_comp)
    assert eff_cycle.COP_cool == pytest.approx(eff_cycle.Q_ref / eff_cycle.P_comp)


def test_heating_cop_exceeds_cooling_cop_by_one(eff_cycle):
    """For an isenthalpic valve, Q_heat - Q_ref = P_comp, so COP_heat = COP_cool + 1."""
    eff_cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    assert eff_cycle.COP_heat == pytest.approx(eff_cycle.COP_cool + 1.0, rel=1e-6)


def test_cop_values_physically_reasonable(eff_cycle):
    eff_cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    assert eff_cycle.Q_heat > 0
    assert eff_cycle.Q_ref > 0
    assert eff_cycle.P_comp > 0
    assert eff_cycle.COP_cool > 1  # a working refrigeration cycle


def test_isenthalpic_expansion(eff_cycle):
    """With an efficiency compressor the valve is isenthalpic: h4 == h3."""
    eff_cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    assert eff_cycle.h4 == pytest.approx(eff_cycle.h3)


def test_polynomial_cycle_runs(fluid_ref, poly_coeffs):
    cycle = Cycle(fluid_ref, Te=TE, Tc=TC, sh=5.0, sc=0.0)
    cycle.set_compressor(PolynomialCompressor(fluid_ref, **poly_coeffs))
    cycle.set_condenser(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.set_evaporator(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.calc(m_flow_sec_cond=1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    # Polynomial duties are taken straight from the compressor map.
    assert cycle.Q_heat == pytest.approx(cycle.compressor.Q_heat)
    assert cycle.Q_ref == pytest.approx(cycle.compressor.Q_ref)
    assert cycle.COP_heat > 0


def test_calc_with_heat_exchangers(fluid_ref):
    cycle = Cycle(fluid_ref, Te=TE, Tc=TC, sh=5.0, sc=0.0)
    cycle.set_compressor(EffCompressor(fluid_ref, is_eff=0.8, vol_eff=0.9, N=1500 / 60, Vd=5e-5))
    cycle.set_condenser(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.set_evaporator(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    # The condenser sees the compressor discharge as its refrigerant inlet.
    assert cycle.condenser.state == "cond"
    assert cycle.condenser.Q == pytest.approx(cycle.m_flow * (cycle.h2 - cycle.h3))


def test_str_contains_all_states(eff_cycle):
    eff_cycle.calc(m_flow_sec_cond=0.1, p_sec_cond=1e5, T_sec_in_cond=20+273.15, m_flow_sec_evap=0.1, p_sec_evap=1e5, T_sec_in_evap=10+273.15)
    text = str(eff_cycle)
    assert "Cycle" in text
    for i in (1, 2, 3, 4):
        assert f"State {i}" in text
