"""Tests for the top-level ``Cycle`` orchestrator."""

import pytest
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.BasicHeatExchanger import BasicHeatExchanger
from refrigerationtoolbox.cycle.Cycle import Cycle
from refrigerationtoolbox.cycle.EffCompressor import EffCompressor
from refrigerationtoolbox.cycle.PolynomialCompressor import PolynomialCompressor

TE = 273.15
TC = 303.15
SH = 5.0
SC = 0.0


@pytest.fixture
def eff_cycle(fluid_ref):
    cycle = Cycle(fluid_ref, Te=TE, Tc=TC, sh=SH, sc=SC)
    cycle.set_compressor(EffCompressor(fluid_ref, is_eff=0.8, vol_eff=0.9, N=1500 / 60, Vd=5e-5))
    cycle.set_condenser(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.set_evaporator(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    return cycle


def test_verify_eff_cycle_no_sh_sc(eff_cycle):
    """ "Verify solution against TLK Energy log(p)-h diagram without superheat and subcooling
    for R134a, TE=0°C, TC=30°C, SH=0 K, SC=0 K https://tlk-energy.de/en/phase-diagrams/pressure-enthalpy"""
    eff_cycle.sh = 0
    eff_cycle.sc = 0
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    assert eff_cycle.h1 == pytest.approx(398600, rel=1e-3)
    assert eff_cycle.h2 == pytest.approx(423700, rel=1e-3)
    assert eff_cycle.h3 == pytest.approx(241700, rel=1e-3)
    assert eff_cycle.h4 == pytest.approx(241700, rel=1e-3)


def test_verify_eff_cycle_sh_sc(eff_cycle):
    """ "Verify solution against TLK Energy log(p)-h diagram without superheat and subcooling
    for R134a, TE=0°C, TC=30°C, SH=4 K, SC=6 K https://tlk-energy.de/en/phase-diagrams/pressure-enthalpy"""
    eff_cycle.sh = 4
    eff_cycle.sc = 6
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    assert eff_cycle.h1 == pytest.approx(402200, rel=1e-3)
    assert eff_cycle.h2 == pytest.approx(427800, rel=1e-3)
    assert eff_cycle.h3 == pytest.approx(233100, rel=1e-3)
    assert eff_cycle.h4 == pytest.approx(233100, rel=1e-3)


def test_calc_populates_all_states(eff_cycle):
    for i in (1, 2, 3, 4):
        assert getattr(eff_cycle, f"T{i}") is None
        assert getattr(eff_cycle, f"h{i}") is None
        assert getattr(eff_cycle, f"p{i}") is None
        assert getattr(eff_cycle, f"s{i}") is None
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    for i in (1, 2, 3, 4):
        assert getattr(eff_cycle, f"T{i}") is not None
        assert getattr(eff_cycle, f"h{i}") is not None
        assert getattr(eff_cycle, f"p{i}") is not None
        assert getattr(eff_cycle, f"s{i}") is not None


def test_pressure_levels_ordered(eff_cycle):
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    # Condensing pressure is above evaporating pressure.
    assert eff_cycle.pc > eff_cycle.pe
    assert eff_cycle.p2 == pytest.approx(eff_cycle.pc)
    assert eff_cycle.p1 == pytest.approx(eff_cycle.pe)


def test_cop_definitions(eff_cycle):
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    assert eff_cycle.COP_heat == pytest.approx(eff_cycle.Q_heat / eff_cycle.P_comp)
    assert eff_cycle.COP_cool == pytest.approx(eff_cycle.Q_ref / eff_cycle.P_comp)


def test_heating_cop_exceeds_cooling_cop_by_one(eff_cycle):
    """For an isenthalpic valve, Q_heat - Q_ref = P_comp, so COP_heat = COP_cool + 1."""
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    assert eff_cycle.COP_heat == pytest.approx(eff_cycle.COP_cool + 1.0, rel=1e-6)


def test_cop_values_physically_reasonable(eff_cycle):
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    assert eff_cycle.Q_heat > 0
    assert eff_cycle.Q_ref > 0
    assert eff_cycle.P_comp > 0
    assert eff_cycle.COP_cool > 1  # a working refrigeration cycle


def test_isenthalpic_expansion(eff_cycle):
    """With an efficiency compressor the valve is isenthalpic: h4 == h3."""
    eff_cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    assert eff_cycle.h4 == pytest.approx(eff_cycle.h3)


def test_polynomial_cycle_runs(fluid_ref, poly_coeffs):
    cycle = Cycle(fluid_ref, Te=TE, Tc=TC, sh=5.0, sc=0.0)
    cycle.set_compressor(PolynomialCompressor(fluid_ref, **poly_coeffs))
    cycle.set_condenser(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.set_evaporator(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.calc(
        m_flow_sec_cond=1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    # Polynomial duties are taken straight from the compressor map.
    assert cycle.Q_heat == pytest.approx(cycle.compressor.Q_heat)
    assert cycle.Q_ref == pytest.approx(cycle.compressor.Q_ref)
    assert cycle.COP_heat > 0


def test_calc_with_heat_exchangers(fluid_ref):
    cycle = Cycle(fluid_ref, Te=TE, Tc=TC, sh=5.0, sc=0.0)
    cycle.set_compressor(EffCompressor(fluid_ref, is_eff=0.8, vol_eff=0.9, N=1500 / 60, Vd=5e-5))
    cycle.set_condenser(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.set_evaporator(BasicHeatExchanger(fluid_ref, AbstractState("HEOS", "Water")))
    cycle.calc(
        m_flow_sec_cond=0.1,
        p_sec_cond=1e5,
        T_sec_in_cond=20 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )
    # The condenser sees the compressor discharge as its refrigerant inlet.
    assert cycle.condenser.state == "cond"
    assert cycle.condenser.Q == pytest.approx(cycle.m_flow * (cycle.h2 - cycle.h3))
