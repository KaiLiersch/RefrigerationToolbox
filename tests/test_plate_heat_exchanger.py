"""Tests for the discretized ``PlateHeatExchanger`` model."""

import pytest

from refrigerationtoolbox.cycle.PlateHeatExchanger import PlateHeatExchanger
from tests.conftest import COND_STATE


@pytest.fixture
def phex(fluid_ref, fluid_sec, plate_geom):
    return PlateHeatExchanger(fluid_ref, fluid_sec, plate_geom)


def test_geometry_derived_quantities(phex):
    # Derived areas/diameters must be positive.
    assert phex.A_p > 0
    assert phex.D_h > 0
    assert phex.A_ch > 0
    assert phex.Nt == 100


def test_calc_condenser_runs(phex):
    phex.calc(**COND_STATE)
    assert phex.state == "cond"


def test_heat_duty_uses_absolute_enthalpy_drop(phex):
    phex.calc(**COND_STATE)
    expected = COND_STATE["m_flow_ref"] * abs(COND_STATE["h_ref_in"] - COND_STATE["h_ref_out"])
    assert phex.Q == pytest.approx(expected)


def test_discretization_arrays_sized(phex):
    phex.calc(**COND_STATE)
    assert len(phex.phase_ref_el) == phex.num_elements
    assert len(phex.h_ref_nodes) == phex.num_elements + 1
    assert len(phex.h_sec_nodes) == phex.num_elements + 1


def test_pressure_drops_and_area_positive(phex):
    phex.calc(**COND_STATE)
    assert phex.dp_ref > 0
    assert phex.dp_sec > 0
    assert phex.A_phex > 0
    # Minimum approach temperature should be positive for a feasible condenser.
    assert phex.dT_min > 0


def test_str_reports_plate_count(phex):
    phex.calc(**COND_STATE)
    text = str(phex)
    assert "PlateHeatExchanger" in text
    assert "Nt = 100" in text
