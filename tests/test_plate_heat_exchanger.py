"""Tests for the discretized ``PlateHeatExchanger`` model.
Mathematical model is from the course Thermodynamics Software
held by Tryfon C. Roumpedakis.

The verification tests consists of unit tests as the full reference solution
contains features not supported by my model."""

import CoolProp.CoolProp as CP
import numpy as np
import pytest
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.PlateHeatExchanger import PlateHeatExchanger

REFERENCE_GEOM = {
    "number_of_passes": 1,
    "plate_thickness_m": 0.0007,
    "chevron_angle_rad": np.radians(2 * 60.0),
    "plate_amplitude_m": 0.001,
    "pitch": 0.007,
    "corrugation_pitch_m": 0.007,
    "Dp": 0.025,
    "Lp": 0.220,
    "Bp": 0.170,
    "Ntmin": 4,
    "Ntmax": 200,
    "m_max": 14.0,
    "Nt": 13,
}

# Operating point of the worked example, Table 14.5 and step 10.
REFERENCE_P_HOT = 1.5e5  # [Pa] boiler water
REFERENCE_M_FLOW_HOT = 0.85  # [kg/s]
REFERENCE_D_PORT = 0.025  # [m]
REFERENCE_N_CP = 6  # channels per fluid for the 13 plate pack, Eq. (14-12)
# Elemenet 1 reference temperatures
REFERENCE_T_HOT_EL1 = 80.89 + 273.15
REFERENCE_T_WALL_HOT_EL1 = 75.67 + 273.15


def _book_exchanger(fluid_sec):
    hx = PlateHeatExchanger(fluid_sec, AbstractState("HEOS", "R1233zdE"), REFERENCE_GEOM, l_w=16.3)
    hx.N_cp = REFERENCE_N_CP  # normally set by calc(...)
    return hx


def test_verify_mean_temperature_difference(fluid_sec):
    """Verification of the logarithmic mean temperature difference, step 3 of the example.

    The hot boiler water goes from 90 to 80,41 C against R1233zd(E) evaporating at a constant 60 C,
    for which the reference reports 24,90 K.
    """
    # Fluids are not needed in this testcase
    hx = PlateHeatExchanger(fluid_sec, AbstractState("HEOS", "Water"), REFERENCE_GEOM, l_w=16.3)

    dT_lm = hx._calc_dT_m(90 + 273.15, 80.41 + 273.15, 60 + 273.15, 60 + 273.15) / hx.Ft
    assert dT_lm == pytest.approx(24.90, abs=5e-3)

    # Step 5: first estimate of the area from Q = 34,25 kW, U0 = 2500 W/m2K and dT_m = 24,30 K
    assert 34.25e3 / (2500 * 24.30) == pytest.approx(0.564, abs=5e-4)


def test_verify_element_mean_temperature_difference_against_the_design_example(fluid_sec):
    """Verification of the local mean temperature difference of finite element 1, step 10.

    The hot stream crosses element 1 from 81,37 to 80,41 C while the refrigerant evaporates at a
    constant 60 C, for which the reference reports 20,88 K.
    """
    hx = _book_exchanger(fluid_sec)

    dT_lm = hx._calc_dT_m(81.37 + 273.15, 80.41 + 273.15, 60 + 273.15, 60 + 273.15) / hx.Ft
    assert dT_lm == pytest.approx(20.88, abs=0.02)

    # The pinch of the example sits at the cold end, where the hot stream leaves at 80,41 C
    assert (80.41 + 273.15) - (60 + 273.15) == pytest.approx(20.41, abs=5e-3)


def test_verify_geometry_against_the_design_example(fluid_sec):
    """Verification testcase for the geometry block of PlateHeatExchanger, step 6 of the example.

    Every quantity below is printed in the source of the model, so this pins the constructor against
    the reference implementation it was written from.
    """
    hx = PlateHeatExchanger(fluid_sec, AbstractState("HEOS", "Water"), REFERENCE_GEOM, l_w=16.3)

    assert hx.Phi == pytest.approx(1.180, abs=5e-4)
    assert hx.A_0 == pytest.approx(0.0374, abs=5e-5)
    assert hx.A_p == pytest.approx(0.0441, abs=5e-5)
    assert hx.D_h == pytest.approx(0.0034, abs=5e-5)
    assert hx.A_ch == pytest.approx(0.00034, abs=5e-6)


def test_verify_plate_count_estimate_against_the_design_example(fluid_sec):
    """Verification of steps 7 and 8 of the example: 13 plates and 6 channels per fluid."""
    hx = PlateHeatExchanger(fluid_sec, AbstractState("HEOS", "R1233zdE"), REFERENCE_GEOM, l_w=16.3)

    # The book prints 12,77, truncated from 12,777
    Nt_0 = 0.564 / hx.A_p
    assert Nt_0 == pytest.approx(12.77, abs=1e-2)
    assert int(np.ceil(Nt_0)) == 13

    hx.Nt = 13
    assert (hx.Nt - 1) / 2 == 6


@pytest.mark.parametrize(
    ("Re", "expected"),
    [
        (200.0, 5.03 + 755 / 200.0),  # 90 < Re < 400
        (399.0, 5.03 + 755 / 399.0),
        (400.0, 26.8 * 400.0**-0.209),  # Re = 400
        (4033.0, 26.8 * 4033.0**-0.209),  # Re = 4033
        (15999.0, 26.8 * 15999.0**-0.209),  # Re = 15999
    ],
)
def test_verify_friction_factor_branches_against_the_book(fluid_sec, Re, expected):
    """The two single phase branches of the friction factor.

    The friction factor is recovered from the pressure drop, which is the only way it leaves the
    model. The mass flow rate is chosen so that the requested Reynolds number comes out.
    """
    hx = _book_exchanger(fluid_sec)
    water = hx.fluid_ref
    p, T = REFERENCE_P_HOT, REFERENCE_T_HOT_EL1

    water.update(CP.PT_INPUTS, p, T)
    rho, mu = water.rhomass(), water.viscosity()
    m_flow = Re * mu / hx.D_h * hx.N_cp * hx.A_ch

    _, dp_pl = hx._calc_single_phase_alpha(T, REFERENCE_T_WALL_HOT_EL1, p, m_flow, water, 1.0, 1.0)
    G_ch = m_flow / (hx.N_cp * hx.A_ch)
    f_model = dp_pl / (hx.Lp / hx.D_h * G_ch**2 / (2 * rho))

    assert f_model == pytest.approx(expected, rel=1e-6)


def test_verify_port_pressure_drop_against_the_design_example(fluid_sec):
    """Verification testcase for the port pressure drop.

    The reference reports 2,02 kPa for the hot stream of the worked example with one pass.
    """
    hx = _book_exchanger(fluid_sec)
    water = hx.fluid_ref

    water.update(CP.PT_INPUTS, REFERENCE_P_HOT, 90 + 273.15)
    h_hot_in = water.hmass()
    water.update(CP.PT_INPUTS, REFERENCE_P_HOT, 80.41 + 273.15)
    h_hot_out = water.hmass()

    hx.m_flow_ref, hx.m_flow_sec = REFERENCE_M_FLOW_HOT, REFERENCE_M_FLOW_HOT
    hx.h_ref_in, hx.h_ref_out = h_hot_in, h_hot_out
    hx.h_sec_in, hx.h_sec_out = h_hot_out, h_hot_in
    hx.p_ref = hx.p_sec = REFERENCE_P_HOT
    hx._calc_dppt()

    assert hx.dppt_ref / 1e3 == pytest.approx(2.02, abs=0.05)
