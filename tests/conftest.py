"""Shared pytest fixtures for the refrigerationtoolbox test suite.

CoolProp ``AbstractState`` objects are stateful and are mutated in place during
``calc``. All fluid fixtures are therefore function scoped so every test gets a
fresh, un-mutated state.
"""

import numpy as np
import pytest
from CoolProp.CoolProp import AbstractState


@pytest.fixture
def fluid_ref():
    """Fresh R134a refrigerant state."""
    return AbstractState("HEOS", "R134a")


@pytest.fixture
def fluid_sec():
    """Fresh Water secondary-side state."""
    return AbstractState("HEOS", "Water")


@pytest.fixture
def plate_geom():
    """Plate heat exchanger geometry used in the project examples."""
    return {
        "number_of_passes": 1,
        "plate_thickness_m": 0.0007,
        "chevron_angle_rad": 1.0471975511965976,
        "pitch": 0.0025,
        "plate_amplitude_m": 0.001,
        "corrugation_pitch_m": 0.007,
        "Dp": 0.023,
        "Lp": 0.25,
        "Bp": 0.113,
        "Ntmin": 4,
        "Ntmax": 150,
        "m_max": 14.0,
        "Nt": 100,
    }


# Bitzer ESH730Y polynomial coefficients (already converted from Celsius to
# Kelvin), matching the example in Cycle.py.
@pytest.fixture
def poly_coeffs():
    from refrigerationtoolbox.cycle.Compressor import celsius_to_kelvin_coeffs

    m_flow_coeffs = celsius_to_kelvin_coeffs(
        np.array([
            433.405599495116, 14.9909912535786, -0.332858512339685,
            0.204623586153767, 0.00783475919685306, -0.00550234077032464,
            0.00147385448738919, 3.00018105510409e-05, -0.000180738817371161,
            -7.5930225497561e-07,
        ]) / 3600.0
    )
    Q_ref_coeffs = celsius_to_kelvin_coeffs(
        np.array([
            24505.1298847642, 923.573745214948, -189.016989612147,
            14.4708589480514, -5.40337736007146, -0.16360259809141,
            0.0859846535814322, -0.0908501243025993, -0.0171380876239127,
            -0.000200200163193361,
        ])
    )
    P_coeffs = celsius_to_kelvin_coeffs(
        np.array([
            2096.19728039823, 15.1114673839915, 33.2755840370894,
            -0.0533830987334415, -0.194300599703649, 0.433409284725298,
            -0.00897303325535081, 0.00530884169336271, 1.52588234366472e-05,
            0.00607808111008199,
        ])
    )
    return {
        "m_flow_coeffs": m_flow_coeffs,
        "Q_ref_coeffs": Q_ref_coeffs,
        "P_coeffs": P_coeffs,
    }


# Reference operating point for a single condenser/evaporator ``calc`` call.
# h_ref_in > h_ref_out => condenser (cooling/condensing the refrigerant).
COND_STATE = dict(
    m_flow_ref=0.01,
    m_flow_sec=0.1,
    h_ref_in=433840.0,
    h_ref_out=241720.0,
    p_ref=770200.0,
    p_sec=100000.0,
    T_sec_in=20 + 273.15,
)
