from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

from refrigerationtoolbox.cycle.Compressor import Compressor

class PolynomialCompressor(Compressor):
    """Compressor modelled with manufacturer performance polynomials.

    Uses the ten-coefficient bi-cubic polynomials in the evaporation and condensation
    temperatures (AHRI 540 / EN 12900 form) to predict the mass flow rate, the cooling
    capacity and the input power directly from rating data. The heating capacity is the sum
    of power and cooling capacity, and the outlet enthalpy is derived from the heating
    capacity and the condenser-outlet state. The coefficients are expected in Kelvin; use
    :func:`~refrigerationtoolbox.cycle.Compressor.celsius_to_kelvin_coeffs` to convert
    Celsius-based rating polynomials.

    Args:
        fluid: CoolProp state object for the refrigerant.
        m_flow_coeffs: Ten polynomial coefficients for the mass flow rate [kg/s].
        Q_ref_coeffs: Ten polynomial coefficients for the cooling capacity [W].
        P_coeffs: Ten polynomial coefficients for the input power [W].
    """

    def __init__(self, fluid : AbstractState, m_flow_coeffs : np.ndarray, Q_ref_coeffs : np.ndarray,
                 P_coeffs : np.ndarray) -> None:
        super().__init__(fluid)

        self.m_flow_coeffs = m_flow_coeffs
        self.Q_ref_coeffs = Q_ref_coeffs
        self.P_coeffs = P_coeffs

    def calc(self, Te : float, Tc : float, sh : float, sc : float) -> None:
        """Solve the compressor from the rating polynomials.

        Evaluates the mass flow, cooling capacity and power polynomials at the given
        temperatures, forms the heating capacity as their sum, and derives the discharge
        enthalpy from the heating capacity and the (optionally subcooled) condenser-outlet
        state. The suction state is set at the evaporation pressure and superheated inlet
        temperature, and the discharge temperature/entropy/density follow from the discharge
        enthalpy at the condenser pressure.

        Args:
            Te: Evaporation temperature [K].
            Tc: Condensation temperature [K].
            sh: Superheat at the suction inlet [K].
            sc: Subcooling at the condenser outlet [K].
        """
        vars_poly = np.array([1, Te, Tc, Te**2, Te*Tc, Tc**2, Te**3, Te**2*Tc, Te*Tc**2, Tc**3])
        self.m_flow = np.dot(self.m_flow_coeffs, vars_poly).item()
        self.Q_ref = np.dot(self.Q_ref_coeffs, vars_poly).item()
        self.P = np.dot(self.P_coeffs, vars_poly).item()

        self.Q_heat = self.P + self.Q_ref
        q_heat = self.Q_heat / self.m_flow

        self.fluid.update(CP.QT_INPUTS, 0.0, Tc)
        p_c = self.fluid.p()
        h3 = self.fluid.hmass()
        # Only update the state using pressure and temperature if subcooling is significant to avoid
        # numerical issues with CoolProp.
        if sc > 1e-4:
            self.fluid.update(CP.PT_INPUTS, p_c, Tc - sc)
            h3 = self.fluid.hmass()
        
        self.h_out = h3 + q_heat

        self.T_in = Te + sh
        self.fluid.update(CP.QT_INPUTS, 1.0, Te)
        self.p_in = self.fluid.p()
        # Calculate the enthalpy and entropy at the evaporator pressure only if superheat is significant to avoid numerical issues with CoolProp
        if sh > 1e-4:
            self.fluid.update(CP.PT_INPUTS, self.p_in, Te + sh)

        self.h_in = self.fluid.hmass()
        self.s_in = self.fluid.smass()
        self.d_in = self.fluid.rhomass()

        self.fluid.update(CP.QT_INPUTS, 1.0, Tc)
        self.p_out = self.fluid.p()

        self.fluid.update(CP.HmassP_INPUTS, self.h_out, self.p_out)
        self.T_out = self.fluid.T()
        self.s_out = self.fluid.smass()
        self.d_out = self.fluid.rhomass()

    def __str__(self) -> str:
        """Return the base compressor summary extended with the heating and cooling capacities."""
        def fmt(val : float | None, scale : float = 1) -> str:
            return "N/A" if val is None else f"{val / scale:.2f}"

        return (
            f"{super().__str__()}\n"
            f"  Q_heat={fmt(self.Q_heat, 1e3)} kW, Q_ref={fmt(self.Q_ref, 1e3)} kW"
        )