from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

from refrigerationtoolbox.cycle.Compressor import Compressor

class PolynomialCompressor(Compressor):
    def __init__(self, fluid : AbstractState, m_flow_coeffs : list, Q_ref_coeffs : list, P_coeffs : list):
        super().__init__(fluid)

        self.m_flow_coeffs = m_flow_coeffs
        self.Q_ref_coeffs = Q_ref_coeffs
        self.P_coeffs = P_coeffs

    def calc(self, Te : float, Tc : float, sh : float, sc : float):
        vars_poly = np.array([1, Te, Tc, Te**2, Te*Tc, Tc**2, Te**3, Te**2*Tc, Te*Tc**2, Tc**3])
        self.m_flow = np.dot(self.m_flow_coeffs, vars_poly).item()
        self.Q_ref = np.dot(self.Q_ref_coeffs, vars_poly).item()
        self.P = np.dot(self.P_coeffs, vars_poly).item()

        self.Q_heat = self.P + self.Q_ref
        q_heat = self.Q_heat / self.m_flow

        self.fluid.update(CP.QT_INPUTS, 0.0, Tc)
        h3 = self.fluid.hmass()
        if sc > 1e-4:
            self.fluid.update(CP.HmassT_INPUTS, h3, Tc - sc)
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

    def __str__(self):
        def fmt(val, scale=1):
            return "N/A" if val is None else f"{val / scale:.2f}"

        return (
            f"{super().__str__()}\n"
            f"  Q_heat={fmt(self.Q_heat, 1e3)} kW, Q_ref={fmt(self.Q_ref, 1e3)} kW"
        )