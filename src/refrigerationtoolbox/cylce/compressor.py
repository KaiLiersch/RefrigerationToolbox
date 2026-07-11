
from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

class Compressor(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def calc(self, fluid : AbstractState, Te : float, Tc : float, sh : float, sc : float):
        pass


class EffCompressor(Compressor):
    def __init__(self, is_eff : float, vol_eff : float, N : float, Vd : float):
        super().__init__()
        self.is_eff = is_eff
        self.vol_eff = vol_eff
        self.N = N
        self.Vd = Vd

    def calc(self, fluid : AbstractState, Te : float, Tc : float, sh : float, sc : float):        
        # Calculate the isentropic enthalpy at the condenser pressure
        fluid.update(CP.QT_INPUTS, 1.0, Te)
        pe = fluid.p()
        # Calculate the enthalpy and entropy at the evaporator pressure only if superheat is significant to avoid numerical issues with CoolProp
        if sh > 1e-4:
            fluid.update(CP.PT_INPUTS, pe, Te + sh)
        h_in = fluid.hmass()
        s_in = fluid.smass()
        d_in = fluid.rhomass()
        
        fluid.update(CP.QT_INPUTS, 1.0, Tc)
        pc = fluid.p()

        fluid.update(CP.PSmass_INPUTS, pc, s_in)
        h_isentropic = fluid.hmass()

        # Calculate the actual enthalpy using the isentropic efficiency
        h_out = h_in + (h_isentropic - h_in) / self.is_eff

        m_flow = self.vol_eff * self.N * self.Vd * d_in
        P = m_flow * (h_out - h_in)
        return h_out, m_flow, P

class PolynomialCompressor(Compressor):
    def __init__(self, m_flow_coeffs : list, Q_ref_coeffs : list, P_coeffs : list):
        super().__init__()

        self.m_flow_coeffs = m_flow_coeffs
        self.Q_ref_coeffs = Q_ref_coeffs
        self.P_coeffs = P_coeffs

    def calc(self, fluid : AbstractState, Te : float, Tc : float, sh : float, sc : float):
        vars_poly = np.array([1, Te, Tc, Te**2, Te*Tc, Tc**2, Te**3, Te**2*Tc, Te*Tc**2, Tc**3])
        m_flow = np.dot(self.m_flow_coeffs, vars_poly).item()
        Q_ref = np.dot(self.Q_ref_coeffs, vars_poly).item()
        P = np.dot(self.P_coeffs, vars_poly).item()

        Q_heat = P + Q_ref
        q_heat = Q_heat / m_flow

        print(Q_heat, Q_ref, P, m_flow, q_heat)

        fluid.update(CP.QT_INPUTS, 0.0, Tc)
        hc = fluid.hmass()
        if sc > 1e-4:
            fluid.update(CP.HmassT_INPUTS, hc, Tc - sc)
            hc = fluid.hmass()
        
        h_out = hc + q_heat
        return h_out, m_flow, P, Q_heat, Q_ref

    

def celsius_to_kelvin_coeffs(coeffs):
    """
    Convert AHRI540 / EN12900 10-coefficient polynomial
    from Celsius to Kelvin.

    Coefficient order:
        1,
        Te,
        Tc,
        Te²,
        Te Tc,
        Tc²,
        Te³,
        Te² Tc,
        Te Tc²,
        Tc³
    """

    K = 273.15

    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = coeffs

    b1 = (a1 - K*(a2 + a3) + K**2*(a4 + a5 + a6) - K**3*(a7 + a8 + a9 + a10))
    b2 = (a2 - 2*K*a4 - K*a5 + 3*K**2*a7 + 2*K**2*a8 + K**2*a9)
    b3 = (a3 - K*a5 - 2*K*a6 + K**2*a8 + 2*K**2*a9 + 3*K**2*a10)
    b4 = (a4 - 3*K*a7 - K*a8)
    b5 = (a5 - 2*K*a8 - 2*K*a9)
    b6 = (a6 - K*a9 - 3*K*a10)
    b7 = a7
    b8 = a8
    b9 = a9
    b10 = a10

    return np.array([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10])
    
