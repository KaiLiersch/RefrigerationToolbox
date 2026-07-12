
from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

class Compressor(ABC):
    def __init__(self, fluid : AbstractState):
        self.fluid = fluid

        self.T_in = None
        self.p_in = None
        self.h_in = None
        self.s_in = None
        self.d_in = None

        self.T_out = None
        self.p_out = None
        self.h_out = None
        self.s_out = None
        self.d_out = None

        self.m_flow = None
        self.P = None

    @abstractmethod
    def calc(self, Te : float, Tc : float, sh : float, sc : float):
        pass
        
    def __str__(self):
        def fmt(val, scale=1):
            return "N/A" if val is None else f"{val / scale:.2f}"

        return (
            f"{self.__class__.__name__} [{self.fluid.fluid_names()[0]}] (m_flow = {fmt(self.m_flow)} kg/s, P = {fmt(self.P, 1e3)} kW)\n"
            f"  Inlet : T={fmt(self.T_in)} K, p={fmt(self.p_in, 1e3)} kPa, h={fmt(self.h_in, 1e3)} kJ/kg, s={fmt(self.s_in, 1e3)} kJ/kgK, d={fmt(self.d_in, 1)} kg/m³\n"
            f"  Outlet: T={fmt(self.T_out)} K, p={fmt(self.p_out, 1e3)} kPa, h={fmt(self.h_out, 1e3)} kJ/kg, s={fmt(self.s_out, 1e3)} kJ/kgK, d={fmt(self.d_out, 1)} kg/m³"
        )
    
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
    
