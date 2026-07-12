from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

class HeatExchanger(ABC):
    def __init__(self, fluid_ref : AbstractState, fluid_sec : AbstractState):
        self.fluid_ref = fluid_ref
        self.fluid_sec = fluid_sec

        self.state = "Undefined"

        self.T_ref_in = None
        self.h_ref_in = None
        self.s_ref_in = None
        self.d_ref_in = None
        self.p_ref = None

        self.T_ref_out = None
        self.h_ref_out = None
        self.s_ref_out = None
        self.d_ref_out = None

        self.T_sec_in = None
        self.h_sec_in = None
        self.s_sec_in = None
        self.d_sec_in = None
        self.p_sec = None

        self.T_sec_out = None
        self.h_sec_out = None
        self.s_sec_out = None
        self.d_sec_out = None
        self.f_sec_out = None

        self.m_flow_ref = None
        self.m_flow_sec = None
        self.Q = None

    @abstractmethod
    def calc(self, m_flow_ref : float, m_flow_sec : float, h_ref_in : float,
                h_ref_out : float, p_ref :float, p_sec : float, T_sec_in : float):
        self.m_flow_ref = m_flow_ref
        self.m_flow_sec = m_flow_sec
        self.h_ref_in = h_ref_in
        self.h_ref_out = h_ref_out
        self.p_ref = p_ref
        self.p_sec = p_sec
        self.T_sec_in = T_sec_in

        self.state = "cond" if h_ref_in > h_ref_out else "evap"

    def __str__(self):
        def fmt(val, scale=1):
            return "N/A" if val is None else f"{val / scale:.2f}"

        return (
            f"{self.__class__.__name__} ({self.state}, Q = {fmt(self.Q, 1e3)} kW)\n"
            f"  Refrigerant [{self.fluid_ref.fluid_names()[0]}] (m_flow = {fmt(self.m_flow_ref)} kg/s p={fmt(self.p_ref, 1e3)} kPa):\n"
            f"    Inlet : T={fmt(self.T_ref_in)} K, h={fmt(self.h_ref_in, 1e3)} kJ/kg, s={fmt(self.s_ref_in, 1e3)} kJ/kgK, d={fmt(self.d_ref_in, 1)} kg/m³\n"
            f"    Outlet: T={fmt(self.T_ref_out)} K, h={fmt(self.h_ref_out, 1e3)} kJ/kg, s={fmt(self.s_ref_out, 1e3)} kJ/kgKd={fmt(self.d_ref_out, 1)} kg/m³\n"
            f"  Secondary fluid [{self.fluid_sec.fluid_names()[0]}] (m_flow = {fmt(self.m_flow_sec)} kg/s, p={fmt(self.p_sec, 1e3)} kPa):\n"
            f"    Inlet : T={fmt(self.T_sec_in)} K, h={fmt(self.h_sec_in, 1e3)} kJ/kg, s={fmt(self.s_sec_in, 1e3)} kJ/kgK, d={fmt(self.d_sec_in, 1)} kg/m³\n"
            f"    Outlet: T={fmt(self.T_sec_out)} K, h={fmt(self.h_sec_out, 1e3)} kJ/kg, s={fmt(self.s_sec_out, 1e3)} kJ/kgK, d={fmt(self.d_sec_out, 1)} kg/m³"
        )
