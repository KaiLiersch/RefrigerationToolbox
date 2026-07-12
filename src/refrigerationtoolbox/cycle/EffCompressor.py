import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

from refrigerationtoolbox.cycle.Compressor import Compressor

class EffCompressor(Compressor):
    def __init__(self, fluid : AbstractState, is_eff : float, vol_eff : float, N : float, Vd : float):
        super().__init__(fluid)
        self.is_eff = is_eff
        self.vol_eff = vol_eff
        self.N = N
        self.Vd = Vd

    def calc(self, Te : float, Tc : float, sh : float, sc : float):        
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

        self.fluid.update(CP.PSmass_INPUTS, self.p_out, self.s_in)
        # Calculate the isentropic enthalpy at the condenser pressure
        h_isentropic = self.fluid.hmass()

        # Calculate the actual enthalpy using the isentropic efficiency
        self.h_out = self.h_in + (h_isentropic - self.h_in) / self.is_eff
        self.fluid.update(CP.HmassP_INPUTS, self.h_out, self.p_out)
        self.T_out = self.fluid.T()
        self.s_out = self.fluid.smass()
        self.d_out = self.fluid.rhomass()

        self.m_flow = self.vol_eff * self.N * self.Vd * self.d_in
        self.P = self.m_flow * (self.h_out - self.h_in)
