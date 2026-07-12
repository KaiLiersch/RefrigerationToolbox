import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

from refrigerationtoolbox.cycle.HeatExchanger import HeatExchanger

class BasicHeatExchanger(HeatExchanger):
    def __init__(self, fluid_ref : AbstractState, fluid_sec : AbstractState):
        super().__init__(fluid_ref, fluid_sec)

    def calc(self, m_flow_ref : float, m_flow_sec : float,
                h_ref_in : float, h_ref_out : float, p_ref : float, p_sec : float, T_sec_in : float):
        super().calc(m_flow_ref, m_flow_sec, h_ref_in, h_ref_out, p_ref, p_sec, T_sec_in)

        # Caluclating Q and all Ts and hs for the input and output of the refrigeration and secondary flow
        self.Q = m_flow_ref * (h_ref_in - h_ref_out)

        self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_in, p_ref)
        self.T_ref_in = self.fluid_ref.T()
        self.s_ref_in = self.fluid_ref.smass()
        self.d_ref_in = self.fluid_ref.rhomass()
        self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_out, p_ref)
        self.T_ref_out = self.fluid_ref.T()
        self.s_ref_out = self.fluid_ref.smass()
        self.d_ref_out = self.fluid_ref.rhomass()
        self.fluid_sec.update(CP.PT_INPUTS, p_sec, T_sec_in)
        self.h_sec_in = self.fluid_sec.hmass()
        self.s_sec_in = self.fluid_sec.smass()
        self.d_sec_in = self.fluid_sec.rhomass()
        self.h_sec_out = self.h_sec_in + self.Q / self.m_flow_sec
        self.fluid_sec.update(CP.HmassP_INPUTS, self.h_sec_out, p_sec)
        self.T_sec_out = self.fluid_sec.T()
        self.s_sec_out = self.fluid_sec.smass()
        self.d_sec_out = self.fluid_sec.rhomass()
