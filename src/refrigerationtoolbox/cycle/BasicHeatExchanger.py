import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.HeatExchanger import HeatExchanger


class BasicHeatExchanger(HeatExchanger):
    """Ideal two-stream heat exchanger with no geometry or pressure drop.

    The simplest concrete :class:`HeatExchanger`: it applies an energy balance only. The
    duty is fixed by the refrigerant enthalpy change, the secondary outlet enthalpy follows
    from that duty, and both streams keep their inlet pressure (no pressure drop, no sizing
    and no feasibility check). Useful as a lightweight stand-in when the detailed
    plate-exchanger model is not needed.
    """

    def __init__(self, fluid_ref: AbstractState, fluid_sec: AbstractState) -> None:
        super().__init__(fluid_ref, fluid_sec)

    def calc(
        self,
        m_flow_ref: float,
        m_flow_sec: float,
        h_ref_in: float,
        h_ref_out: float,
        p_ref: float,
        p_sec: float,
        T_sec_in: float,
    ) -> None:
        """Solve the exchanger from an energy balance alone.

        The duty is the refrigerant enthalpy change times its mass flow. The refrigerant
        inlet/outlet states follow from the given enthalpies and pressure; the secondary
        outlet enthalpy is obtained from the same duty and its outlet state resolved from
        it. Both streams keep their inlet pressure (no pressure drop is modelled).

        Args:
            m_flow_ref: Refrigerant mass flow rate [kg/s].
            m_flow_sec: Secondary fluid mass flow rate [kg/s].
            h_ref_in: Refrigerant specific enthalpy at the inlet [J/kg].
            h_ref_out: Refrigerant specific enthalpy at the outlet [J/kg].
            p_ref: Refrigerant pressure [Pa].
            p_sec: Secondary fluid pressure [Pa].
            T_sec_in: Secondary fluid inlet temperature [K].
        """
        super().calc(m_flow_ref, m_flow_sec, h_ref_in, h_ref_out, p_ref, p_sec, T_sec_in)

        # Duty and both streams' inlet/outlet states from an energy balance
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
