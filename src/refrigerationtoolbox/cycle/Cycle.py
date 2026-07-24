import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.Compressor import Compressor
from refrigerationtoolbox.cycle.HeatExchanger import HeatExchanger
from refrigerationtoolbox.cycle.PolynomialCompressor import PolynomialCompressor


class Cycle:
    """Single-stage vapour-compression refrigeration / heat-pump cycle.

    Ties together a compressor and (optionally) an evaporator and a condenser and solves
    the four corner states of the cycle for a given evaporation and condensation
    temperature. The states are numbered in the direction of the refrigerant flow:

        1. compressor inlet / evaporator outlet (superheated suction gas),
        2. compressor outlet / condenser inlet (discharge gas),
        3. condenser outlet / expansion-valve inlet (subcooled liquid),
        4. expansion-valve outlet / evaporator inlet (two-phase mixture).

    Expansion is treated as isenthalpic (state 4 inherits the enthalpy of state 3), except
    for a :class:`PolynomialCompressor`, where state 4 is instead derived from the
    polynomial cooling capacity to stay consistent with the rating data. From the state
    enthalpies it computes the condenser and evaporator duties, the compressor power and
    the heating and cooling coefficients of performance. If heat-exchanger objects are
    attached they are solved as well so their sizing and pressure drops can be checked.

    Args:
        refrigerant: CoolProp state object for the working fluid.
        Te: Evaporation temperature [K].
        Tc: Condensation temperature [K].
        sh: Superheat at the evaporator outlet [K].
        sc: Subcooling at the condenser outlet [K].
        compressor: Compressor model; may also be set later via :meth:`set_compressor`.
        evaporator: Evaporator model; may also be set later via :meth:`set_evaporator`.
        condenser: Condenser model; may also be set later via :meth:`set_condenser`.
    """

    def __init__(
        self,
        refrigerant: AbstractState,
        Te: float,
        Tc: float,
        sh: float = 0.0,
        sc: float = 0.0,
        compressor: Compressor | None = None,
        evaporator: HeatExchanger | None = None,
        condenser: HeatExchanger | None = None,
    ) -> None:
        self.fluid = refrigerant
        self.compressor = compressor
        self.evaporator = evaporator
        self.condenser = condenser

        self.Te = Te  # K
        self.Tc = Tc  # K
        self.sh = sh  # K
        self.sc = sc  # K

        self.Q_heat = None
        self.Q_ref = None
        self.P_comp = None
        self.m_flow = None
        self.COP_heat = None
        self.COP_cool = None

        self.T1 = None
        self.p1 = None
        self.h1 = None
        self.s1 = None
        self.d1 = None

        self.T2 = None
        self.p2 = None
        self.h2 = None
        self.s2 = None
        self.d2 = None

        self.T3 = None
        self.p3 = None
        self.h3 = None
        self.s3 = None
        self.d3 = None

        self.T4 = None
        self.p4 = None
        self.h4 = None
        self.s4 = None
        self.d4 = None

    def set_compressor(self, compressor: Compressor) -> None:
        """Attach the compressor model used when solving the cycle."""
        self.compressor = compressor

    def set_evaporator(self, evaporator: HeatExchanger) -> None:
        """Attach the evaporator model used when solving the cycle."""
        self.evaporator = evaporator

    def set_condenser(self, condenser: HeatExchanger) -> None:
        """Attach the condenser model used when solving the cycle."""
        self.condenser = condenser

    def calc(
        self,
        m_flow_sec_cond: float,
        p_sec_cond: float,
        T_sec_in_cond: float,
        m_flow_sec_evap: float,
        p_sec_evap: float,
        T_sec_in_evap: float,
    ) -> None:
        """Solve the four cycle states and the performance figures for the current Te/Tc.

        Runs the compressor to get the mass flow, power and states 1 and 2, sets the
        condenser-outlet state 3 from the subcooling, and the evaporator-inlet state 4 from
        an isenthalpic expansion (or, for a :class:`PolynomialCompressor`, from its cooling
        capacity). It then evaluates the condenser and evaporator duties, the compressor
        power and the heating and cooling COPs, and, if attached, solves both heat
        exchangers with the corresponding secondary-side boundary conditions.

        Args:
            m_flow_sec_cond: Secondary (sink) mass flow rate through the condenser [kg/s].
            p_sec_cond: Secondary fluid pressure in the condenser [Pa].
            T_sec_in_cond: Secondary fluid inlet temperature to the condenser [K].
            m_flow_sec_evap: Secondary (source) mass flow rate through the evaporator [kg/s].
            p_sec_evap: Secondary fluid pressure in the evaporator [Pa].
            T_sec_in_evap: Secondary fluid inlet temperature to the evaporator [K].

        Raises:
            RuntimeError: If no compressor or no evaporator has been set.
        """
        if self.compressor is None:
            raise RuntimeError("self.compressor is None. Please set a compressor using cycle.set_compressor(...) or using the constructor")
        if self.evaporator is None:
            raise RuntimeError("self.evaporator is None. Please set an evaporator using cycle.set_evaporator(...) or using the constructor")

        self.fluid.update(CP.QT_INPUTS, 1.0, self.Te)
        self.pe = self.fluid.p()
        self.fluid.update(CP.QT_INPUTS, 1.0, self.Tc)
        self.pc = self.fluid.p()

        # Calculate Q_heat, Q_ref, and P_comp for polynomial compressor, through polynomial
        # Implementation for Polynomial Compressor should not mix using polynomial heating and cooling capacity
        # with the actual enthalpy and mass flow rate calculations to fulfill conservation of energy.
        self.compressor.calc(self.Te, self.Tc, self.sh, self.sc)
        self.m_flow = self.compressor.m_flow
        self.P_comp = self.compressor.P

        self.T1 = self.compressor.T_in
        self.p1 = self.compressor.p_in
        self.h1 = self.compressor.h_in
        self.s1 = self.compressor.s_in
        self.d1 = self.compressor.d_in

        self.T2 = self.compressor.T_out
        self.p2 = self.compressor.p_out
        self.h2 = self.compressor.h_out
        self.s2 = self.compressor.s_out
        self.d2 = self.compressor.d_out

        self.T3 = self.Tc - self.sc
        self.p3 = self.pc
        # Only update state using pressure and temperature if subcooling is significant to avoid numerical issues with CoolProp
        if self.sc > 1e-4:
            self.fluid.update(CP.PT_INPUTS, self.p3, self.T3)
        else:
            self.fluid.update(CP.QT_INPUTS, 0.0, self.T3)
        self.h3 = self.fluid.hmass()
        self.s3 = self.fluid.smass()
        self.d3 = self.fluid.rhomass()

        self.p4 = self.pe
        # For polynomial compressor, the enthalpy at state 4 is calculated based on the
        # polynomial for the cooling capacity. Implementation for Polynomial Compressor should not mix using polynomial heating and cooling capacity
        # with the actual enthalpy calculations to fulfill conservation of energy.
        if type(self.compressor) is PolynomialCompressor:
            # Approximately isenthalp due to error through polynomial models
            self.h4 = self.h1 - (self.compressor.Q_ref / self.compressor.m_flow)
        else:
            # Isenthalp
            self.h4 = self.h3

        self.fluid.update(CP.HmassP_INPUTS, self.h4, self.p4)
        self.T4 = self.Te
        self.s4 = self.fluid.smass()
        self.d4 = self.fluid.rhomass()

        # Values are already calculated for polynomial compressor, so skip the calculation here
        if type(self.compressor) is PolynomialCompressor:
            self.Q_heat = self.compressor.Q_heat
            self.Q_ref = self.compressor.Q_ref
        else:
            self.Q_heat = (self.h2 - self.h3) * self.m_flow  # Heat rejected in the condenser
            self.Q_ref = (self.h1 - self.h4) * self.m_flow  # Heat absorbed in the evaporator

        if self.condenser is not None:
            self.condenser.calc(
                m_flow_ref=self.m_flow,
                m_flow_sec=m_flow_sec_cond,
                h_ref_in=self.h2,
                h_ref_out=self.h3,
                p_ref=self.pc,
                p_sec=p_sec_cond,
                T_sec_in=T_sec_in_cond,
            )

        if self.evaporator is not None:
            self.evaporator.calc(
                m_flow_ref=self.m_flow,
                m_flow_sec=m_flow_sec_evap,
                h_ref_in=self.h4,
                h_ref_out=self.h1,
                p_ref=self.pe,
                p_sec=p_sec_evap,
                T_sec_in=T_sec_in_evap,
            )

        self.COP_heat = self.Q_heat / self.P_comp  # Coefficient of performance for heating
        self.COP_cool = self.Q_ref / self.P_comp  # Coefficient of

    def __str__(self) -> str:
        """Return a multi-line summary of the performance figures and the four cycle states."""

        def fmt(val: float | None, scale: float = 1) -> str:
            return "N/A" if val is None else f"{val / scale:.2f}"

        str_rep = f"Cycle (refrigerant={self.fluid.fluid_names()[0]}, Te={self.Te}, Tc={self.Tc}, sh={self.sh}, sc={self.sc})"
        str_rep += "\n"
        str_rep += f"   Q_heat: {fmt(self.Q_heat, 1e3)} kW, Q_ref: {fmt(self.Q_ref, 1e3)} kW, "
        str_rep += f"P_comp: {fmt(self.P_comp, 1e3)} kW, COP_heat={fmt(self.COP_heat)}, COP_cool= {fmt(self.COP_cool)}, m_flow={fmt(self.m_flow, 1)} kg/s\n"
        str_rep += f"   State 1: T={fmt(self.T1)}, P={fmt(self.p1, 1e3)} kPa, h={fmt(self.h1, 1e3)} kJ/kg, s={fmt(self.s1, 1e3)} kJ/(kgK), rho={fmt(self.d1, 1)} kg/m³\n"
        str_rep += f"   State 2: T={fmt(self.T2)}, P={fmt(self.p2, 1e3)} kPa, h={fmt(self.h2, 1e3)} kJ/kg, s={fmt(self.s2, 1e3)} kJ/(kgK), rho={fmt(self.d2, 1)} kg/m³\n"
        str_rep += f"   State 3: T={fmt(self.T3)}, P={fmt(self.p3, 1e3)} kPa, h={fmt(self.h3, 1e3)} kJ/kg, s={fmt(self.s3, 1e3)} kJ/(kgK), rho={fmt(self.d3, 1)} kg/m³\n"
        str_rep += f"   State 4: T={fmt(self.T4)}, P={fmt(self.p4, 1e3)} kPa, h={fmt(self.h4, 1e3)} kJ/kg, s={fmt(self.s4, 1e3)} kJ/(kgK), rho={fmt(self.d4, 1)} kg/m³"
        return str_rep
