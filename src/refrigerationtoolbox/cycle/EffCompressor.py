import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.Compressor import Compressor


class EffCompressor(Compressor):
    """Compressor modelled through isentropic and volumetric efficiencies.

    The outlet state is obtained by first computing the isentropic (constant-entropy)
    compression to the condenser pressure and then correcting it with the isentropic
    efficiency to get the real enthalpy rise. The mass flow rate follows from the
    volumetric efficiency, the running speed, the displacement volume and the suction-gas
    density, and the power from the enthalpy rise times the mass flow.

    Args:
        fluid: CoolProp state object for the refrigerant.
        is_eff: Isentropic efficiency [-].
        vol_eff: Volumetric efficiency [-].
        N: Rotational speed [1/s].
        Vd: Displacement volume per revolution [m³].
    """

    def __init__(self, fluid: AbstractState, is_eff: float, vol_eff: float, N: float, Vd: float) -> None:
        super().__init__(fluid)
        self.is_eff = is_eff
        self.vol_eff = vol_eff
        self.N = N
        self.Vd = Vd

    def calc(self, Te: float, Tc: float, sh: float, sc: float) -> None:
        """Solve the compressor using the isentropic and volumetric efficiencies.

        The suction state is set at the evaporation pressure and the superheated inlet
        temperature. The discharge enthalpy is the isentropic value at the condenser
        pressure corrected by the isentropic efficiency. The mass flow follows from the
        volumetric efficiency, speed, displacement and suction density, and the power from
        the enthalpy rise. The subcooling ``sc`` is not used by this model.

        Args:
            Te: Evaporation temperature [K].
            Tc: Condensation temperature [K].
            sh: Superheat at the suction inlet [K].
            sc: Subcooling at the condenser outlet [K] (unused here).
        """
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
