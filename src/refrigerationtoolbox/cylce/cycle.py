import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

from refrigerationtoolbox.cylce.compressor import EffCompressor, PolynomialCompressor, celsius_to_kelvin_coeffs

class Cycle:
    def __init__(self, refrigerant: str, Te: float, Tc: float, sh=0.0, sc=0.0,
                 compressor=None, evaporator=None, condenser=None, expansion_valve=None):
        self.refrigerant = refrigerant
        self.fluid = AbstractState("HEOS", refrigerant)
        self.compressor = compressor
        self.evaporator = evaporator
        self.condenser = condenser
        self.expansion_valve = expansion_valve

        self.Te = Te      # K
        self.Tc = Tc      # K
        self.sh = sh      # K
        self.sc = sc      # K

    def set_compressor(self, compressor):
        self.compressor = compressor

    def set_evaporator(self, evaporator):
        self.evaporator = evaporator

    def set_condenser(self, condenser):
        self.condenser = condenser
    
    def set_expansion_valve(self, expansion_valve):
        self.expansion_valve = expansion_valve

    def calc(self):
        self.fluid.update(CP.QT_INPUTS, 1.0, self.Te)
        self.pe = self.fluid.p()
        self.fluid.update(CP.QT_INPUTS, 1.0, self.Tc)
        self.pc = self.fluid.p()

        # Calculate Q_heat, Q_ref, and P_comp for polynomial compressor, through polynomial
        # Implementation for Polynomial Compressor should not mix using polynomial heating and cooling capacity
        # with the actual enthalpy and mass flow rate calculations to fulfill conservation of energy.
        if type(self.compressor) == PolynomialCompressor:
            self.h_out_comp, self.m_flow, self.P_comp, self.Q_heat, self.Q_ref = self.compressor.calc(self.fluid, self.Te, self.Tc, self.sh, self.sc)
        else:
            self.h_out_comp, self.m_flow, self.P_comp = self.compressor.calc(self.fluid, self.Te, self.Tc, self.sh, self.sc)

        self._state1()
        self._state2()
        self._state3()
        self._state4()

        # Values are already calculated for polynomial compressor, so skip the calculation here
        if type(self.compressor) != PolynomialCompressor:
            self.Q_heat = (self.h2 - self.h3) * self.m_flow  # Heat rejected in the condenser
            self.Q_ref = (self.h1 - self.h4) * self.m_flow  # Heat absorbed in the evaporator

        self.COP_heat = self.Q_heat / self.P_comp  # Coefficient of performance for heating
        self.COP_cool = self.Q_ref / self.P_comp  # Coefficient of

    def _state1(self):
        self.T1 = self.Te + self.sh

        # Saturated vapor at evaporator outlet
        self.fluid.update(CP.QT_INPUTS, 1.0, self.Te)
        self.p1 = self.pe

        # Actual suction state
        # Only update state if superheat is significant to avoid numerical issues with CoolProp
        if self.sh > 1e-4:
            self.fluid.update(CP.PT_INPUTS, self.p1, self.T1)

        self.h1 = self.fluid.hmass()
        self.s1 = self.fluid.smass()
        self.d1 = self.fluid.rhomass()

    def _state2(self):
        self.p2 = self.pc
        self.h2 = self.h_out_comp
        self.fluid.update(CP.HmassP_INPUTS, self.h2, self.p2)
        self.T2 = self.fluid.T()
        self.s2 = self.fluid.smass()
        self.d2 = self.fluid.rhomass()

    def _state3(self):
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

    def _state4(self):
        self.p4 = self.pe
        # For polynomial compressor, the enthalpy at state 4 is calculated based on the
        # polynomial for the cooling capacity. Implementation for Polynomial Compressor should not mix using polynomial heating and cooling capacity
        # with the actual enthalpy calculations to fulfill conservation of energy.
        if type(self.compressor) == PolynomialCompressor:
            # Approximately isenthalp due to error through polynomial models
            self.h4 = self.h1 - (self.Q_ref / self.m_flow)
        else:
            # Isenthalp
            self.h4 = self.h3
        self.fluid.update(CP.HmassP_INPUTS, self.h4, self.p4)
        self.T4 = self.Te
        self.s4 = self.fluid.smass()
        self.d4 = self.fluid.rhomass()

    def __str__(self):
        str_rep = f"Cycle(refrigerant={self.refrigerant}, Te={self.Te}, Tc={self.Tc}, sh={self.sh}, sc={self.sc})"
        if hasattr(self, 'T4'):
            str_rep += "\n"
            str_rep += f"Q_heat: {self.Q_heat/1e3:.2f} kW, Q_ref: {self.Q_ref/1e3:.2f} kW, P_comp: {self.P_comp/1e3:.2f} kW COP_heat={self.COP_heat:.2f}, COP_cool= {self.COP_cool:.2f}\n"
            str_rep += f"State 1: T={self.T1:.2f}, P={self.p1/1e3:.2f} kPa, h={self.h1/1e3:.2f} kJ/kg, s={self.s1/1e3:.2f} kJ/(kgK), rho={self.d1:.2f} kg/m³\n"
            str_rep += f"State 2: T={self.T2:.2f}, P={self.p2/1e3:.2f} kPa, h={self.h2/1e3:.2f} kJ/kg, s={self.s2/1e3:.2f} kJ/(kgK), rho={self.d2:.2f} kg/m³\n"
            str_rep += f"State 3: T={self.T3:.2f}, P={self.p3/1e3:.2f} kPa, h={self.h3/1e3:.2f} kJ/kg, s={self.s3/1e3:.2f} kJ/(kgK), rho={self.d3:.2f} kg/m³\n"
            str_rep += f"State 4: T={self.T4:.2f}, P={self.p4/1e3:.2f} kPa, h={self.h4/1e3:.2f} kJ/kg, s={self.s4/1e3:.2f} kJ/(kgK), rho={self.d4:.2f} kg/m³"
        return str_rep

if __name__ == "__main__":


    cycle = Cycle(refrigerant="R134a", Te=273.15, Tc=303.15, sh=5.0, sc=0.0)
    # Bitzer software R143a ESH730Y Te = 0 °C, Tc = 30 °C, sh = 5K, sc = 0 K
    #m [kg/h];433,40559949511600000000;14,99099125357860000000;-0,33285851233968500000;0,20462358615376700000;0,00783475919685306000;-0,00550234077032464000;0,00147385448738919000;0,00003000181055104090;-0,00018073881737116100;-0,00000075930225497561;
    m_flow_coeffs=celsius_to_kelvin_coeffs(np.array([433.405599495116, 14.9909912535786, -0.332858512339685, 0.204623586153767, 0.00783475919685306, -0.00550234077032464, 0.00147385448738919, 3.00018105510409e-05, -0.000180738817371161, -7.5930225497561e-07])/3600)  # Convert from kg/h to kg/s
    #Q [W];24505,12988476420000000000;923,57374521494800000000;-189,01698961214700000000;14,47085894805140000000;-5,40337736007146000000;-0,16360259809141000000;0,08598465358143220000;-0,09085012430259930000;-0,01713808762391270000;-0,00020020016319336100;
    Q_ref_coeffs=celsius_to_kelvin_coeffs(np.array([24505.1298847642, 923.573745214948, -189.016989612147, 14.4708589480514, -5.40337736007146, -0.16360259809141, 0.0859846535814322, -0.0908501243025993, -0.0171380876239127, -0.000200200163193361]))
    #P[W]; 2096,19728039823000000000;15,11146738399150000000;33,27558403708940000000;-0,05338309873344150000;-0,19430059970364900000;0,43340928472529800000;-0,00897303325535081000;0,00530884169336271000;0,00001525882343664720;0,00607808111008199000;
    P_coeffs=celsius_to_kelvin_coeffs(np.array([2096.19728039823, 15.1114673839915, 33.2755840370894, -0.0533830987334415, -0.194300599703649, 0.433409284725298, -0.00897303325535081, 0.00530884169336271, 1.52588234366472e-05, 0.00607808111008199]))

    poly_comp = PolynomialCompressor(m_flow_coeffs=m_flow_coeffs, Q_ref_coeffs=Q_ref_coeffs, P_coeffs=P_coeffs)
    cycle.set_compressor(poly_comp)
    cycle.calc()
    print(cycle)

    eff_comp = EffCompressor(is_eff=0.8, vol_eff=0.9, N=1500/60, Vd=5e-5)
    cycle.set_compressor(eff_comp)
    cycle.calc()

    print(cycle)