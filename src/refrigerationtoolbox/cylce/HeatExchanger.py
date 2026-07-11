from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
import numpy as np

class HeatExchanger(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def calc(self, fluid_ref : AbstractState, fluid_cool : AbstractState, m_flow_ref : float, m_flow_sec : float,
                h_ref_in, h_ref_out, p_ref, p_sec, T_sec_in):
        pass

class PlateHeatExchanger(HeatExchanger):
    def __init__(self, geom : dict, num_elements=10, Ft=0.9, Rf_ref=0.0001 , Rf_sec=0.00003, l_w=20, dp_ref_max=2e4, dp_sec_max=2e4):
        super().__init__()
        self.num_elements = num_elements
        self.Ft = Ft
        self.Rf_ref = Rf_ref
        self.Rf_sec = Rf_sec
        self.l_w = l_w
        self.dp_ref_max = dp_ref_max
        self.dp_sec_max = dp_sec_max

        self.Ntmax = geom["Ntmax"]
        self.phi = geom["chevron_angle_rad"] / 2 # lower case phi - angle against main flow direction
        self.a = geom["plate_amplitude_m"]
        self.L = geom["pitch"]
        self.t = geom["plate_thickness_m"]
        self.Nt = geom["Nt"]
        self.X = 2 * np.pi * self.a / self.L
        self.Phi = 1 / 6 * (1 + np.sqrt(1+self.X**2) + 4 * np.sqrt(1 + 0.5*self.X**2)) 

        self.Dp = geom["Dp"]
        self.Bp = geom["Bp"]
        self.Lp = geom["Lp"]
        self.A_0 = self.Lp * self.Bp
        self.A_p = self.Phi * self.A_0
        self.D_h = 4 * self.a / self.Phi
        self.A_ch = 2 * self.a * self.Bp

    def calc(self, fluid_ref : AbstractState, fluid_sec : AbstractState, m_flow_ref : float, m_flow_sec : float,
                h_ref_in, h_ref_out, p_ref, p_sec, T_sec_in):
        self.fluid_ref = fluid_ref
        self.fluid_sec = fluid_sec
        self.h_ref_in = h_ref_in
        self.h_ref_out = h_ref_out
        self.p_ref = p_ref
        self.p_sec = p_sec
        self.T_sec_in = T_sec_in
        self.m_flow_ref = m_flow_ref
        self.m_flow_sec = m_flow_sec

        self.ref_fluid = "R134a"
        self.sec_fluid = "Water"

        self.state = "cond" if h_ref_in > h_ref_out else "evap"

        # Caluclating Q and all Ts and hs for the input and output of the refrigeration and secondary flow
        self.Q = m_flow_ref * abs(h_ref_in - h_ref_out)

        self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_in, p_ref)
        self.T_ref_in = self.fluid_ref.T()
        self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_out, p_ref)
        self.T_ref_out = self.fluid_ref.T()
        self.fluid_sec.update(CP.PT_INPUTS, p_sec, T_sec_in)
        self.h_sec_in = self.fluid_sec.hmass()
        self.h_sec_out = self.h_sec_in + self.Q / self.m_flow_sec
        self.fluid_sec.update(CP.HmassP_INPUTS, self.h_sec_out, p_sec)
        self.T_sec_out = self.fluid_sec.T()

        dT_m = self._calc_dT_m(self.T_ref_in, self.T_ref_out, self.T_sec_in, self.T_sec_out)

        self.N_cp = (self.Nt - 1) / 2 
        if self.N_cp < 1: self.N_cp = 1

        A_available = self.Nt * self.A_p

        # Generate discretization with non uniform elements
        # Elements conform to boundaries of phase regions. E. g. an element can not contain both a superheated gas and two phase region
        if self.state == "cond":
            phase_ref_el, h_ref_nodes, h_sec_nodes = self._discretize_condensor()
        elif self.state == "evap":
            phase_ref_el, h_ref_nodes, h_sec_nodes = self._discretize_evaporator()

        # Solve on discretization
        A_phex, dppl_ref, dppl_sec, dT_min = self._calc_fe(phase_ref_el, h_ref_nodes, h_sec_nodes)

        # Calculate pressure drop in heat exchanger
        dppt_ref, dppt_sec = self._calc_dppt()

        self.dp_ref = dppt_ref + dppl_ref
        self.dp_sec = dppt_sec + dppl_sec

        # Check if heat exchanger is sufficiently sized
        if A_phex <= A_available:
            # Check if pressure drops are sufficiently small
            self.feas = True

        sol_dict = {
            "Q": self.Q,
            "dp_ref": self.dp_ref,
            "dp_sec": self.dp_sec,
            "Nt": self.Nt,
            "dT_min": dT_min,
            "T_Sec_out": self.T_sec_out
        }
        return sol_dict
    
    def _discretize_condensor(self):
        # find two phase region
        self.fluid_ref.update(CP.PQ_INPUTS, self.p_ref, 0.0)
        h_sat0_ref = self.fluid_ref.hmass()
        self.fluid_ref.update(CP.PQ_INPUTS, self.p_ref, 1.0)
        h_sat1_ref = self.fluid_ref.hmass()

        phase_ref_el = np.zeros(self.num_elements)
        h_ref_nodes = np.zeros(self.num_elements + 1)
        h_sec_nodes = np.zeros(self.num_elements + 1)

        # Discretization placing nodes on phase boundaries to to increase accuracy
        tot_diff_h_ref = self.h_ref_in - self.h_ref_out
        tot_diff_h_sec = self.h_sec_out - self.h_sec_in
        # Calculate fractions of each phase region
        enthalpy_frac_gas = (self.h_ref_in - h_sat1_ref) / tot_diff_h_ref 
        enthalpy_frac_liquid = (h_sat0_ref - self.h_ref_out) / tot_diff_h_ref
        enthaply_frac_vle = 1 - enthalpy_frac_gas - enthalpy_frac_liquid

        # Calculate the number of elements for each phase region
        num_el_gas = int(np.ceil(enthalpy_frac_gas * self.num_elements))
        num_el_liquid = int(np.ceil(enthalpy_frac_liquid * self.num_elements))
        num_el_vle = int(self.num_elements - num_el_gas - num_el_liquid)
        # Save which phase region each element has
        phase_ref_el[0:num_el_gas] = 0
        phase_ref_el[num_el_gas:num_el_gas+num_el_vle] = 1
        phase_ref_el[self.num_elements-num_el_liquid:] = 2

        # Calculate enthalpy step sizes for each phase region for refrigeration fluid
        dh_ref_gas = (self.h_ref_in - h_sat1_ref) / num_el_gas
        dh_ref_vle = (h_sat1_ref - h_sat0_ref) / num_el_vle
        dh_ref_liquid = (h_sat0_ref - self.h_ref_out) / num_el_liquid

        # Calculate enthalpy step size for each phase region for each secondary fluid
        dh_sec_1 = dh_ref_gas / tot_diff_h_ref * tot_diff_h_sec
        dh_sec_2 = dh_ref_vle / tot_diff_h_ref * tot_diff_h_sec
        dh_sec_3 = dh_ref_liquid / tot_diff_h_ref * tot_diff_h_sec

        # Calculate the enthalpies at each node depending on the phase region
        h_ref_nodes[0] = self.h_ref_in
        h_sec_nodes[0] = self.h_sec_out
        # The first few elements always contain superheated gas
        for i in range(1, num_el_gas + 1):
            h_ref_nodes[i] = h_ref_nodes[i-1] - dh_ref_gas    
            h_sec_nodes[i] = h_sec_nodes[i-1] - dh_sec_1
        # Then two phase region elements
        for i in range(num_el_gas + 1, num_el_gas + num_el_vle + 1):
            h_ref_nodes[i] = h_ref_nodes[i-1] - dh_ref_vle
            h_sec_nodes[i] = h_sec_nodes[i-1] - dh_sec_2
        # Finally, subcooled liquid elements
        for i in range(num_el_gas + num_el_vle + 1, self.num_elements + 1):
            h_ref_nodes[i] = h_ref_nodes[i-1] - dh_ref_liquid
            h_sec_nodes[i] = h_sec_nodes[i-1] - dh_sec_3

        return phase_ref_el, h_ref_nodes, h_sec_nodes

    def _discretize_evaporator(self):
        # Not required in exercise 5

        # find two phase region
        self.fluid_ref.update(CP.PQ_INPUTS, self.p_ref, 0.0)
        h_sat0_ref = self.fluid_ref.hmass()
        self.fluid_ref.update(CP.PQ_INPUTS, self.p_ref, 1.0)
        h_sat1_ref = self.fluid_ref.hmass()

        phase_ref_el = np.zeros(self.num_elements)
        h_ref_nodes = np.zeros(self.num_elements + 1)
        h_sec_nodes = np.zeros(self.num_elements + 1)


        # Discretization placing nodes on phase boundaries to to increase accuracy
        tot_diff_h_ref = self.h_ref_out - self.h_ref_in
        tot_diff_h_sec = self.h_sec_in - self.h_sec_out
        enthalpy_frac_gas = (self.h_ref_out - h_sat1_ref) / tot_diff_h_ref 
        enthalpy_frac_liquid = (h_sat0_ref - self.h_ref_in) / tot_diff_h_ref
        if enthalpy_frac_gas < 0: enthalpy_frac_gas = 0
        if enthalpy_frac_liquid < 0: enthalpy_frac_liquid = 0
        enthaply_frac_vle = 1 - enthalpy_frac_gas - enthalpy_frac_liquid

        num_el_gas = int(np.ceil(enthalpy_frac_gas * self.num_elements))
        num_el_liquid = int(np.ceil(enthalpy_frac_liquid * self.num_elements))
        num_el_vle = int(self.num_elements - num_el_gas - num_el_liquid)
        phase_ref_el[0:num_el_liquid] = 2
        phase_ref_el[num_el_liquid:num_el_liquid+num_el_vle] = 1
        phase_ref_el[self.num_elements-num_el_gas:] = 0

        dh_ref_gas = 0 if num_el_gas == 0 else (self.h_ref_out - h_sat1_ref) / num_el_gas
        dh_ref_vle = 0 if num_el_vle == 0 else (h_sat1_ref - max(h_sat0_ref, self.h_ref_in)) / num_el_vle
        dh_ref_liquid = 0 if num_el_liquid == 0 else (h_sat0_ref - self.h_ref_in) / num_el_liquid

        dh_sec_1 = dh_ref_gas / tot_diff_h_ref * tot_diff_h_sec
        dh_sec_2 = dh_ref_vle / tot_diff_h_ref * tot_diff_h_sec
        dh_sec_3 = dh_ref_liquid / tot_diff_h_ref * tot_diff_h_sec

        h_ref_nodes[0] = self.h_ref_in
        h_sec_nodes[0] = self.h_sec_out
        for i in range(1, num_el_liquid + 1):
            h_ref_nodes[i] = h_ref_nodes[i-1] + dh_ref_liquid    
            h_sec_nodes[i] = h_sec_nodes[i-1] + dh_sec_3
        for i in range(num_el_liquid + 1, num_el_liquid + num_el_vle + 1):
            h_ref_nodes[i] = h_ref_nodes[i-1] + dh_ref_vle
            h_sec_nodes[i] = h_sec_nodes[i-1] + dh_sec_2
        for i in range(num_el_liquid + num_el_vle + 1, self.num_elements + 1):
            h_ref_nodes[i] = h_ref_nodes[i-1] + dh_ref_gas
            h_sec_nodes[i] = h_sec_nodes[i-1] + dh_sec_1

        return phase_ref_el, h_ref_nodes, h_sec_nodes

    def _calc_fe(self, phase_ref_el, h_ref_nodes, h_sec_nodes):
        # Define element arrays
        A_phex_elements = np.zeros(self.num_elements)
        dppl_ref_elements = np.zeros(self.num_elements)
        dppl_sec_elements = np.zeros(self.num_elements)
        dT_elements = np.zeros(self.num_elements)

        # Total enthalpy difference between inlet and outlet of refrigeration fluid
        tot_dh = abs(h_ref_nodes[0] - h_ref_nodes[-1]) 

        a_ref_elements = np.zeros(self.num_elements)
        a_sec_elements = np.zeros(self.num_elements)

        # Loop over all elements
        for i, phase in enumerate(phase_ref_el):
            # Enthalpy of element is the average of enthalpy at nodes
            h_ref_el = (h_ref_nodes[i] + h_ref_nodes[i+1]) / 2
            h_sec_el = (h_sec_nodes[i] + h_sec_nodes[i+1]) / 2

            # Calculations are based on the enthalpies at the elements center
            T_ref_el = CP.PropsSI("T", "H", h_ref_el, "P", self.p_ref, self.ref_fluid)
            T_sec_el = CP.PropsSI("T", "H", h_sec_el, "P", self.p_sec, self.sec_fluid)
            T_plate_el = (T_ref_el + T_sec_el) / 2
            Tw_ref_el = (T_ref_el + T_plate_el) / 2
            Tw_sec_el = (T_sec_el + T_plate_el) / 2

            T0_ref = CP.PropsSI("T", "H", h_ref_nodes[i], "P", self.p_ref, self.ref_fluid)
            T1_ref = CP.PropsSI("T", "H", h_ref_nodes[i+1], "P", self.p_ref, self.ref_fluid)
            T0_sec = CP.PropsSI("T", "H", h_sec_nodes[i], "P", self.p_sec, self.sec_fluid)
            T1_sec = CP.PropsSI("T", "H", h_sec_nodes[i+1], "P", self.p_sec, self.sec_fluid)

            # Calculation of temperature difference between hot and cold side
            # Takes minimum value of difference of both nodes and in the center of the element
            if self.state == "cond":
                dT_elements[i] = min([T_ref_el - T_sec_el, T0_ref - T0_sec, T1_ref - T1_sec])
            if self.state == "evap":
                dT_elements[i] = min([T_sec_el - T_ref_el, T0_sec - T0_ref, T1_sec - T1_ref])

            # Enthalpy change of the element
            dh_ref = abs(h_ref_nodes[i] - h_ref_nodes[i + 1])
            dh_sec = abs(h_sec_nodes[i] - h_sec_nodes[i + 1])

            # Secondary fluid is always a single phase fluid
            a_sec, dppl_sec_el = self._calc_single_phase_alpha(T_sec_el, Tw_sec_el, self.p_sec, self.m_flow_sec, self.fluid_sec, dh_ref, tot_dh)
            
            # Refrigeration fluid might be two phase or single phase
            # 0 - gas; 1 - vle; 2 - liquid
            if phase == 0 or phase ==2:
                a_ref, dppl_ref_el = self._calc_single_phase_alpha(T_ref_el, Tw_ref_el, self.p_ref, self.m_flow_ref, self.fluid_ref, dh_ref, tot_dh, h_ref_el)
            elif phase == 1:
                if self.state == "cond":
                    a_ref, dppl_ref_el = self._calc_condensation_alpha(T_ref_el, self.p_ref, h_ref_el, self.m_flow_ref, self.fluid_ref, dh_ref, tot_dh)
                elif self.state == "evap":
                    # Not required for this exercise
                    a_ref, dppl_ref_el = self._calc_evaporation_alpha(T_ref_el, self.p_ref, h_ref_el, self.m_flow_ref, self.fluid_ref, dh_ref, tot_dh)
            else:
                print("WARNING: The phase of the hot fluid could not be properly determinded!")
            a_ref_elements[i] = a_ref
            a_sec_elements[i] = a_sec

            dppl_ref_elements[i] = dppl_ref_el
            dppl_sec_elements[i] = dppl_sec_el

            # Calculating U, A_phex for each element
            U_el = 1 / (1 / a_ref + 1 / a_sec + self.Rf_ref + self.Rf_sec + self.t / self.l_w)

            Q_ref = dh_ref * self.m_flow_ref
            Q_sec = dh_sec * self.m_flow_sec

            self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_nodes[i], self.p_ref)
            T0_ref = self.fluid_ref.T()
            self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_nodes[i+1], self.p_ref)
            T1_ref = self.fluid_ref.T()
            self.fluid_sec.update(CP.HmassP_INPUTS, h_sec_nodes[i], self.p_sec)
            T0_sec = self.fluid_sec.T()
            self.fluid_sec.update(CP.HmassP_INPUTS, h_sec_nodes[i+1], self.p_sec)
            T1_sec = self.fluid_sec.T()

            dT_m_i = self._calc_dT_m(T0_ref, T1_ref, T0_sec, T1_sec)
            
            A_phex_el = Q_ref / (U_el * dT_m_i)
            A_phex_elements[i] = A_phex_el

        dppl_ref = sum(dppl_ref_elements)
        dppl_sec = sum(dppl_sec_elements)
        A_phex = sum(A_phex_elements)
        dT_min = min(dT_elements)

        return A_phex, dppl_ref, dppl_sec, dT_min

    def _calc_dppt(self):
        # Refrigeration cycle
        self.fluid_ref.update(CP.HmassP_INPUTS, self.h_ref_in, self.p_ref)
        rho_ref_in = self.fluid_ref.rhomass() #CP.PropsSI("D", "H", self.h_ref_in, "P", self.p_ref, self.ref_fluid) # Density dependend on inlet temperature here
        upt_ref_in = 4 * self.m_flow_ref / (rho_ref_in * np.pi * self.Dp**2)
        dppt_ref_in = 1.3 * 0.5 * rho_ref_in * upt_ref_in**2

        self.fluid_ref.update(CP.HmassP_INPUTS, self.h_ref_out, self.p_ref)
        rho_ref_out = self.fluid_ref.rhomass()#CP.PropsSI("D", "H", self.h_ref_out, "P", self.p_ref, self.ref_fluid)
        upt_ref_out = 4 * self.m_flow_ref / (rho_ref_out * np.pi * self.Dp**2)
        dppt_ref_out = 1.3 * 0.5 * rho_ref_out * upt_ref_out**2

        dppt_ref = dppt_ref_in + dppt_ref_out

        # Secondary cycle
        self.fluid_sec.update(CP.HmassP_INPUTS, self.h_sec_in, self.p_sec)
        rho_sec_in = self.fluid_sec.rhomass()#CP.PropsSI("D", "H", self.h_sec_in, "P", self.p_sec, self.sec_fluid) # Density dependend on inlet temperature here
        upt_sec_in = 4 * self.m_flow_sec / (rho_sec_in * np.pi * self.Dp**2)
        dppt_sec_in = 1.3 * 0.5 * rho_sec_in * upt_sec_in**2 

        self.fluid_sec.update(CP.HmassP_INPUTS, self.h_sec_out, self.p_sec)
        rho_sec_out = self.fluid_sec.rhomass()#CP.PropsSI("D", "H", self.h_sec_out, "P", self.p_sec, self.sec_fluid)
        upt_sec_out = 4 * self.m_flow_sec / (rho_sec_out * np.pi * self.Dp**2)
        dppt_sec_out = 1.3 * 0.5 * rho_sec_out * upt_sec_out**2 

        dppt_sec = dppt_sec_in + dppt_sec_out
        return dppt_ref, dppt_sec

    def _calc_single_phase_alpha(self, T, Tw, p, m_flow, fluid, dh, tot_dh, h=None):
        # States are calculated using h, p instead of T, p if h != None, required for phase boundaries
        if h == None:
            fluid.update(CP.PT_INPUTS, p, T)
            #CP.PropsSI("V", "T", T, "P", p, fluid)
            #CP.PropsSI("C", "T", T, "P", p, fluid)
            #CP.PropsSI("L", "T", T, "P", p, fluid) # here lambda
            #CP.PropsSI("D", "T", T, "P", p, fluid)
        else:
            fluid.update(CP.HmassP_INPUTS, h, p)
            #mu = CP.PropsSI("V", "H", h, "P", p, fluid)
            #cp = CP.PropsSI("C", "H", h - 1e4, "P", p, fluid)
            #l = CP.PropsSI("L", "H", h, "P", p, fluid) # here lambda
            #rho = CP.PropsSI("D", "H", h, "P", p, fluid)

        mu = fluid.viscosity()
        cp = fluid.cpmass()
        l = fluid.conductivity()
        rho = fluid.rhomass()

        fluid.update(CP.PT_INPUTS, p, Tw)
        muw = fluid.viscosity() #CP.PropsSI("V", "T", Tw, "P", p, fluid)

        G_ch = m_flow / (self.N_cp * self.A_ch)

        Re = G_ch * self.D_h / mu
        Pr = cp * mu / l

        if Re < 2000:
            xi_0 = 64 / Re
            xi_1 = 597 / Re + 3.85
        else:
            xi_0 = (1.8*np.log(Re) - 1.5)**(-2)
            xi_1 = 39 / (Re)**0.289
        # 1 / sqrt(xi)
        inter_xi = np.cos(self.phi) / np.sqrt(0.18 * np.tan(self.phi) + 0.36 * np.sin(self.phi) + xi_0 / np.cos(self.phi)) + (1 - np.cos(self.phi)) / np.sqrt(3.8 * xi_1)
        xi = 1 / (inter_xi * inter_xi)

        Nu_c = 0.122 * Pr**(1/3) * (mu / muw)**(1/6) * (xi * Re*Re * np.sin(2 * self.phi))**0.374
        a = Nu_c * l / self.D_h

        if Re < 1600:
            f = 26.8 * Re**(-0.209)
        else:
            f = 26.8 * 1600**(-0.209)
        # scale by element contribution dh/dh_tot, not physicly correct but I don't have a better idea here
        # dh/dh_tot instead of 1/num_elements because elements may not have the same size
        dppl = f * dh / tot_dh * self.Lp / self.D_h * (G_ch**2 / (2 * rho))
        
        return a, dppl

    def _calc_condensation_alpha(self, T, p, h, m_flow, fluid, dh, tot_dh):
        fluid.update(CP.QT_INPUTS, 0, T)
        mu_L = fluid.viscosity() #CP.PropsSI("V", "T", T, "Q", 0, fluid)
        cp_L = fluid.cpmass() #CP.PropsSI("C", "T", T, "Q", 0, fluid)
        lambda_L = fluid.conductivity() #CP.PropsSI("L", "T", T, "Q", 0, fluid)
        h_l = fluid.hmass() #CP.PropsSI("H", "T", T, "Q", 0, fluid)
        rho_l = fluid.rhomass() #CP.PropsSI("D", "T", T, "Q", 0, fluid)
        fluid.update(CP.QT_INPUTS, 1, T)
        h_g = fluid.hmass()#CP.PropsSI("H", "T", T, "Q", 1, fluid)
        rho_g = fluid.rhomass() #CP.PropsSI("D", "T", T, "Q", 1, fluid)
        p_crit = fluid.p_critical() #CP.PropsSI("Pcrit", "", 0, "", 0, fluid)
        fluid.update(CP.HmassP_INPUTS, h, p)
        x = fluid.Q() #CP.PropsSI("Q", "H", h, "P", p, fluid)
        mu = fluid.viscosity() #CP.PropsSI("V", "H", h, "P", p, fluid)
        rho = fluid.rhomass() #CP.PropsSI("D", "H", h, "P", p, fluid)

        p_red = p / p_crit
        G_ch = m_flow / (self.N_cp * self.A_ch)

        Re = G_ch * self.D_h / mu
        Re_L = G_ch * self.D_h / mu_L
        Pr_L = cp_L * mu_L / lambda_L
        
        Nu_h = 0.023 * Re_L**0.8 * Pr_L**0.4*((1 - x)**0.8 + (3.8 * x**0.76 * (1 - x)**0.04) / (p_red**0.38)) # CHANGED one sign from (1-x)**0.8 * (3.8 ...) to (1-x)**0.8 + (3.8 ...)
    
        a = Nu_h * lambda_L / self.D_h

        Re_eq = G_ch * ((1 - x) + x * (rho_l / rho_g)**0.5)*self.D_h / mu_L

        q = m_flow * tot_dh / (self.Nt * self.A_p) #dh * m_flow / self.A_p # No idea how I should calculate the q
        B0 = q / (G_ch*(h_g - h_l))
        f = 94.75 * Re_eq**(-0.0467)*Re**(-0.4)*B0**(0.5)*p_red**0.8
        dppl = f * dh / tot_dh * self.Lp / self.D_h * (G_ch*G_ch / (2 * rho))
        return a, dppl

    def _calc_evaporation_alpha(self, T, p, h, m_flow, fluid, dh, tot_dh):
        fluid.update(CP.QT_INPUTS, 0, T)
        mu_L = fluid.viscosity() #CP.PropsSI("V", "T", T, "Q", 0, fluid)
        cp_L = fluid.cpmass() #CP.PropsSI("C", "T", T, "Q", 0, fluid)
        lambda_L = fluid.conductivity() #CP.PropsSI("L", "T", T, "Q", 0, fluid)
        h_l = fluid.hmass() #CP.PropsSI("H", "T", T, "Q", 0, fluid)
        rho_l = fluid.rhomass() #CP.PropsSI("D", "T", T, "Q", 0, fluid)
        fluid.update(CP.QT_INPUTS, 1, T)
        h_g = fluid.hmass()#CP.PropsSI("H", "T", T, "Q", 1, fluid)
        rho_g = fluid.rhomass() #CP.PropsSI("D", "T", T, "Q", 1, fluid)
        p_crit = fluid.p_critical() #CP.PropsSI("Pcrit", "", 0, "", 0, fluid)
        fluid.update(CP.HmassP_INPUTS, h, p)
        x = fluid.Q() #CP.PropsSI("Q", "H", h, "P", p, fluid)
        mu = fluid.viscosity() #CP.PropsSI("V", "H", h, "P", p, fluid)
        rho = fluid.rhomass() #CP.PropsSI("D", "H", h, "P", p, fluid)

        p_red = p / p_crit
        G_ch = m_flow / (self.N_cp * self.A_ch)
        G_ch_eq = G_ch * (1- x + x * (rho_l / rho_g)**0.5)

        Re_eq = G_ch_eq * self.D_h / mu_L
        Re_L = G_ch * self.D_h / mu_L
        Pr_L = cp_L * mu_L / lambda_L

        q = m_flow * tot_dh / (self.N_t0 * self.A_p)#dh * m_flow / self.A_p # No idea how I should calculate the Q
        B0_eq = q / (G_ch_eq*(h_g - h_l))

        Nu_h = 19.26 * Re_L**0.5 * B0_eq**0.3 * Pr_L ** (1/3)
        a = Nu_h * lambda_L / self.D_h

        if Re_eq < 6000:
            f = 6.947*10**5 * Re_L**(-0.5) * Re_eq**(-1.109)
        else:
            f = 31.21 * Re_L**(-0.5) * Re_eq**(0.04557)

        dppl = f * dh / tot_dh * self.Lp / self.D_h * (G_ch*G_ch / (2 * rho))
        return a, dppl
    

    def _calc_dT_m(self, T_ref_in, T_ref_out, T_sec_in, T_sec_out):
        # Assuming counter flow heat exchanger
        dT_A = T_ref_in - T_sec_out
        dT_B = T_ref_out - T_sec_in
        dT_lm = (dT_A - dT_B) / np.log(dT_A / dT_B)
        dT_m = abs(self.Ft * dT_lm)

        return dT_m
            
if __name__ == "__main__":
    # Example usage of the PlateHeatExchanger class
    geom = {
        "number_of_passes": 1,
        "plate_thickness_m": 0.0007,
        "chevron_angle_rad": 1.0471975511965976,
        "pitch": 0.0025,
        "plate_amplitude_m": 0.001,
        "corrugation_pitch_m": 0.007,
        "Dp": 0.023,
        "Lp": 0.25,
        "Bp": 0.113,
        "Ntmin": 4,
        "Ntmax": 150,
        "m_max": 14.0,
        "Nt": 100
    }

    fluid_ref = AbstractState("HEOS", "R134a")
    fluid_sec = AbstractState("HEOS", "Water")

    heat_exchanger = PlateHeatExchanger(geom)
    sol = heat_exchanger.calc(fluid_ref=fluid_ref,
                        fluid_sec=fluid_sec,
                        m_flow_ref=0.01,
                        m_flow_sec=0.1,
                        h_ref_in=433840,
                        h_ref_out=241720,
                        p_ref=770200,
                        p_sec=100000,
                        T_sec_in=20 + 273.15)
    print(sol)