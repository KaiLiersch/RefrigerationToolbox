import CoolProp.CoolProp as CP
import numpy as np
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.HeatExchanger import HeatExchanger


class PlateHeatExchanger(HeatExchanger):
    """Thermal model and sizing routine for a chevron (herringbone) plate heat exchanger.

    The exchanger transfers heat between a refrigerant stream (``fluid_ref``, which may
    change phase) and a single-phase secondary stream (``fluid_sec``). It can operate as
    a condenser or as an evaporator; which one is decided from the refrigerant enthalpy
    change. Adjacent plates form narrow channels through which the two streams flow in
    counter-current arrangement, and the corrugations enlarge the surface area and
    intensify heat transfer.

    Sizing follows a discretized approach: each plate is split into a series of finite
    elements along the flow direction. Within every element the local heat-transfer
    coefficient of each stream is evaluated with a correlation appropriate to its phase
    (single phase, condensation or evaporation), combined into a local overall
    coefficient, and turned into a required heat-transfer area. Summing the element areas
    gives the total area needed for the duty, which is compared against the area that the
    chosen number of plates actually provides. The frictional and port pressure drops are
    accumulated in the same loop so the design can also be checked against pressure limits.
    """

    def __init__(
        self,
        fluid_ref: AbstractState,
        fluid_sec: AbstractState,
        geom: dict,
        Rf_ref: float = 0.0001,
        Rf_sec: float = 0.00003,
        num_elements: int = 10,
        Ft: float = 1.0,
        l_w: float = 20,
        dp_ref_max: float = 2e4,
        dp_sec_max: float = 2e4,
        dT_pinch: float = 2,
        feas_tol: float = 1e-6,
    ) -> None:
        """Set up the exchanger from operating parameters and a plate geometry.

        Args:
            fluid_ref: CoolProp state object for the refrigerant (primary) stream.
            fluid_sec: CoolProp state object for the secondary (single-phase) stream.
            geom: Plate geometry dictionary. Expected keys are described below where the
                corresponding attributes are assigned.
            Rf_ref: Fouling (deposit) thermal resistance on the refrigerant side [m²K/W].
            Rf_sec: Fouling thermal resistance on the secondary side [m²K/W].
            num_elements: Number of finite elements used to discretize the plate.
            Ft: Correction factor applied to the logarithmic mean temperature difference.
            l_w: Thermal conductivity of the plate material [W/mK].
            dp_ref_max: Maximum allowed pressure drop on the refrigerant side [Pa].
            dp_sec_max: Maximum allowed pressure drop on the secondary side [Pa].
            dT_pinch: Minimum allowed temperature difference between the streams (pinch) [K].
            feas_tol: Relative tolerance by which each limit is relaxed in the feasibility
                check so that a design sitting exactly on a limit is not rejected on rounding.
        """
        super().__init__(fluid_ref, fluid_sec)
        self.num_elements = num_elements
        self.Ft = Ft
        self.l_w = l_w
        self.Rf_ref = Rf_ref
        self.Rf_sec = Rf_sec
        self.dp_ref_max = dp_ref_max
        self.dp_sec_max = dp_sec_max
        self.dT_pinch = dT_pinch
        # Relative tolerance by which each limit is relaxed in the feasibility check
        self.feas_tol = feas_tol

        self.Ntmin = geom["Ntmin"]  # Smallest plate count considered when searching for a design
        self.Ntmax = geom["Ntmax"]  # Largest plate count considered when searching for a design
        # Corrugation inclination angle measured against the main flow direction [rad]. The stored
        # value is half of the full chevron angle passed in through the geometry dictionary.
        self.phi = geom["chevron_angle_rad"] / 2
        # Corrugation amplitude [m], i.e. half the gap between two adjacent plates across the stack.
        self.a = geom["plate_amplitude_m"]
        self.L = geom["corrugation_pitch_m"]  # Corrugation pitch (wavelength of the pattern) [m]
        self.t = geom["plate_thickness_m"]  # Plate wall thickness [m]
        self.Nt = geom["Nt"]  # Number of plates in the pack
        # Dimensionless wave number of the corrugation, used to size the enlargement factor.
        self.X = 2 * np.pi * self.a / self.L
        # Surface enlargement factor: ratio of the corrugated area to the flat projected area.
        self.Phi = 1 / 6 * (1 + np.sqrt(1 + self.X**2) + 4 * np.sqrt(1 + 0.5 * self.X**2))

        self.Dp = geom["Dp"]  # Port (inlet/outlet nozzle) diameter [m]
        self.Bp = geom["Bp"]  # Width of the heat-transfer surface of a plate [m]
        self.Lp = geom["Lp"]  # Height (length) of the heat-transfer surface of a plate [m]
        self.A_0 = self.Lp * self.Bp  # Flat projected area of one plate (ignoring corrugation) [m²]
        self.A_p = self.Phi * self.A_0  # Effective heat-transfer area of one plate incl. corrugation [m²]
        self.D_h = 4 * self.a / self.Phi  # Hydraulic diameter of a channel between two plates [m]
        self.A_ch = 2 * self.a * self.Bp  # Free flow cross-section of one channel between two plates [m²]
        self.A_available = self.Nt * self.A_p  # Total heat-transfer area provided by the plate pack [m²]

        self.feas = False  # Whether the current design satisfies the area, pressure and pinch limits

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
        """Rate the exchanger for a fixed plate count and record whether the design is feasible.

        Given the operating point (mass flows, refrigerant inlet/outlet enthalpies, both
        pressures and the secondary inlet temperature) this computes the duty, the mean
        temperature difference, discretizes the plate, solves every element for its local
        heat-transfer coefficient and required area, and evaluates the pressure drops. The
        design is feasible when the required area fits within the installed area and both
        pressure drops and the pinch stay within their limits. Results are stored on the
        instance (``A_phex``, ``dp_ref``, ``dp_sec``, ``dT_min``, ``feas``).

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
        self.feas = False
        # Heat duty transferred between the streams [W]
        self.Q = m_flow_ref * abs(h_ref_in - h_ref_out)
        # Resolve temperature, entropy and density at both inlets and outlets of both streams
        self._determine_boundary_states()

        # dT_m = self._calc_dT_m(self.T_ref_in, self.T_ref_out, self.T_sec_in, self.T_sec_out)

        # Number of channels carrying each stream. The plate pack alternates streams, so one
        # stream sees roughly half of the gaps between the Nt plates. Clamped to at least one.
        self.N_cp = (self.Nt - 1) / 2
        if self.N_cp < 1:
            self.N_cp = 1

        self.A_available = self.Nt * self.A_p

        # Generate discretization with non uniform elements
        # Elements conform to boundaries of phase regions. E. g. an element can not contain both a superheated gas and two phase region
        if self.state == "cond":
            self._discretize_condenser()
        elif self.state == "evap":
            self._discretize_evaporator()

        # Solve on discretization
        self._calc_fe(self.phase_ref_el, self.h_ref_nodes, self.h_sec_nodes)

        # Calculate pressure drop in heat exchanger
        self._calc_dppt()

        self.dp_ref = self.dppt_ref + self.dppl_ref
        self.dp_sec = self.dppt_sec + self.dppl_sec

        # Check if heat exchanger is sufficiently sized
        # Every limit is relaxed by the relative tolerance feas_tol, otherwise a design that sits
        # exactly on one of the limits, for example the result of an optimization, is rejected
        # because of rounding errors in the last digits
        if self.A_phex <= self.A_available * (1 + self.feas_tol):
            if self.dp_ref <= self.dp_ref_max * (1 + self.feas_tol) and self.dp_sec <= self.dp_sec_max * (1 + self.feas_tol) and self.dT_min >= self.dT_pinch * (1 - self.feas_tol):
                # TODO: Add warnings
                self.feas = True

    def determine_min_n_plates(
        self,
        U0: float,
        m_flow_ref: float,
        m_flow_sec: float,
        h_ref_in: float,
        h_ref_out: float,
        p_ref: float,
        p_sec: float,
        T_sec_in: float,
    ) -> None:
        """Find the smallest plate count that yields a feasible design for the given duty.

        A starting plate count is estimated from an assumed overall heat-transfer
        coefficient ``U0``, then plates are added one at a time. For each trial count the
        exchanger is fully solved (discretization, heat transfer, pressure drop) and the
        area, pressure-drop and pinch limits are checked. The search stops at the first
        count that satisfies all limits, leaving ``self.Nt`` and ``self.feas`` set to that
        result, or at ``Ntmax`` if none is found.

        Args:
            U0: Assumed overall heat-transfer coefficient for the first area estimate [W/m²K].
            m_flow_ref: Refrigerant mass flow rate [kg/s].
            m_flow_sec: Secondary fluid mass flow rate [kg/s].
            h_ref_in: Refrigerant specific enthalpy at the inlet [J/kg].
            h_ref_out: Refrigerant specific enthalpy at the outlet [J/kg].
            p_ref: Refrigerant pressure [Pa].
            p_sec: Secondary fluid pressure [Pa].
            T_sec_in: Secondary fluid inlet temperature [K].
        """
        self.m_flow_ref = m_flow_ref
        self.m_flow_sec = m_flow_sec
        self.h_ref_in = h_ref_in
        self.h_ref_out = h_ref_out
        self.p_ref = p_ref
        self.p_sec = p_sec
        self.T_sec_in = T_sec_in

        self.state = "cond" if h_ref_in > h_ref_out else "evap"

        # Heat duty transferred between the streams [W]
        self.Q = m_flow_ref * abs(h_ref_in - h_ref_out)
        # Resolve temperature, entropy and density at both inlets and outlets of both streams
        self._determine_boundary_states()

        dT_m = self._calc_dT_m(self.T_ref_in, self.T_ref_out, self.T_sec_in, self.T_sec_out)

        # First area estimate from the assumed overall coefficient, and the plate count it implies
        self.A_phex_0 = self.Q / (U0 * dT_m)
        Nt0 = max(1, int(self.A_phex_0 / self.A_p))

        for added_plates in range(0, self.Ntmax - Nt0 + 1):
            self.Nt = Nt0 + added_plates

            self.N_cp = (self.Nt - 1) / 2
            if self.N_cp < 1:
                self.N_cp = 1

            A_available = self.Nt * self.A_p

            # Generate discretization with non uniform elements
            # Elements conform to boundaries of phase regions. E. g. an element can not contain both a superheated gas and two phase region
            if self.state == "cond":
                self._discretize_condenser()
            elif self.state == "evap":
                self._discretize_evaporator()

            # Solve on discretization
            self._calc_fe(self.phase_ref_el, self.h_ref_nodes, self.h_sec_nodes)

            # Calculate pressure drop in heat exchanger
            self._calc_dppt()

            self.dp_ref = self.dppt_ref + self.dppl_ref
            self.dp_sec = self.dppt_sec + self.dppl_sec

            # print(self.Nt, self.dp_sec, self.dp_ref, self.dT_min, self.A_phex, A_available)

            # Check if heat exchanger is sufficiently sized
            if self.A_phex <= A_available:
                # Check if pressure drops are sufficiently small
                if self.dp_ref <= self.dp_ref_max and self.dp_sec <= self.dp_sec_max:
                    if self.dT_min >= self.dT_pinch:
                        self.feas = True
                        break

    def _determine_boundary_states(self) -> None:
        """Fill in the temperature, entropy and density at the four stream boundaries.

        The refrigerant inlet and outlet states follow from its known enthalpies and
        pressure. The secondary inlet state follows from its temperature and pressure; its
        outlet enthalpy is then obtained from an energy balance on the duty ``Q`` and the
        outlet state resolved from that enthalpy.
        """
        self.fluid_ref.update(CP.HmassP_INPUTS, self.h_ref_in, self.p_ref)
        self.T_ref_in = self.fluid_ref.T()
        self.s_ref_in = self.fluid_ref.smass()
        self.d_ref_in = self.fluid_ref.rhomass()
        self.fluid_ref.update(CP.HmassP_INPUTS, self.h_ref_out, self.p_ref)
        self.T_ref_out = self.fluid_ref.T()
        self.s_ref_out = self.fluid_ref.smass()
        self.d_ref_out = self.fluid_ref.rhomass()
        self.fluid_sec.update(CP.PT_INPUTS, self.p_sec, self.T_sec_in)
        self.h_sec_in = self.fluid_sec.hmass()
        self.s_sec_in = self.fluid_sec.smass()
        self.d_sec_in = self.fluid_sec.rhomass()
        self.h_sec_out = self.h_sec_in + self.Q / self.m_flow_sec
        self.fluid_sec.update(CP.HmassP_INPUTS, self.h_sec_out, self.p_sec)
        self.T_sec_out = self.fluid_sec.T()
        self.s_sec_out = self.fluid_sec.smass()
        self.d_sec_out = self.fluid_sec.rhomass()

    def _discretize_condenser(self) -> None:
        """Build the finite-element grid for a condensing refrigerant stream.

        The refrigerant enters as (possibly superheated) vapour and leaves as (possibly
        subcooled) liquid, passing through the two-phase region in between. Element
        boundaries are aligned with the saturation points so that no element straddles two
        phase regions, which keeps the correlations valid element by element. Elements are
        distributed between the three regions in proportion to their share of the enthalpy
        change, and the node enthalpies of both streams are laid out from the refrigerant
        inlet (hottest end) to the outlet. Results are stored in ``phase_ref_el`` (0 = gas,
        1 = two-phase, 2 = liquid) and the node enthalpy arrays ``h_ref_nodes``/``h_sec_nodes``.
        """
        # Saturated liquid and saturated vapour enthalpies at the refrigerant pressure
        self.fluid_ref.update(CP.PQ_INPUTS, self.p_ref, 0.0)
        h_sat0_ref = self.fluid_ref.hmass()
        self.fluid_ref.update(CP.PQ_INPUTS, self.p_ref, 1.0)
        h_sat1_ref = self.fluid_ref.hmass()

        phase_ref_el = np.zeros(self.num_elements)  # Phase label per element (0 gas, 1 two-phase, 2 liquid)
        h_ref_nodes = np.zeros(self.num_elements + 1)  # Refrigerant enthalpy at each element boundary [J/kg]
        h_sec_nodes = np.zeros(self.num_elements + 1)  # Secondary enthalpy at each element boundary [J/kg]

        # Discretization placing nodes on phase boundaries to to increase accuracy
        tot_diff_h_ref = self.h_ref_in - self.h_ref_out  # Total refrigerant enthalpy change [J/kg]
        tot_diff_h_sec = self.h_sec_out - self.h_sec_in  # Total secondary enthalpy change [J/kg]
        # Share of the total refrigerant enthalpy change spent in each phase region.
        # Without superheat at the inlet or without subcooling at the outlet the corresponding phase
        # region does not exist and its fraction is zero
        enthalpy_frac_gas = max(0.0, (self.h_ref_in - h_sat1_ref) / tot_diff_h_ref)
        enthalpy_frac_liquid = max(0.0, (h_sat0_ref - self.h_ref_out) / tot_diff_h_ref)
        # enthaply_frac_vle = 1 - enthalpy_frac_gas - enthalpy_frac_liquid

        # Calculate the number of elements for each phase region
        num_el_gas = int(np.ceil(enthalpy_frac_gas * self.num_elements))
        num_el_liquid = int(np.ceil(enthalpy_frac_liquid * self.num_elements))
        num_el_vle = max(0, int(self.num_elements - num_el_gas - num_el_liquid))
        # Save which phase region each element has
        phase_ref_el[0:num_el_gas] = 0
        phase_ref_el[num_el_gas : num_el_gas + num_el_vle] = 1
        phase_ref_el[self.num_elements - num_el_liquid :] = 2

        # Calculate enthalpy step sizes for each phase region for refrigeration fluid
        # A phase region that does not exist gets no elements, so it also gets no step size instead
        # of a division by an element count of zero
        dh_ref_gas = 0 if num_el_gas == 0 else (self.h_ref_in - h_sat1_ref) / num_el_gas
        dh_ref_vle = 0 if num_el_vle == 0 else (h_sat1_ref - h_sat0_ref) / num_el_vle
        dh_ref_liquid = 0 if num_el_liquid == 0 else (h_sat0_ref - self.h_ref_out) / num_el_liquid

        # Matching secondary-fluid enthalpy step per region, scaled so both streams cover the
        # same fraction of their respective total enthalpy change across each element
        dh_sec_1 = dh_ref_gas / tot_diff_h_ref * tot_diff_h_sec  # step opposite the gas region
        dh_sec_2 = dh_ref_vle / tot_diff_h_ref * tot_diff_h_sec  # step opposite the two-phase region
        dh_sec_3 = dh_ref_liquid / tot_diff_h_ref * tot_diff_h_sec  # step opposite the liquid region

        # Calculate the enthalpies at each node depending on the phase region
        h_ref_nodes[0] = self.h_ref_in
        h_sec_nodes[0] = self.h_sec_out
        # The first few elements always contain superheated gas
        for i in range(1, num_el_gas + 1):
            h_ref_nodes[i] = h_ref_nodes[i - 1] - dh_ref_gas
            h_sec_nodes[i] = h_sec_nodes[i - 1] - dh_sec_1
        # Then two phase region elements
        for i in range(num_el_gas + 1, num_el_gas + num_el_vle + 1):
            h_ref_nodes[i] = h_ref_nodes[i - 1] - dh_ref_vle
            h_sec_nodes[i] = h_sec_nodes[i - 1] - dh_sec_2
        # Finally, subcooled liquid elements
        for i in range(num_el_gas + num_el_vle + 1, self.num_elements + 1):
            h_ref_nodes[i] = h_ref_nodes[i - 1] - dh_ref_liquid
            h_sec_nodes[i] = h_sec_nodes[i - 1] - dh_sec_3

        self.phase_ref_el = phase_ref_el
        self.h_ref_nodes = h_ref_nodes
        self.h_sec_nodes = h_sec_nodes

    def _discretize_evaporator(self) -> None:
        """Build the finite-element grid for an evaporating refrigerant stream.

        Mirror of :meth:`_discretize_condenser` for the reverse process: the refrigerant
        enters as (possibly subcooled) liquid, boils through the two-phase region and
        leaves as (possibly superheated) vapour, with its enthalpy increasing along the
        flow. Element boundaries again sit on the saturation points, elements are shared
        between the liquid, two-phase and gas regions, and the node enthalpies are laid out
        from the refrigerant inlet to the outlet. Results are stored in the same
        ``phase_ref_el`` / ``h_ref_nodes`` / ``h_sec_nodes`` attributes.
        """
        # Saturated liquid and saturated vapour enthalpies at the refrigerant pressure
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
        if enthalpy_frac_gas < 0:
            enthalpy_frac_gas = 0
        if enthalpy_frac_liquid < 0:
            enthalpy_frac_liquid = 0
        # enthaply_frac_vle = 1 - enthalpy_frac_gas - enthalpy_frac_liquid

        num_el_gas = int(np.ceil(enthalpy_frac_gas * self.num_elements))
        num_el_liquid = int(np.ceil(enthalpy_frac_liquid * self.num_elements))
        num_el_vle = int(self.num_elements - num_el_gas - num_el_liquid)
        phase_ref_el[0:num_el_liquid] = 2
        phase_ref_el[num_el_liquid : num_el_liquid + num_el_vle] = 1
        phase_ref_el[self.num_elements - num_el_gas :] = 0

        dh_ref_gas = 0 if num_el_gas == 0 else (self.h_ref_out - h_sat1_ref) / num_el_gas
        dh_ref_vle = 0 if num_el_vle == 0 else (h_sat1_ref - max(h_sat0_ref, self.h_ref_in)) / num_el_vle
        dh_ref_liquid = 0 if num_el_liquid == 0 else (h_sat0_ref - self.h_ref_in) / num_el_liquid

        dh_sec_1 = dh_ref_gas / tot_diff_h_ref * tot_diff_h_sec
        dh_sec_2 = dh_ref_vle / tot_diff_h_ref * tot_diff_h_sec
        dh_sec_3 = dh_ref_liquid / tot_diff_h_ref * tot_diff_h_sec

        h_ref_nodes[0] = self.h_ref_in
        h_sec_nodes[0] = self.h_sec_out
        for i in range(1, num_el_liquid + 1):
            h_ref_nodes[i] = h_ref_nodes[i - 1] + dh_ref_liquid
            h_sec_nodes[i] = h_sec_nodes[i - 1] + dh_sec_3
        for i in range(num_el_liquid + 1, num_el_liquid + num_el_vle + 1):
            h_ref_nodes[i] = h_ref_nodes[i - 1] + dh_ref_vle
            h_sec_nodes[i] = h_sec_nodes[i - 1] + dh_sec_2
        for i in range(num_el_liquid + num_el_vle + 1, self.num_elements + 1):
            h_ref_nodes[i] = h_ref_nodes[i - 1] + dh_ref_gas
            h_sec_nodes[i] = h_sec_nodes[i - 1] + dh_sec_1

        self.phase_ref_el = phase_ref_el
        self.h_ref_nodes = h_ref_nodes
        self.h_sec_nodes = h_sec_nodes

    def _calc_fe(self, phase_ref_el: np.ndarray, h_ref_nodes: np.ndarray, h_sec_nodes: np.ndarray) -> None:
        """Solve the exchanger element by element and aggregate area, pressure drop and pinch.

        For every finite element this evaluates the stream temperatures (at both nodes and
        at the centre) and the wall temperatures, selects a heat-transfer correlation for
        each side according to its phase, forms the local overall coefficient (including
        both fouling resistances and the plate wall), and converts the local duty into the
        local area required. Frictional plate pressure drops are accumulated at the same
        time. On return the summed results are stored as ``A_phex`` (total required area),
        ``dppl_ref`` / ``dppl_sec`` (plate friction pressure drop per stream) and
        ``dT_min`` (smallest stream-to-stream temperature difference, i.e. the pinch).

        Args:
            phase_ref_el: Phase label of each element (0 gas, 1 two-phase, 2 liquid).
            h_ref_nodes: Refrigerant enthalpy at each element boundary [J/kg].
            h_sec_nodes: Secondary fluid enthalpy at each element boundary [J/kg].
        """
        # Per-element result buffers
        A_phex_elements = np.zeros(self.num_elements)  # Required heat-transfer area per element [m²]
        dppl_ref_elements = np.zeros(self.num_elements)  # Refrigerant plate friction pressure drop per element [Pa]
        dppl_sec_elements = np.zeros(self.num_elements)  # Secondary plate friction pressure drop per element [Pa]
        dT_elements = np.zeros(self.num_elements)  # Stream-to-stream temperature difference per element [K]

        # Total enthalpy difference between inlet and outlet of refrigeration fluid
        tot_dh = abs(h_ref_nodes[0] - h_ref_nodes[-1])

        a_ref_elements = np.zeros(self.num_elements)  # Refrigerant heat-transfer coefficient per element [W/m²K]
        a_sec_elements = np.zeros(self.num_elements)  # Secondary heat-transfer coefficient per element [W/m²K]

        # Loop over all elements
        for i, phase in enumerate(phase_ref_el):
            # Enthalpy of element is the average of enthalpy at nodes
            h_ref_el = (h_ref_nodes[i] + h_ref_nodes[i + 1]) / 2
            h_sec_el = (h_sec_nodes[i] + h_sec_nodes[i + 1]) / 2

            # Calculations are based on the enthalpies at the elements center
            self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_el, self.p_ref)
            T_ref_el = self.fluid_ref.T()  # CP.PropsSI("T", "H", h_ref_el, "P", self.p_ref, self.ref_fluid)
            self.fluid_sec.update(CP.HmassP_INPUTS, h_sec_el, self.p_sec)
            T_sec_el = self.fluid_sec.T()  # CP.PropsSI("T", "H", h_sec_el, "P", self.p_sec, self.sec_fluid)
            # Estimate the plate temperature as the mean of the two stream temperatures, and each
            # wall (film) temperature as the mean of its stream and the plate.
            T_plate_el = (T_ref_el + T_sec_el) / 2
            Tw_ref_el = (T_ref_el + T_plate_el) / 2  # Wall temperature seen by the refrigerant film [K]
            Tw_sec_el = (T_sec_el + T_plate_el) / 2  # Wall temperature seen by the secondary film [K]

            self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_nodes[i], self.p_ref)
            T0_ref = self.fluid_ref.T()
            self.fluid_ref.update(CP.HmassP_INPUTS, h_ref_nodes[i + 1], self.p_ref)
            T1_ref = self.fluid_ref.T()
            self.fluid_sec.update(CP.HmassP_INPUTS, h_sec_nodes[i], self.p_sec)
            T0_sec = self.fluid_sec.T()
            self.fluid_sec.update(CP.HmassP_INPUTS, h_sec_nodes[i + 1], self.p_sec)
            T1_sec = self.fluid_sec.T()

            # Calculation of temperature difference between hot and cold side
            # Takes minimum value of difference of both nodes and in the center of the element
            if self.state == "cond":
                dT_elements[i] = min([T_ref_el - T_sec_el, T0_ref - T0_sec, T1_ref - T1_sec])
            if self.state == "evap":
                dT_elements[i] = min([T_sec_el - T_ref_el, T0_sec - T0_ref, T1_sec - T1_ref])

            # Enthalpy change of the element
            dh_ref = abs(h_ref_nodes[i] - h_ref_nodes[i + 1])
            # dh_sec = abs(h_sec_nodes[i] - h_sec_nodes[i + 1])

            # Secondary fluid is always a single phase fluid
            a_sec, dppl_sec_el = self._calc_single_phase_alpha(T_sec_el, Tw_sec_el, self.p_sec, self.m_flow_sec, self.fluid_sec, dh_ref, tot_dh)

            # Refrigeration fluid might be two phase or single phase
            # 0 - gas; 1 - vle; 2 - liquid
            if phase == 0 or phase == 2:
                a_ref, dppl_ref_el = self._calc_single_phase_alpha(T_ref_el, Tw_ref_el, self.p_ref, self.m_flow_ref, self.fluid_ref, dh_ref, tot_dh, h_ref_el)
            elif phase == 1:
                if self.state == "cond":
                    a_ref, dppl_ref_el = self._calc_condensation_alpha(T_ref_el, self.p_ref, h_ref_el, self.m_flow_ref, self.fluid_ref, dh_ref, tot_dh)
                elif self.state == "evap":
                    # Not required for this exercise
                    a_ref, dppl_ref_el = self._calc_evaporation_alpha(T_ref_el, self.p_ref, h_ref_el, self.m_flow_ref, self.fluid_ref, dh_ref, tot_dh)
            a_ref_elements[i] = a_ref
            a_sec_elements[i] = a_sec

            dppl_ref_elements[i] = dppl_ref_el
            dppl_sec_elements[i] = dppl_sec_el

            # Local overall coefficient from the two film coefficients in series with the two
            # fouling resistances and the plate wall conduction [W/m²K]
            U_el = 1 / (1 / a_ref + 1 / a_sec + self.Rf_ref + self.Rf_sec + self.t / self.l_w)

            Q_ref = dh_ref * self.m_flow_ref  # Duty of this element from the refrigerant side [W]
            # Q_sec = dh_sec * self.m_flow_sec  # Duty of this element from the secondary side [W]

            dT_m_i = self._calc_dT_m(T0_ref, T1_ref, T0_sec, T1_sec)  # Mean temperature difference of the element [K]

            # Area needed to transfer this element's duty at its local coefficient and driving force
            A_phex_el = Q_ref / (U_el * dT_m_i)
            A_phex_elements[i] = A_phex_el

        self.dppl_ref = sum(dppl_ref_elements)
        self.dppl_sec = sum(dppl_sec_elements)
        self.A_phex = sum(A_phex_elements)
        self.dT_min = min(dT_elements)

    def _calc_dppt(self) -> None:
        """Compute the port (nozzle) pressure drop for both streams and store them.

        Sets ``dppt_ref`` and ``dppt_sec``, the pressure losses at the inlet/outlet ports
        that add to the frictional plate losses to give the total pressure drop per stream.
        """
        # Refrigerant stream
        self.dppt_ref = self._calc_dppt_stream(self.fluid_ref, self.m_flow_ref, self.h_ref_in, self.h_ref_out, self.p_ref)
        # Secondary stream
        self.dppt_sec = self._calc_dppt_stream(self.fluid_sec, self.m_flow_sec, self.h_sec_in, self.h_sec_out, self.p_sec)

    def _calc_dppt_stream(self, fluid: AbstractState, m_flow: float, h_in: float, h_out: float, p: float) -> float:
        """Port pressure drop of a single stream from its port mass velocity.

        The loss is taken as a fixed multiple of the port velocity head, using the inlet
        specific volume.

        Args:
            fluid: CoolProp state object for the stream.
            m_flow: Stream mass flow rate [kg/s].
            h_in: Inlet specific enthalpy [J/kg].
            h_out: Outlet specific enthalpy [J/kg] (currently unused; see note below).
            p: Stream pressure [Pa].

        Returns:
            Port pressure drop [Pa].
        """
        fluid.update(CP.HmassP_INPUTS, h_in, p)
        v_in = 1 / fluid.rhomass()  # Specific volume at the inlet [m³/kg]
        fluid.update(CP.HmassP_INPUTS, h_out, p)
        # Only the inlet specific volume is used here; the outlet/mean-volume variant is kept
        # commented out for reference.
        # v_out = 1 / fluid.rhomass()
        # v_mean = 0.5 * (v_in + v_out)

        # Mass flow rate per port cross section [kg/m²s]
        G_pt = 4 * m_flow / (np.pi * self.Dp**2)

        return 1.3 * 0.5 * G_pt**2 * v_in

    def _calc_single_phase_alpha(
        self,
        T: float,
        Tw: float,
        p: float,
        m_flow: float,
        fluid: AbstractState,
        dh: float,
        tot_dh: float,
        h: float | None = None,
    ) -> tuple[float, float]:
        """Single-phase heat-transfer coefficient and plate friction drop of one element.

        Uses a chevron-plate Nusselt correlation (with a wall-viscosity correction and a
        Darcy friction factor derived from the corrugation geometry) to get the film
        coefficient, and a separate single-phase friction correlation for the frictional
        pressure loss over the element's share of the plate length.

        Args:
            T: Bulk temperature of the stream in the element [K].
            Tw: Wall (film) temperature used for the viscosity correction [K].
            p: Stream pressure [Pa].
            m_flow: Stream mass flow rate [kg/s].
            fluid: CoolProp state object for the stream.
            dh: Enthalpy change across this element [J/kg].
            tot_dh: Total enthalpy change of the stream across the exchanger [J/kg].
            h: Optional specific enthalpy [J/kg]; when given the state is fixed from (h, p)
                instead of (p, T), which is needed right on a saturation boundary where the
                temperature alone is ambiguous.

        Returns:
            Tuple of the film heat-transfer coefficient [W/m²K] and the element's plate
            friction pressure drop [Pa].
        """
        # States are calculated using h, p instead of T, p if h != None, required for phase boundaries
        if h is None:
            fluid.update(CP.PT_INPUTS, p, T)
        else:
            fluid.update(CP.HmassP_INPUTS, h, p)

        mu = fluid.viscosity()  # Dynamic viscosity at bulk conditions [Pa·s]
        cp = fluid.cpmass()  # Specific heat capacity [J/kgK]
        l_therm = fluid.conductivity()  # Thermal conductivity [W/mK]
        rho = fluid.rhomass()  # Density [kg/m³]

        fluid.update(CP.PT_INPUTS, p, Tw)
        muw = fluid.viscosity()  # Dynamic viscosity at the wall temperature [Pa·s]

        G_ch = m_flow / (self.N_cp * self.A_ch)  # Mass velocity in one channel [kg/m²s]

        Re = G_ch * self.D_h / mu  # Reynolds number
        Pr = cp * mu / l_therm  # Prandtl number

        # Friction factor split into a smooth-channel term (xi_0) and a corrugation term (xi_1),
        # each with a laminar (Re < 2000) and a turbulent branch
        if Re < 2000:
            xi_0 = 64 / Re
            xi_1 = 597 / Re + 3.85
        else:
            xi_0 = (1.8 * np.log(Re) - 1.5) ** (-2)
            xi_1 = 39 / (Re) ** 0.289
        # inter_xi is 1 / sqrt(xi); the chevron angle blends the two contributions
        inter_xi = np.cos(self.phi) / np.sqrt(0.18 * np.tan(self.phi) + 0.36 * np.sin(self.phi) + xi_0 / np.cos(self.phi)) + (1 - np.cos(self.phi)) / np.sqrt(3.8 * xi_1)
        xi = 1 / (inter_xi * inter_xi)  # Darcy friction factor of the corrugated channel

        Nu_c = 0.122 * Pr ** (1 / 3) * (mu / muw) ** (1 / 6) * (xi * Re * Re * np.sin(2 * self.phi)) ** 0.374  # Nusselt number
        a = Nu_c * l_therm / self.D_h  # Film heat-transfer coefficient [W/m²K]

        # Fanning-type friction factor for the pressure-drop correlation, branched on Reynolds
        # number and clamped at the ends of its validity range
        if Re < 90:
            f = 5.03 + 755 / 90
        elif Re < 400:
            f = 5.03 + 755 / Re  # 26.8 * Re**(-0.209)
        elif Re < 16000:
            f = 26.8 * Re ** (-0.209)
        else:
            f = 26.8 * 16000 ** (-0.209)
        # Element friction drop over its portion of the plate. The dh/tot_dh weight (instead of
        # 1/num_elements) shares the plate length by enthalpy change, since elements differ in size.
        dppl = f * dh / tot_dh * self.Lp / self.D_h * (G_ch**2 / (2 * rho))

        return a, dppl

    def _calc_condensation_alpha(self, T: float, p: float, h: float, m_flow: float, fluid: AbstractState, dh: float, tot_dh: float) -> tuple[float, float]:
        """Condensation heat-transfer coefficient and plate friction drop of one element.

        Applies a two-phase condensation Nusselt correlation evaluated at saturated-liquid
        properties and modulated by the local vapour quality and reduced pressure, together
        with a two-phase friction correlation for the frictional pressure loss. Fluid
        properties are taken at the saturation states bracketing the element.

        Args:
            T: Saturation temperature of the condensing stream in the element [K].
            p: Stream pressure [Pa].
            h: Specific enthalpy at the element centre [J/kg], used to obtain the quality.
            m_flow: Stream mass flow rate [kg/s].
            fluid: CoolProp state object for the condensing stream.
            dh: Enthalpy change across this element [J/kg].
            tot_dh: Total enthalpy change of the stream across the exchanger [J/kg].

        Returns:
            Tuple of the film heat-transfer coefficient [W/m²K] and the element's plate
            friction pressure drop [Pa].
        """
        fluid.update(CP.QT_INPUTS, 0, T)
        mu_L = fluid.viscosity()  # Saturated-liquid viscosity [Pa·s]
        cp_L = fluid.cpmass()  # Saturated-liquid specific heat [J/kgK]
        lambda_L = fluid.conductivity()  # Saturated-liquid thermal conductivity [W/mK]
        h_l = fluid.hmass()  # Saturated-liquid enthalpy [J/kg]
        rho_l = fluid.rhomass()  # Saturated-liquid density [kg/m³]
        fluid.update(CP.QT_INPUTS, 1, T)
        h_g = fluid.hmass()  # Saturated-vapour enthalpy [J/kg]
        rho_g = fluid.rhomass()  # Saturated-vapour density [kg/m³]
        p_crit = fluid.p_critical()  # Critical pressure of the fluid [Pa]
        fluid.update(CP.HmassP_INPUTS, h, p)
        x = fluid.Q()  # Local vapour quality [-]
        mu = fluid.viscosity()  # Two-phase mixture viscosity [Pa·s]
        rho = fluid.rhomass()  # Two-phase mixture density [kg/m³]

        p_red = p / p_crit  # Reduced pressure [-]
        G_ch = m_flow / (self.N_cp * self.A_ch)  # Mass velocity in one channel [kg/m²s]

        Re = G_ch * self.D_h / mu  # Two-phase Reynolds number
        Re_L = G_ch * self.D_h / mu_L  # Reynolds number as if all mass flowed as saturated liquid
        Pr_L = cp_L * mu_L / lambda_L  # Saturated-liquid Prandtl number

        Nu_h = (
            0.023 * Re_L**0.8 * Pr_L**0.4 * ((1 - x) ** 0.8 + (3.8 * x**0.76 * (1 - x) ** 0.04) / (p_red**0.38))
        )  # CHANGED one sign from (1-x)**0.8 * (3.8 ...) to (1-x)**0.8 + (3.8 ...)

        a = Nu_h * lambda_L / self.D_h  # Film heat-transfer coefficient [W/m²K]

        Re_eq = G_ch * ((1 - x) + x * (rho_l / rho_g) ** 0.5) * self.D_h / mu_L  # Equivalent all-liquid Reynolds number

        q = m_flow * tot_dh / (self.Nt * self.A_p)  # Average heat flux over the plate area [W/m²]
        B0 = q / (G_ch * (h_g - h_l))  # Boiling number [-]

        # Two-phase friction factor, branched on Reynolds number with the ends of its range clamped
        if Re < 500:
            f = 94.75 * Re_eq ** (-0.0467) * 500 ** (-0.4) * B0 ** (0.5) * p_red**0.8
        elif Re < 10000:
            f = 94.75 * Re_eq ** (-0.0467) * Re ** (-0.4) * B0 ** (0.5) * p_red**0.8
        else:
            f = 94.75 * Re_eq ** (-0.0467) * 1e5 ** (-0.4) * B0 ** (0.5) * p_red**0.8

        # Element friction drop over its enthalpy-weighted share of the plate length
        dppl = f * dh / tot_dh * self.Lp / self.D_h * (G_ch * G_ch / (2 * rho))
        return a, dppl

    def _calc_evaporation_alpha(self, T: float, p: float, h: float, m_flow: float, fluid: AbstractState, dh: float, tot_dh: float) -> tuple[float, float]:
        """Evaporation heat-transfer coefficient and plate friction drop of one element.

        Applies a two-phase boiling Nusselt correlation based on the all-liquid Reynolds
        number, a boiling number and the saturated-liquid Prandtl number, together with a
        two-phase friction correlation using an equivalent (all-liquid-equivalent) mass
        velocity.

        Args:
            T: Saturation temperature of the evaporating stream in the element [K].
            p: Stream pressure [Pa].
            h: Specific enthalpy at the element centre [J/kg], used to obtain the quality.
            m_flow: Stream mass flow rate [kg/s].
            fluid: CoolProp state object for the evaporating stream.
            dh: Enthalpy change across this element [J/kg].
            tot_dh: Total enthalpy change of the stream across the exchanger [J/kg].

        Returns:
            Tuple of the film heat-transfer coefficient [W/m²K] and the element's plate
            friction pressure drop [Pa].
        """
        fluid.update(CP.QT_INPUTS, 0, T)
        mu_L = fluid.viscosity()  # Saturated-liquid viscosity [Pa·s]
        cp_L = fluid.cpmass()  # Saturated-liquid specific heat [J/kgK]
        lambda_L = fluid.conductivity()  # Saturated-liquid thermal conductivity [W/mK]
        h_l = fluid.hmass()  # Saturated-liquid enthalpy [J/kg]
        rho_l = fluid.rhomass()  # Saturated-liquid density [kg/m³]
        fluid.update(CP.QT_INPUTS, 1, T)
        h_g = fluid.hmass()  # Saturated-vapour enthalpy [J/kg]
        rho_g = fluid.rhomass()  # Saturated-vapour density [kg/m³]
        # p_crit = fluid.p_critical()  # Critical pressure of the fluid [Pa]
        fluid.update(CP.HmassP_INPUTS, h, p)
        x = fluid.Q()  # Local vapour quality [-]
        # mu = fluid.viscosity()  # Two-phase mixture viscosity [Pa·s]
        rho = fluid.rhomass()  # Two-phase mixture density [kg/m³]

        # p_red = p / p_crit  # Reduced pressure [-]
        G_ch = m_flow / (self.N_cp * self.A_ch)  # Mass velocity in one channel [kg/m²s]
        G_ch_eq = G_ch * (1 - x + x * (rho_l / rho_g) ** 0.5)  # Equivalent all-liquid mass velocity [kg/m²s]

        Re_eq = G_ch_eq * self.D_h / mu_L  # Equivalent all-liquid Reynolds number
        Re_L = G_ch * self.D_h / mu_L  # Reynolds number as if all mass flowed as saturated liquid
        Pr_L = cp_L * mu_L / lambda_L  # Saturated-liquid Prandtl number

        q = m_flow * tot_dh / (self.Nt * self.A_p)  # dh * m_flow / self.A_p # No idea how I should calculate the Q
        B0_eq = q / (G_ch_eq * (h_g - h_l))  # Equivalent boiling number [-]

        Nu_h = 19.26 * Re_L**0.5 * B0_eq**0.3 * Pr_L ** (1 / 3)  # Nusselt number
        a = Nu_h * lambda_L / self.D_h  # heat-transfer coefficient [W/m²K]

        # Two-phase friction factor, branched on the equivalent Reynolds number
        if Re_eq < 6000:
            f = 6.947 * 10**5 * Re_L ** (-0.5) * Re_eq ** (-1.109)
        else:
            f = 31.21 * Re_L ** (-0.5) * Re_eq ** (0.04557)

        # Element friction drop over its enthalpy-weighted share of the plate length
        dppl = f * dh / tot_dh * self.Lp / self.D_h * (G_ch * G_ch / (2 * rho))
        return a, dppl

    def _calc_dT_m(self, T_ref_in: float, T_ref_out: float, T_sec_in: float, T_sec_out: float) -> float:
        """Corrected logarithmic mean temperature difference for a counter-flow arrangement.

        Returns the logarithmic mean temperature difference between the two streams scaled
        by the correction factor ``Ft``. When both end differences are (near) equal the log
        mean is undefined, so the common end value is used directly to avoid a division by
        zero on zero-width elements.

        Args:
            T_ref_in: Refrigerant-side temperature at the inlet of the interval [K].
            T_ref_out: Refrigerant-side temperature at the outlet of the interval [K].
            T_sec_in: Secondary-side temperature at the inlet of the interval [K].
            T_sec_out: Secondary-side temperature at the outlet of the interval [K].

        Returns:
            Corrected mean temperature difference [K].
        """
        # Assuming counter flow heat exchanger
        dT_A = T_ref_in - T_sec_out  # Temperature difference at one end of the interval [K]
        dT_B = T_ref_out - T_sec_in  # Temperature difference at the other end of the interval [K]
        # Fall back to the plain end difference when both ends are equal (zero width elements)
        if abs(dT_A - dT_B) < 1e-9 * max(abs(dT_A), abs(dT_B), 1e-9):
            return abs(self.Ft * dT_A)
        dT_lm = (dT_A - dT_B) / np.log(dT_A / dT_B)
        dT_m = abs(self.Ft * dT_lm)

        return dT_m

    def __str__(self) -> str:
        """Human-readable summary extending the base description with plate-specific results.

        Adds the feasibility flag, plate count and the two stream pressure drops and pinch,
        which are only available after :meth:`calc` or :meth:`determine_min_n_plates` has run.
        """

        def fmt(val: float | None, scale: float = 1) -> str:
            return "N/A" if val is None else f"{val / scale:.2f}"

        dp_ref = getattr(self, "dp_ref", None)
        dp_sec = getattr(self, "dp_sec", None)
        return (
            f"{super().__str__()}\n Feasable = {self.feas}, Nt = {self.Nt}, dp_ref = {fmt(dp_ref, 1e3)} kPa, dp_sec = {fmt(dp_sec, 1e3)} kPa, self.dT_min = {fmt(self.dT_min, 1)} K"
        )
