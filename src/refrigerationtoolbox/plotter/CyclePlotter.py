import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.BasicHeatExchanger import BasicHeatExchanger
from refrigerationtoolbox.cycle.Compressor import celsius_to_kelvin_coeffs
from refrigerationtoolbox.cycle.Cycle import Cycle
from refrigerationtoolbox.cycle.EffCompressor import EffCompressor
from refrigerationtoolbox.cycle.PlateHeatExchanger import PlateHeatExchanger
from refrigerationtoolbox.cycle.PolynomialCompressor import PolynomialCompressor

RED = "#D81B60"
BLUE = "#1E88E5"
YELLOW = "#FFC107"
GREEN = "#004D40"


class CyclePlotter:
    # Default colours for the (dashed) isotherm / isobar families.
    _iso_colors = ["tab:green", "tab:orange", "tab:red"]

    def __init__(self):
        pass

    # --- log p-h diagram --- #

    def plot_ph(self, ax, cycle : Cycle):
        fluid = cycle.fluid

        p_min = cycle.pe * 0.4               # a little below evaporating pressure
        p_max = fluid.p_critical() * 1.05    # show the closed top of the dome

        self.plot_ph_dome(ax, fluid, color="black")
        self.plot_ph_cycle(ax, cycle)
        # Enthalpy extent from the dome and cycle, captured before the isotherm
        # tails widen it so they get clipped to a sensible range.
        h_lo, h_hi = ax.get_xlim()
        self.plot_ph_isotherms(ax, fluid, [cycle.Te, cycle.Tc, cycle.T2],
                               p_min=p_min, p_max=p_max, colors=[RED, RED, RED])

        ax.set_yscale("log")
        ax.set_xlabel("Specific enthalpy $h$ (kJ/kg)")
        ax.set_ylabel("Pressure $p$ (kPa)")
        ax.set_title(f"log p-h diagram ({fluid.fluid_names()[0]})")
        ax.set_xlim(h_lo, h_hi)
        ax.set_ylim(p_min / 1e3, p_max / 1e3)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(loc="best", fontsize=8)
        return ax

    def plot_ph_dome(self, ax, fluid : AbstractState, color="black", lw=1.5, label="Saturation dome",
                     linestyle="-"):
        """Plot the two-phase saturation dome on p-h axes (h in kJ/kg, p in kPa)."""
        Tcrit = fluid.T_critical()
        T_sat = np.linspace(fluid.Ttriple() + 0.1, Tcrit - 0.1, 300)
        h_liq, p_liq, h_vap, p_vap = [], [], [], []
        for T in T_sat:
            fluid.update(CP.QT_INPUTS, 0.0, T)   # saturated liquid
            h_liq.append(fluid.hmass())
            p_liq.append(fluid.p())
            fluid.update(CP.QT_INPUTS, 1.0, T)   # saturated vapour
            h_vap.append(fluid.hmass())
            p_vap.append(fluid.p())

        # One continuous curve: up the liquid line, back down the vapour line.
        h_dome = np.concatenate([h_liq, h_vap[::-1]])
        p_dome = np.concatenate([p_liq, p_vap[::-1]])
        ax.plot(h_dome / 1e3, p_dome / 1e3, color=color, lw=lw,
                label=label, zorder=2, linestyle=linestyle)
        return ax

    def plot_ph_isotherms(self, ax, fluid : AbstractState, temps,
                          p_min=None, p_max=None, colors=None):
        Tcrit = fluid.T_critical()
        pcrit = fluid.p_critical()
        p_min = pcrit * 0.02 if p_min is None else p_min
        p_max = pcrit * 1.05 if p_max is None else p_max
        colors = self._iso_colors if colors is None else colors
        for i, T in enumerate(temps):
            h_iso, p_iso = self._isotherm(fluid, T, p_min, p_max, Tcrit)
            ax.plot(h_iso / 1e3, p_iso / 1e3, ls="--", lw=1.0,
                    color=colors[i % len(colors)],
                    label=f"T = {T - 273.15:.0f} °C", zorder=3)
        return ax

    def plot_ph_cycle(self, ax, cycle : Cycle, color=BLUE, lw=2.0,
                       textcoords="offset points", xytext=(7, 7), fontweight="bold",
                       label="Cycle", markersize=6):
        h_states = np.array([cycle.h1, cycle.h2, cycle.h3, cycle.h4, cycle.h1])
        p_states = np.array([cycle.p1, cycle.p2, cycle.p3, cycle.p4, cycle.p1])
        ax.plot(h_states / 1e3, p_states / 1e3, "-o", color=color, lw=lw,
                markersize=markersize, label=label, zorder=5)
        for i, (h, p) in enumerate(zip(h_states[:4], p_states[:4], strict=True), start=1):
            ax.annotate(str(i), (h / 1e3, p / 1e3), textcoords=textcoords,
                        xytext=xytext, fontweight="bold", zorder=6)
        return ax

    # --- T-s diagram --- #

    def plot_ts(self, ax, cycle : Cycle):
        fluid = cycle.fluid

        T_min = cycle.Te - 40.0                # a little below evaporating temp.
        T_max = fluid.T_critical() + 15.0      # a little above the critical temp.

        self.plot_ts_dome(ax, fluid, T_min=T_min)
        self.plot_ts_cycle(ax, cycle)
        # Entropy extent from the dome and cycle, captured before the isobar
        # tails widen it so they get clipped to a sensible range.
        s_lo, s_hi = ax.get_xlim()
        p_high = min(cycle.pc * 2.0, fluid.p_critical() * 0.9)
        self.plot_ts_isobars(ax, fluid, [cycle.pe, cycle.pc, p_high],
                             T_min=T_min, T_max=T_max, colors=[RED, RED, RED])

        ax.set_xlabel("Specific entropy $s$ (kJ/kgK)")
        ax.set_ylabel("Temperature $T$ (K)")
        ax.set_title(f"T-s diagram ({fluid.fluid_names()[0]})")
        ax.set_xlim(s_lo, s_hi)
        ax.set_ylim(T_min, T_max)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(loc="best", fontsize=8)
        return ax

    def plot_ts_dome(self, ax, fluid : AbstractState, T_min=None, color="black", lw=1.5, label="Saturation dome"):
        Tcrit = fluid.T_critical()
        T_min = fluid.Ttriple() + 0.1 if T_min is None else T_min
        T_sat = np.linspace(T_min, Tcrit - 0.1, 300)
        s_liq, T_liq, s_vap, T_vap = [], [], [], []
        for T in T_sat:
            fluid.update(CP.QT_INPUTS, 0.0, T)   # saturated liquid
            s_liq.append(fluid.smass())
            T_liq.append(T)
            fluid.update(CP.QT_INPUTS, 1.0, T)   # saturated vapour
            s_vap.append(fluid.smass())
            T_vap.append(T)

        # One continuous curve: up the liquid line, back down the vapour line.
        s_dome = np.concatenate([s_liq, s_vap[::-1]])
        T_dome = np.concatenate([T_liq, T_vap[::-1]])
        ax.plot(s_dome / 1e3, T_dome, color=color, lw=lw,
                label=label, zorder=2)
        return ax

    def plot_ts_isobars(self, ax, fluid : AbstractState, pressures,
                        T_min=None, T_max=None, colors=None, ls="--", lw=1.0):
        Tcrit = fluid.T_critical()
        pcrit = fluid.p_critical()
        T_min = fluid.Ttriple() + 0.1 if T_min is None else T_min
        T_max = Tcrit + 15.0 if T_max is None else T_max
        colors = self._iso_colors if colors is None else colors
        for i, p in enumerate(pressures):
            s_iso, T_iso = self._isobar(fluid, p, T_min, T_max, pcrit)
            ax.plot(s_iso / 1e3, T_iso, ls=ls, lw=lw,
                    color=colors[i % len(colors)],
                    label=f"p = {p / 1e3:.0f} kPa", zorder=3)
        return ax

    def plot_ts_cycle(self, ax, cycle : Cycle, color=BLUE, textcoords="offset points",
                        xytext=(-14, 7), fontweight="bold", linestyle="-", 
                        lw=2.0, label="Cycle"):
        fluid = cycle.fluid
        s_23, T_23 = self._isobar_segment(fluid, cycle.pc, cycle.s2, cycle.s3)
        s_41, T_41 = self._isobar_segment(fluid, cycle.pe, cycle.s4, cycle.s1)
        path_s = np.concatenate([[cycle.s1, cycle.s2], s_23, [cycle.s4], s_41])
        path_T = np.concatenate([[cycle.T1, cycle.T2], T_23, [cycle.T4], T_41])
        ax.plot(path_s / 1e3, path_T, linestyle=linestyle, color=color, lw=lw, label=label, zorder=5)

        # State points and their labels.
        s_states = np.array([cycle.s1, cycle.s2, cycle.s3, cycle.s4])
        T_states = np.array([cycle.T1, cycle.T2, cycle.T3, cycle.T4])
        ax.plot(s_states / 1e3, T_states, "o", color=color, markersize=6, zorder=6)
        for i, (s, T) in enumerate(zip(s_states, T_states, strict=True), start=1):
            ax.annotate(str(i), (s / 1e3, T), textcoords=textcoords,
                        xytext=xytext, fontweight=fontweight, zorder=6)
        return ax

    @staticmethod
    def _isotherm(fluid, T, p_min, p_max, Tcrit):
        def h_at(p):
            fluid.update(CP.PT_INPUTS, p, T)
            return fluid.hmass()

        if T < Tcrit:
            # Saturation pressure and the enthalpies bounding the two-phase gap.
            fluid.update(CP.QT_INPUTS, 0.0, T)
            p_sat, h_l = fluid.p(), fluid.hmass()
            fluid.update(CP.QT_INPUTS, 1.0, T)
            h_v = fluid.hmass()

            # Subcooled liquid: high pressure down to just above saturation.
            p_sub = np.linspace(p_max, p_sat * 1.001, 40)
            h_sub = np.array([h_at(p) for p in p_sub])
            # Superheated vapour: just below saturation down to p_min.
            p_sup = np.linspace(p_sat * 0.999, p_min, 40)
            h_sup = np.array([h_at(p) for p in p_sup])

            h = np.concatenate([h_sub, [h_l, h_v], h_sup])
            p = np.concatenate([p_sub, [p_sat, p_sat], p_sup])
        else:
            # Supercritical: a single smooth single-phase line.
            p = np.linspace(p_max, p_min, 80)
            h = np.array([h_at(pi) for pi in p])
        return h, p

    @staticmethod
    def _isobar(fluid, p, T_min, T_max, pcrit):
        def s_at(T):
            fluid.update(CP.PT_INPUTS, p, T)
            return fluid.smass()

        if p < pcrit:
            # Saturation temperature and the entropies bounding the two-phase gap.
            fluid.update(CP.PQ_INPUTS, p, 0.0)
            T_sat, s_l = fluid.T(), fluid.smass()
            fluid.update(CP.PQ_INPUTS, p, 1.0)
            s_v = fluid.smass()

            segments_s, segments_T = [], []
            # Subcooled liquid: T_min up to just below saturation.
            if T_sat - 0.1 > T_min:
                T_sub = np.linspace(T_min, T_sat - 0.1, 40)
                segments_s.append(np.array([s_at(T) for T in T_sub]))
                segments_T.append(T_sub)
            # Horizontal jump across the dome at the saturation temperature.
            segments_s.append(np.array([s_l, s_v]))
            segments_T.append(np.array([T_sat, T_sat]))
            # Superheated vapour: just above saturation up to T_max.
            T_sup = np.linspace(T_sat + 0.1, T_max, 40)
            segments_s.append(np.array([s_at(T) for T in T_sup]))
            segments_T.append(T_sup)

            s = np.concatenate(segments_s)
            T = np.concatenate(segments_T)
        else:
            # Supercritical: a single smooth single-phase line.
            T = np.linspace(T_min, T_max, 80)
            s = np.array([s_at(Ti) for Ti in T])
        return s, T

    @staticmethod
    def _isobar_segment(fluid, p, s_start, s_end, n=60):
        s = np.linspace(s_start, s_end, n)
        T = np.empty_like(s)
        for i, si in enumerate(s):
            fluid.update(CP.PSmass_INPUTS, p, si)
            T[i] = fluid.T()
        return s, T

if __name__ == "__main__":
    Te = 0 + 273.15  # [K] evaporation temperature on the low pressure side
    Tc = 30 + 273.15 # [K] condensation temperature on the high pressure side
    sh = 5           # [K] superheating at compressor inlet
    sc = 0           # [K] subcooling at condensor outlet

    # Refrigerant used in the refrigeration cycle
    fluid_ref = AbstractState("HEOS", "R134a")

    cycle = Cycle(fluid_ref, Te=0+273.15, Tc=30+273.15, sh=5.0, sc=5.0)

    # using water in the secondary cycles of the evaporators and condensers
    fluid_sec_evap = AbstractState("HEOS", "Water")
    fluid_sec_cond = AbstractState("HEOS", "Water")

    # Efficiency-based compressor parameters
    is_eff = 0.8  # Isentropic efficiency
    vol_eff = 0.9 # Volumetric efficiency
    N = 1500 / 60 # [Hz] Compressor frequency (1500 rpm / 60)
    Vd = 5e-5     # [m³] Displaced volume per revolution

    # Boundary condition at heat exchangers
    m_flow_sec_cond = 0.1       # [kg/s] Mass flow rate at the inlet of the secondary cycle of the condenser
    p_sec_cond = 1e5            # [Pa] Pressure in the secondary cycle of the condenser
    T_sec_in_cond = 20+273.15   # [K] Temperature at the inlet of the secondary cycle of the condenser
    m_flow_sec_evap = 0.1       # [kg/s] Mass flow rate at the inlet of the secondary cycle of the evaporator
    p_sec_evap = 1e5            # [Pa] Pressure in the secondary cycle of the evaporator
    T_sec_in_evap = 10+273.15   # [K] Temperature at the inlet of the secondary cycle of the evaporator

    eff_comp = EffCompressor(fluid_ref, is_eff=is_eff, vol_eff=vol_eff, N=N, Vd=Vd)
    evaporator = BasicHeatExchanger(fluid_ref, fluid_sec_evap)
    condenser = BasicHeatExchanger(fluid_ref, fluid_sec_cond)

    eff_comp = EffCompressor(fluid_ref, is_eff=is_eff, vol_eff=vol_eff, N=N, Vd=Vd)
    evaporator = BasicHeatExchanger(fluid_ref, fluid_sec_evap)
    condenser = BasicHeatExchanger(fluid_ref, fluid_sec_cond)

    cycle.set_compressor(eff_comp)
    cycle.set_condenser(condenser)
    cycle.set_evaporator(evaporator)

    cycle.calc(m_flow_sec_cond=m_flow_sec_cond, p_sec_cond=p_sec_cond, T_sec_in_cond=T_sec_in_cond, m_flow_sec_evap=m_flow_sec_evap, p_sec_evap=p_sec_evap, T_sec_in_evap=T_sec_in_evap)

    plotter = CyclePlotter()

    fig_ph = plt.figure(figsize=(8, 6))
    ax_ph = plt.axes()
    plotter.plot_ph(ax_ph, cycle)
    fig_ph.tight_layout()

    fig_ts = plt.figure(figsize=(8, 6))
    ax_ts = plt.axes()
    plotter.plot_ts(ax_ts, cycle)
    fig_ts.tight_layout()

    plt.show()
