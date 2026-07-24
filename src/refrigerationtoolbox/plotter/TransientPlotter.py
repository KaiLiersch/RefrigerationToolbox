from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
from CoolProp.CoolProp import AbstractState
from matplotlib.axes import Axes

from refrigerationtoolbox.cycle.EffCompressor import EffCompressor
from refrigerationtoolbox.cycle.PlateHeatExchanger import PlateHeatExchanger
from refrigerationtoolbox.cycle.TransientCycle import TransientCycle
from refrigerationtoolbox.optimizer.PlateHeatExchangerCycleOptimizer import PlateHeatExchangerCycleOptimizer

RED = "#D81B60"
BLUE = "#1E88E5"
YELLOW = "#FFC107"
GREEN = "#004D40"


class TransientPlotter:
    """Matplotlib drawing helper for a :class:`TransientCycle` time series.

    Renders the per-step results of a solved transient onto caller-supplied axes, such as
    the evaporation and condensation temperatures over time, and marks the points in time
    where the heat-exchanger constraints are violated. The class holds no state of its own;
    each method takes the axes and the transient to draw and returns the axes.
    """

    def __init__(self) -> None:
        pass

    def plot_temperatures(
        self,
        ax: Axes,
        transient: TransientCycle,
        t_scale: float = 3600.0,
        t_label: str = "Time $t$ (h)",
        celsius: bool = True,
        lw: float = 2.0,
        markersize: float = 5,
    ) -> Axes:
        """Plot the evaporation and the condensation temperature over time."""
        t = self._time(transient, t_scale)
        offset = 273.15 if celsius else 0.0
        unit = "°C" if celsius else "K"

        ax.plot(
            t,
            transient.Tc - offset,
            "-o",
            color=RED,
            lw=lw,
            markersize=markersize,
            label="Condensation $T_c$",
            zorder=4,
        )
        ax.plot(
            t,
            transient.Te - offset,
            "-o",
            color=BLUE,
            lw=lw,
            markersize=markersize,
            label="Evaporation $T_e$",
            zorder=4,
        )
        self._mark_infeasible(
            ax,
            t,
            np.concatenate([transient.Te, transient.Tc]) - offset,
            np.concatenate([transient.feas, transient.feas]),
        )

        ax.set_xlabel(t_label)
        ax.set_ylabel(f"Temperature $T$ ({unit})")
        ax.set_title("Evaporation and condensation temperature")
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(loc="best", fontsize=8)
        return ax

    def plot_cop(
        self,
        ax: Axes,
        transient: TransientCycle,
        t_scale: float = 3600.0,
        t_label: str = "Time $t$ (h)",
        cooling: bool = True,
        lw: float = 2.0,
        markersize: float = 5,
    ) -> Axes:
        """Plot the coefficient of performance over time, heating and optionally cooling."""
        t = self._time(transient, t_scale)

        ax.plot(
            t,
            transient.COP_heat,
            "-o",
            color=GREEN,
            lw=lw,
            markersize=markersize,
            label="Heating $COP_{heat}$",
            zorder=4,
        )
        values = [transient.COP_heat]
        if cooling:
            ax.plot(
                t,
                transient.COP_cool,
                "--o",
                color=YELLOW,
                lw=lw,
                markersize=markersize,
                label="Cooling $COP_{cool}$",
                zorder=4,
            )
            values.append(transient.COP_cool)
        self._mark_infeasible(ax, t, np.concatenate(values), np.concatenate([transient.feas] * len(values)))

        ax.set_xlabel(t_label)
        ax.set_ylabel("Coefficient of performance $COP$")
        ax.set_title("Coefficient of performance")
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(loc="best", fontsize=8)
        return ax

    # --- both stacked on a shared time axis --- #

    def plot_transient(
        self,
        axes: Sequence[Axes] | np.ndarray,
        transient: TransientCycle,
        t_scale: float = 3600.0,
        t_label: str = "Time $t$ (h)",
    ) -> Sequence[Axes] | np.ndarray:
        """Plot the temperatures and the coefficient of performance on two axes sharing the time axis."""
        if len(axes) < 2:
            raise ValueError(f"plot_transient needs two axes, got {len(axes)}")

        self.plot_temperatures(axes[0], transient, t_scale=t_scale, t_label=t_label)
        self.plot_cop(axes[1], transient, t_scale=t_scale, t_label=t_label)
        # The upper plot shares the time axis with the lower one, so only the lower one is labelled
        axes[0].set_xlabel("")
        return axes

    @staticmethod
    def _time(transient: TransientCycle, t_scale: float) -> np.ndarray:
        """Return the transient's time array divided by ``t_scale`` (e.g. 3600 for hours).

        Raises:
            RuntimeError: If the transient has not been solved yet.
        """
        if transient.t is None:
            raise RuntimeError("The transient holds no results. Please call TransientCycle.calc(...) first")
        return transient.t / t_scale

    @staticmethod
    def _mark_infeasible(ax: Axes, t: np.ndarray, values: np.ndarray, feas: np.ndarray, color: str = "black", markersize: float = 9) -> Axes:
        """Cross out the points in time at which the heat exchangers violate their constraints."""
        if np.all(feas):
            return ax
        t_all = np.tile(t, values.size // t.size)
        ax.plot(
            t_all[~feas],
            values[~feas],
            "x",
            color=color,
            markersize=markersize,
            label="Infeasible",
            linestyle="none",
            zorder=6,
        )
        return ax


if __name__ == "__main__":
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
        "Nt": 100,
    }

    fluid_ref = AbstractState("HEOS", "R134a")
    fluid_sec = AbstractState("HEOS", "Water")

    condenser = PlateHeatExchanger(fluid_ref, fluid_sec, geom)
    condenser.determine_min_n_plates(
        U0=2500,
        m_flow_ref=0.01,
        m_flow_sec=1,
        h_ref_in=432500,
        h_ref_out=234500,
        p_ref=770000,
        p_sec=100000,
        T_sec_in=15 + 273.15,
    )

    evaporator = PlateHeatExchanger(fluid_ref, fluid_sec, geom)
    evaporator.determine_min_n_plates(
        U0=2500,
        m_flow_ref=0.01,
        m_flow_sec=0.1,
        h_ref_in=234500,
        h_ref_out=403100,
        p_ref=293000,
        p_sec=100000,
        T_sec_in=10 + 273.15,
    )

    compressor = EffCompressor(fluid_ref, 0.8, 0.9, 30, 0.000025)

    optimizer = PlateHeatExchangerCycleOptimizer(
        fluid_ref,
        Te_start=273 + 0,
        Tc_start=273 + 30,
        sh=5,
        sc=5,
        compressor=compressor,
        evaporator=evaporator,
        condenser=condenser,
    )

    # One day with the source temperature of the evaporator following the ambient temperature
    t = np.linspace(0, 24 * 3600, 25)
    T_sec_in_evap = 10 + 273.15 + 5.0 * np.sin(2 * np.pi * (t / (24 * 3600) - 0.25))

    transient = TransientCycle(optimizer)
    transient.calc(
        t=t,
        m_flow_sec_cond=1,
        p_sec_cond=1e5,
        T_sec_in_cond=15 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=T_sec_in_evap,
    )

    plotter = TransientPlotter()

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    plotter.plot_transient(axes, transient)
    fig.tight_layout()

    plt.show()
