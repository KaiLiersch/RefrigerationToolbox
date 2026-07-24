from collections.abc import Sequence
from typing import Self

import numpy as np
from CoolProp.CoolProp import AbstractState

from refrigerationtoolbox.cycle.EffCompressor import EffCompressor
from refrigerationtoolbox.cycle.PlateHeatExchanger import PlateHeatExchanger
from refrigerationtoolbox.optimizer.PlateHeatExchangerCycleOptimizer import PlateHeatExchangerCycleOptimizer

# A boundary condition is either constant over the whole transient or given per point in time
BoundaryCondition = float | Sequence[float] | np.ndarray


class TransientCycle:
    """Quasi-steady time series of an optimized vapour-compression cycle.

    Sweeps a set of boundary conditions over time and, at every point in time, re-optimizes
    the cycle with a :class:`PlateHeatExchangerCycleOptimizer` to find the evaporation and
    condensation temperatures that maximize the coefficient of performance while respecting
    the heat-exchanger limits. The transient is treated as a succession of steady states
    (no thermal storage or dynamics), each warm-started from the previous solution. Each
    boundary condition may be a single value held constant over the whole run or a series
    with one value per point in time. Per-step results (``Te``, ``Tc``, the COPs, duties,
    power, mass flow and a feasibility flag) are stored as arrays aligned with ``t``; steps
    that fail to solve are recorded as ``NaN`` / ``False``.

    Args:
        optimizer: The cycle optimizer used to solve each steady state; its attached
            compressor and heat exchangers define the system being simulated.
    """

    def __init__(self, optimizer : PlateHeatExchangerCycleOptimizer) -> None:
        self.optimizer = optimizer

        self.t = None
        self.Te = None
        self.Tc = None
        self.COP_heat = None
        self.COP_cool = None
        self.Q_heat = None
        self.Q_ref = None
        self.P_comp = None
        self.m_flow = None
        self.feas = None

    def calc(self, t : Sequence[float] | np.ndarray,
             m_flow_sec_cond : BoundaryCondition, p_sec_cond : BoundaryCondition,
             T_sec_in_cond : BoundaryCondition, m_flow_sec_evap : BoundaryCondition,
             p_sec_evap : BoundaryCondition, T_sec_in_evap : BoundaryCondition,
             verbose : bool = False) -> Self:
        """Solve the cycle at every point in time and store the resulting time series.

        Each boundary condition may be a scalar (held constant) or a sequence with one
        value per point in time. For every time step the optimizer is run to find the COP
        optimum for that step's boundary conditions, warm-started from the previous step.
        Steps whose steady state cannot be solved are recorded as ``NaN`` / not feasible and
        do not warm-start the next step. Results are stored on the instance as arrays
        aligned with ``t`` (``Te``, ``Tc``, ``COP_heat``, ``COP_cool``, ``Q_heat``,
        ``Q_ref``, ``P_comp``, ``m_flow``, ``feas``).

        Args:
            t: One-dimensional array of time points [s].
            m_flow_sec_cond: Condenser secondary mass flow rate [kg/s], scalar or per step.
            p_sec_cond: Condenser secondary pressure [Pa], scalar or per step.
            T_sec_in_cond: Condenser secondary inlet temperature [K], scalar or per step.
            m_flow_sec_evap: Evaporator secondary mass flow rate [kg/s], scalar or per step.
            p_sec_evap: Evaporator secondary pressure [Pa], scalar or per step.
            T_sec_in_evap: Evaporator secondary inlet temperature [K], scalar or per step.
            verbose: Forwarded to the optimizer for diagnostic output.

        Returns:
            This instance, so the call can be chained.

        Raises:
            ValueError: If ``t`` is not one-dimensional or a boundary-condition sequence
                does not match the length of ``t``.
        """
        t = np.asarray(t, dtype=float)
        if t.ndim != 1:
            raise ValueError(f"t must be one dimensional, got shape {t.shape}")
        n = t.size

        bc_series = {
            "m_flow_sec_cond": self._as_series(m_flow_sec_cond, n, "m_flow_sec_cond"),
            "p_sec_cond": self._as_series(p_sec_cond, n, "p_sec_cond"),
            "T_sec_in_cond": self._as_series(T_sec_in_cond, n, "T_sec_in_cond"),
            "m_flow_sec_evap": self._as_series(m_flow_sec_evap, n, "m_flow_sec_evap"),
            "p_sec_evap": self._as_series(p_sec_evap, n, "p_sec_evap"),
            "T_sec_in_evap": self._as_series(T_sec_in_evap, n, "T_sec_in_evap"),
        }

        self.t = t
        self.Te = np.full(n, np.nan)
        self.Tc = np.full(n, np.nan)
        self.COP_heat = np.full(n, np.nan)
        self.COP_cool = np.full(n, np.nan)
        self.Q_heat = np.full(n, np.nan)
        self.Q_ref = np.full(n, np.nan)
        self.P_comp = np.full(n, np.nan)
        self.m_flow = np.full(n, np.nan)
        self.feas = np.zeros(n, dtype=bool)

        for i in range(n):
            bc = {name: series[i] for name, series in bc_series.items()}

            Te_last, Tc_last = self.optimizer.Te, self.optimizer.Tc
            try:
                cycle = self.optimizer.optimize_cop(**bc, verbose=verbose)
            except Exception as exc:
                print(f"WARNING: The steady state at t = {t[i]} could not be solved ({exc}). "
                      "The step is stored as NaN.")
                # Restore the last usable temperatures, otherwise the next step is warm started
                # from whatever the failed solve left behind
                self.optimizer.Te, self.optimizer.Tc = Te_last, Tc_last
                continue

            self.Te[i] = self.optimizer.Te
            self.Tc[i] = self.optimizer.Tc
            self.COP_heat[i] = cycle.COP_heat
            self.COP_cool[i] = cycle.COP_cool
            self.Q_heat[i] = cycle.Q_heat
            self.Q_ref[i] = cycle.Q_ref
            self.P_comp[i] = cycle.P_comp
            self.m_flow[i] = cycle.m_flow
            self.feas[i] = self.optimizer.feas

            if not self.optimizer.feas:
                print(f"WARNING: The steady state at t = {t[i]} violates the constraints of the "
                      "heat exchangers. See TransientCycle.feas.")

        return self

    @staticmethod
    def _as_series(value : BoundaryCondition, n : int, name : str) -> np.ndarray:
        """Expand a scalar boundary condition to n points in time or check the length of a sequence."""
        series = np.asarray(value, dtype=float)
        if series.ndim == 0:
            return np.full(n, float(series))
        if series.size != n:
            raise ValueError(f"{name} has {series.size} values, but t has {n} points in time")
        return series

    def __str__(self) -> str:
        """Return a table of the per-step temperatures, heating COP, duty and feasibility."""
        if self.t is None:
            return "TransientCycle (calc(...) has not been called yet)"

        str_rep = f"TransientCycle ({self.t.size} points in time, "
        str_rep += f"t = {self.t[0]:.1f} ... {self.t[-1]:.1f}, {np.count_nonzero(self.feas)} feasible)\n"
        str_rep += f"   {'t':>10} {'Te [°C]':>10} {'Tc [°C]':>10} {'COP_heat':>10} "
        str_rep += f"{'Q_heat [kW]':>12} {'feasible':>9}\n"
        for i in range(self.t.size):
            str_rep += (f"   {self.t[i]:10.1f} {self.Te[i] - 273.15:10.2f} {self.Tc[i] - 273.15:10.2f} "
                        f"{self.COP_heat[i]:10.2f} {self.Q_heat[i] / 1e3:12.2f} {str(self.feas[i]):>9}\n")
        return str_rep.rstrip()


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

    optimizer = PlateHeatExchangerCycleOptimizer(fluid_ref, Te_start=273 + 0, Tc_start=273 + 30, sh=5, sc=5,
                                                 compressor=compressor, evaporator=evaporator, condenser=condenser)

    # One day with the source temperature of the evaporator following the ambient temperature
    t = np.linspace(0, 24 * 3600, 13)
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

    print(transient)
