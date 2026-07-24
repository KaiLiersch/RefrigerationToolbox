from typing import Any

import numpy as np
from CoolProp.CoolProp import AbstractState
from scipy.optimize import OptimizeResult, minimize

from refrigerationtoolbox.cycle.Compressor import Compressor
from refrigerationtoolbox.cycle.Cycle import Cycle
from refrigerationtoolbox.cycle.EffCompressor import EffCompressor
from refrigerationtoolbox.cycle.PlateHeatExchanger import PlateHeatExchanger

# The secondary side boundary conditions, in the order they are passed to Cycle.calc(...):
# m_flow_sec_cond, p_sec_cond, T_sec_in_cond, m_flow_sec_evap, p_sec_evap, T_sec_in_evap
BoundaryConditions = tuple[float, float, float, float, float, float]


class PlateHeatExchangerCycleOptimizer:
    """Maximizes the heating COP of a cycle over its evaporation and condensation temperatures.

    Wraps a :class:`Cycle` built from a compressor and two fixed-geometry plate heat
    exchangers and searches for the ``(Te, Tc)`` pair that gives the highest heating
    coefficient of performance for a given set of secondary-side boundary conditions. The
    search uses SLSQP subject to inequality constraints that keep both exchangers feasible:
    each pinch temperature difference must stay above its limit and each stream pressure
    drop and required area must stay within the installed values. The pressure-drop and
    area constraints are normalized by their limits so they share one tolerance with the
    temperature pinches. Because the exchanger discretization changes the element count in
    steps, the objective and constraints are only piecewise smooth, which is handled with a
    finite-difference step large enough to step over that noise and with caching of the
    most recent evaluation. After a solve the optimum is recomputed once so the attached
    compressor and exchangers are left in a consistent state.

    Args:
        refrigerant: CoolProp state object for the working fluid.
        Te_start: Initial guess for the evaporation temperature [K]; overwritten with the
            solution after :meth:`optimize_cop` runs.
        Tc_start: Initial guess for the condensation temperature [K]; overwritten likewise.
        sh: Superheat at the evaporator outlet [K].
        sc: Subcooling at the condenser outlet [K].
        compressor: Compressor model for the cycle.
        evaporator: Plate heat exchanger acting as the evaporator.
        condenser: Plate heat exchanger acting as the condenser.
    """

    def __init__(
        self,
        refrigerant: AbstractState,
        Te_start: float,
        Tc_start: float,
        sh: float = 0.0,
        sc: float = 0.0,
        compressor: Compressor | None = None,
        evaporator: PlateHeatExchanger | None = None,
        condenser: PlateHeatExchanger | None = None,
    ) -> None:

        self.sh = sh
        self.sc = sc

        self.fluid = refrigerant
        self.compressor = compressor
        self.evaporator = evaporator
        self.condenser = condenser

        # Starting values for the optimization, overwritten with the solution by optimize_cop(...)
        self.Te = Te_start
        self.Tc = Tc_start

        self.cycle: Cycle | None = None
        self.result: OptimizeResult | None = None
        self.constraint_values: dict[str, float] | None = None
        self.feas: bool = False

        self._cache_key: tuple[float, float] | None = None
        self._cache: dict[str, Any] | None = None

    def optimize_cop(
        self,
        m_flow_sec_cond: float,
        p_sec_cond: float,
        T_sec_in_cond: float,
        m_flow_sec_evap: float,
        p_sec_evap: float,
        T_sec_in_evap: float,
        Te_bounds: tuple[float, float] | None = None,
        Tc_bounds: tuple[float, float] | None = None,
        eps: float = 5e-2,
        ftol: float = 1e-8,
        feas_tol: float = 1e-6,
        verbose: bool = False,
    ) -> Cycle:
        """Find the evaporation and condensation temperatures that maximize the heating COP.

        Runs an SLSQP search over ``(Te, Tc)``, starting from the current stored values and
        subject to the exchanger feasibility constraints. On return ``Te``, ``Tc``,
        ``cycle``, ``result``, ``constraint_values`` and ``feas`` are updated, and the
        optimum is recomputed once so the compressor and exchangers are left consistent.

        Args:
            m_flow_sec_cond: Condenser secondary mass flow rate [kg/s].
            p_sec_cond: Condenser secondary pressure [Pa].
            T_sec_in_cond: Condenser secondary inlet temperature [K].
            m_flow_sec_evap: Evaporator secondary mass flow rate [kg/s].
            p_sec_evap: Evaporator secondary pressure [Pa].
            T_sec_in_evap: Evaporator secondary inlet temperature [K].
            Te_bounds: Optional (lower, upper) bounds for Te [K]; derived from the secondary
                inlet temperatures and pinch when omitted.
            Tc_bounds: Optional (lower, upper) bounds for Tc [K]; derived likewise.
            eps: Absolute finite-difference step for the optimizer [K], large enough to step
                over the discretization's piecewise-constant noise.
            ftol: Convergence tolerance passed to SLSQP.
            feas_tol: Tolerance by which each normalized constraint may be negative and the
                design still counts as feasible.
            verbose: Enables diagnostic output during evaluation.

        Returns:
            The recomputed cycle at the optimum.
        """
        self._check_components()

        bc = (m_flow_sec_cond, p_sec_cond, T_sec_in_cond, m_flow_sec_evap, p_sec_evap, T_sec_in_evap)
        lower, upper = self._bounds(bc, Te_bounds, Tc_bounds)

        x0 = np.clip(np.array([self.Te, self.Tc]), lower + 1e-3, upper - 1e-3)

        self.result = minimize(
            lambda x: -self._evaluate(x, bc, verbose)["COP_heat"],
            x0,
            method="SLSQP",
            bounds=list(zip(lower, upper, strict=True)),
            constraints=[{"type": "ineq", "fun": lambda x: self._constraints(x, bc, verbose)}],
            # Absolute finite difference step in K. The number of elements per phase region in the
            # heat exchanger discretization changes in steps, so objective and constraints are only
            # piecewise smooth and a step of 0.05 K is needed to step over that noise.
            options={"eps": eps, "ftol": ftol},
        )

        self.Te, self.Tc = self.result.x

        # The last evaluation of minimize is not necessarily the optimum, so recalculate there to
        # leave the compressor and the two heat exchangers in a consistent state
        self.cycle = self._calc_cycle(self.Te, self.Tc, bc, verbose)

        constraints = self._constraints(self.result.x, bc, verbose)
        self.constraint_values = dict(zip(self._constraint_names(), constraints, strict=True))
        # Determined from the constraints rather than from PlateHeatExchanger.feas, because the
        # constraints are normalized and therefore share one tolerance across temperatures and pressures
        self.feas = bool(np.all(constraints >= -feas_tol))

        return self.cycle

    def _constraints(self, x: np.ndarray, bc: BoundaryConditions, verbose: bool) -> np.ndarray:
        """Feasibility constraints of both heat exchangers, all of them >= 0 when satisfied.

        The pressure drops are normalized with their limit so that all constraints are of the same
        order of magnitude, otherwise the pressure drops in Pa would dominate the pinch points in K.
        """
        return self._evaluate(x, bc, verbose)["constraints"]

    @staticmethod
    def _constraint_names() -> tuple[str, ...]:
        """Return the labels of the constraint vector, in the same order as :meth:`_evaluate`."""
        return (
            "condenser dT_min",
            "evaporator dT_min",
            "condenser dp_ref",
            "condenser dp_sec",
            "evaporator dp_ref",
            "evaporator dp_sec",
            "condenser A_phex",
            "evaporator A_phex",
        )

    def _evaluate(self, x: np.ndarray, bc: BoundaryConditions, verbose: bool) -> dict[str, Any]:
        """Solve the cycle at a candidate (Te, Tc) and return its COP and constraint vector.

        Wraps :meth:`_calc_cycle` and packages the heating COP together with the eight
        exchanger constraints (two pinch differences, four normalized pressure-drop margins,
        two normalized area margins). The result is cached under the (Te, Tc) key so the
        optimizer's repeated objective and constraint calls at the same point reuse one
        solve. A non-finite or failed evaluation (for example a temperature cross) is turned
        into a strongly infeasible sentinel so the optimizer steps away from it.

        Args:
            x: Candidate temperatures ``[Te, Tc]`` [K].
            bc: The six secondary-side boundary conditions.
            verbose: Enables diagnostic output.

        Returns:
            Dict with ``"COP_heat"`` and the ``"constraints"`` array.
        """
        key = (float(x[0]), float(x[1]))
        if self._cache_key == key and self._cache is not None:
            return self._cache

        try:
            # A temperature cross produces NaN rather than an exception, which numpy reports as a
            # warning. The result is checked for finite values below, so the warning is only noise.
            with np.errstate(invalid="ignore", divide="ignore"):
                cycle = self._calc_cycle(key[0], key[1], bc, verbose)
            values = {
                "COP_heat": cycle.COP_heat,
                "constraints": np.array(
                    [
                        self.condenser.dT_min - self.condenser.dT_pinch,
                        self.evaporator.dT_min - self.evaporator.dT_pinch,
                        (self.condenser.dp_ref_max - self.condenser.dp_ref) / self.condenser.dp_ref_max,
                        (self.condenser.dp_sec_max - self.condenser.dp_sec) / self.condenser.dp_sec_max,
                        (self.evaporator.dp_ref_max - self.evaporator.dp_ref) / self.evaporator.dp_ref_max,
                        (self.evaporator.dp_sec_max - self.evaporator.dp_sec) / self.evaporator.dp_sec_max,
                        (self.condenser.A_available - self.condenser.A_phex) / self.condenser.A_available,
                        (self.evaporator.A_available - self.evaporator.A_phex) / self.evaporator.A_available,
                    ]
                ),
            }
            if not (np.isfinite(values["COP_heat"]) and np.all(np.isfinite(values["constraints"]))):
                raise ArithmeticError("the cycle could not be evaluated to finite values")
        except Exception:
            values = {"COP_heat": 0.0, "constraints": np.full(8, -1e3)}

        self._cache_key = key
        self._cache = values
        return values

    def _check_components(self) -> None:
        """Raise if the compressor, evaporator or condenser has not been set.

        Raises:
            RuntimeError: When any of the three cycle components is missing.
        """
        if self.compressor is None:
            raise RuntimeError("self.compressor is None. Please set a compressor using the constructor")
        if self.evaporator is None:
            raise RuntimeError("self.evaporator is None. Please set an evaporator using the constructor")
        if self.condenser is None:
            raise RuntimeError("self.condenser is None. Please set a condenser using the constructor")

    def _bounds(self, bc: BoundaryConditions, Te_bounds: tuple[float, float] | None, Tc_bounds: tuple[float, float] | None) -> tuple[np.ndarray, np.ndarray]:
        """Build the lower and upper search bounds for (Te, Tc).

        When a bound is not supplied it is derived from the secondary inlet temperatures and
        the exchanger pinch: Te must stay below the evaporator source temperature (and above
        the fluid's triple point), Tc above the condenser sink temperature (and below the
        critical point).

        Args:
            bc: The six secondary-side boundary conditions.
            Te_bounds: Optional explicit (lower, upper) bounds for Te [K].
            Tc_bounds: Optional explicit (lower, upper) bounds for Tc [K].

        Returns:
            Tuple of the lower and upper bound arrays ``[Te, Tc]``.

        Raises:
            ValueError: If the resulting search interval is empty.
        """
        T_sec_in_cond, T_sec_in_evap = bc[2], bc[5]

        if Te_bounds is None:
            Te_bounds = (
                max(T_sec_in_evap - 50.0, self.fluid.Ttriple() + 1.0),
                T_sec_in_evap - self.evaporator.dT_pinch,
            )
        if Tc_bounds is None:
            Tc_bounds = (
                T_sec_in_cond + self.condenser.dT_pinch,
                min(T_sec_in_cond + 60.0, self.fluid.T_critical() - 1.0),
            )

        lower = np.array([Te_bounds[0], Tc_bounds[0]])
        upper = np.array([Te_bounds[1], Tc_bounds[1]])
        if np.any(lower >= upper):
            raise ValueError(f"Empty search interval for Te={Te_bounds} and/or Tc={Tc_bounds}. Check the secondary fluid inlet temperatures and the pinch temperature differences.")

        return lower, upper

    def _calc_cycle(self, Te: float, Tc: float, bc: BoundaryConditions, verbose: bool) -> Cycle:
        """Build and solve a :class:`Cycle` at the given temperatures and boundary conditions.

        Args:
            Te: Evaporation temperature [K].
            Tc: Condensation temperature [K].
            bc: The six secondary-side boundary conditions.
            verbose: Unused placeholder kept for a uniform call signature.

        Returns:
            The solved cycle, sharing this optimizer's compressor and exchangers.
        """
        m_flow_sec_cond, p_sec_cond, T_sec_in_cond, m_flow_sec_evap, p_sec_evap, T_sec_in_evap = bc

        cycle = Cycle(
            self.fluid,
            Te,
            Tc,
            self.sh,
            self.sc,
            compressor=self.compressor,
            evaporator=self.evaporator,
            condenser=self.condenser,
        )

        cycle.calc(
            m_flow_sec_cond=m_flow_sec_cond,
            p_sec_cond=p_sec_cond,
            T_sec_in_cond=T_sec_in_cond,
            m_flow_sec_evap=m_flow_sec_evap,
            p_sec_evap=p_sec_evap,
            T_sec_in_evap=T_sec_in_evap,
        )

        return cycle

    def __str__(self) -> str:
        """Return a summary of the optimum temperatures, convergence, constraints and cycle."""

        def fmt(val: float | None, scale: float = 1) -> str:
            return "N/A" if val is None else f"{val / scale:.2f}"

        str_rep = f"PlateHeatExchangerCycleOptimizer (refrigerant={self.fluid.fluid_names()[0]}, "
        str_rep += f"sh={self.sh}, sc={self.sc})\n"
        str_rep += f"   Te={fmt(self.Te)} K, Tc={fmt(self.Tc)} K"
        if self.result is None:
            str_rep += " (starting values, optimize_cop(...) has not been called yet)"
        else:
            str_rep += f"\n   COP maximization: converged={self.result.success}, COP={-self.result.fun:.2f}, "
            str_rep += f"feasible={self.feas}"
            if not self.feas:
                str_rep += f" ({self.result.message})"
            str_rep += "\n   Constraints (>= 0 is feasible, ~0 is active):"
            for name, value in self.constraint_values.items():
                str_rep += f"\n      {name:<18} {value:+.3e}"
        if self.cycle is not None:
            str_rep += f"\n{self.cycle}"
        return str_rep


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

    cycle = Cycle(fluid_ref, Te=273 + 0, Tc=273 + 30, sh=5, sc=5)

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
    compressor.calc(0 + 273.15, 30 + 273.15, 5, 5)

    cycle.set_condenser(condenser)
    cycle.set_evaporator(evaporator)
    cycle.set_compressor(compressor)
    cycle.calc(1, 1e5, 15 + 273.15, 0.1, 1e5, 10 + 273.15)

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
    optimizer.optimize_cop(
        m_flow_sec_cond=1,
        p_sec_cond=1e5,
        T_sec_in_cond=15 + 273.15,
        m_flow_sec_evap=0.1,
        p_sec_evap=1e5,
        T_sec_in_evap=10 + 273.15,
    )

    print(optimizer)
    print(optimizer.condenser)
    print(optimizer.evaporator)
