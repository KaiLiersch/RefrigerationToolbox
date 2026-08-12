# Polynomial compressor

`PolynomialCompressor` predicts the mass flow $\dot m$, cooling capacity $Q_{\mathrm{e}}$
and input power $P_{\mathrm{el}}$ directly from manufacturer rating polynomials (AHRI 540 / EN 12900).

## Rating polynomials

Each quantity is modelled by a bi-cubic polynomial in the evaporation and condensation temperatures,
with ten coefficients $a_1, \dots, a_{10}$:

$$
f(T_\mathrm{e}, T_\mathrm{c}) = a_1 + a_2 T_\mathrm{e} + a_3 T_\mathrm{c} + a_4 T_\mathrm{e}^2
+ a_5 T_\mathrm{e} T_\mathrm{c} + a_6 T_\mathrm{c}^2 + a_7 T_\mathrm{e}^3
+ a_8 T_\mathrm{e}^2 T_\mathrm{c} + a_9 T_\mathrm{e} T_\mathrm{c}^2 + a_{10} T_\mathrm{c}^3
$$

Three independent coefficient sets evaluated at the same $(T_\mathrm{e}, T_\mathrm{c})$ give the mass flow, cooling capacity and power for the operating point. Coefficients must be in Kelvin. Use [`celsius_to_kelvin_coeffs`](../api.md) to convert Celsius-based polynomials. Note that, the coefficients change for different operating points making this model unfit for the optimization of the operating point.

## Heating capacity

By an energy balance over the compressor,
the heating capacity is the sum of power and cooling capacity:

$$
Q_\mathrm{c} = P + Q_\mathrm{e}
$$

## Suction state

Same as [EffCompressor](efficiency_compressor.md) — superheated vapour at $p_\mathrm{e}$:

$$
T_{\mathrm{in}} = T_\mathrm{e} + \Delta T_{\mathrm{sh}}, \qquad p_{\mathrm{in}} = p_{\mathrm{sat}}(T_\mathrm{e})
$$

## Discharge state

The polynomials do not directly provide the compressor discharge state, so the discharge
enthalpy $h_{\mathrm{out}}$ is reconstructed from the condenser-outlet enthalpy $h_3$ and the
specific heating capacity $q_\mathrm{c} = Q_\mathrm{c} / \dot m$:

$$
h_{\mathrm{out}} = h_3 + q_\mathrm{c}
$$

$h_3$ corresponds to the enthalpy at the condenser outlet in [Cycle](cycle.md). $T_{\mathrm{out}}$, $s_{\mathrm{out}}$ and
$\rho_{\mathrm{out}}$ then follow from $h_{\mathrm{out}}$ at $p_\mathrm{c} = p_{\mathrm{sat}}(T_\mathrm{c})$.

## Effects on overall cycles

Because $\dot m$, $Q_\mathrm{e}$ and $P$ come from three independently curve-fitted
polynomials rather than a single equation of state, the resulting states are only
*approximately* energy-consistent. This is why [Cycle](cycle.md) uses the polynomial duties
directly, and derives the evaporator-inlet enthalpy from $Q_\mathrm{e}$ rather than
assuming a strictly isenthalpic expansion, when a `PolynomialCompressor` is attached.
