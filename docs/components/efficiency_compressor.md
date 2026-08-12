# Efficiency compressor

The performance of the `EffCompressor` is described by two
efficiencies (isentropic and volumetric), a running speed and a displacement volume. It
calculates its inlet and outlet states, the mass flow $\dot m$ and the power $P$ used by [Cycle](cycle.md) for a given operating point defined by $T_e, T_c, \Delta T_{\mathrm{sh}}, \Delta T_{\mathrm{sc}}$.
Subcooling $\Delta T_{\mathrm{sc}}$ is not used here.

## Suction state

The suction state is superheated vapour at pressure $p_e$, analogous to
[Cycle](cycle.md):

$$
T_{\mathrm{in}} = T_\mathrm{e} + \Delta T_{\mathrm{sh}}, \qquad p_{\mathrm{in}} = p_{\mathrm{sat}}(T_\mathrm{e})
$$

## Compression with isentropic efficiency

Compression to the condensing pressure $p_\mathrm{c} = p_{\mathrm{sat}}(T_\mathrm{c})$ is first evaluated
isentropically, i.e. at constant entropy $s_{\mathrm{in}}$:

$$
h_{\mathrm{is}} = h(p_\mathrm{c}, s_{\mathrm{in}})
$$

$h_{\mathrm{is}}$ is the enthalpy the gas would reach if the compressor were reversible and
adiabatic.

A real compressor needs more work than the isentropic case, because of friction, leakage
and other losses. The isentropic efficiency $\eta_{\mathrm{is}}$ scales the ideal
enthalpy rise up to a realistic one:

$$
h_{\mathrm{out}} = h_{\mathrm{in}} + \frac{h_{\mathrm{is}} - h_{\mathrm{in}}}{\eta_{\mathrm{is}}}
$$

## Mass flow with volumetric efficiency

The compressor displaces a volume $V_\mathrm{d}$ per revolution at speed $N$, so the
swept volume flow rate is $N V_d$. Not all of that volume is actually filled with fresh
suction gas. Re-expansion of gas trapped in the clearance volume, leakage past the
pistons/vanes and pressure losses in the valves all reduce it. This is captured by the
volumetric efficiency $\eta_{\mathrm{vol}}$:

$$
\dot m = \eta_{\mathrm{vol}}\, N\, V_\mathrm{d}\, \rho_\mathrm{in}
$$

where $\rho_{\mathrm{in}}$ is the suction-gas density.

## Power

The shaft/electrical power is the enthalpy rise across the compressor times the mass
flow, i.e. the rate of work performed on the refrigerant:

$$
P_{\mathrm{el}} = \dot m\,(h_{\mathrm{out}} - h_{\mathrm{in}})
$$
