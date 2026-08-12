# Basic heat exchanger
`BasicHeatExchanger` only enforces an energy balance between the refrigerant and the secondary fluid. It does not check wether one of the fluids consistently stays hotter then the other side of the fluid, allowing the non-physical flow of heat from a cool to a hot fluid.

![HeatExchnager](../assets/heat_exchanger_render.svg){ width="600" }

## Duty

The duty is fixed directly by the refrigerant enthalpy change and its mass flow. The enthalpies of the refrigerant at the inlet and outlet and th emass flow are provided through boundary conditions

$$
Q = \dot m_{\mathrm{ref}} \,(h_{\mathrm{ref,in}} - h_{\mathrm{ref,out}})
$$

$Q$ is positive when the refrigerant releases heat (condenser, $h_{\mathrm{ref,in}} >
h_{\mathrm{ref,out}}$) and negative when it absorbs heat (evaporator). The base class uses
the same sign to set `state` to `"cond"` or `"evap"`.

## Secondary-side outlet

The secondary fluid receives or supplies exactly that duty. The inlet enthalpy and the mass flow of the secondary fluid are known.

$$
h_{\mathrm{sec,out}} = h_{\mathrm{sec,in}} + \frac{Q}{\dot m_{\mathrm{sec}}}
$$

The secondary fluid heats up when $Q>0$ (condenser) and cools down when $Q<0$ (evaporator).

