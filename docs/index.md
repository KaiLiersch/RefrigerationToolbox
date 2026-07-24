# RefrigerationToolbox

Toolbox for modelling refrigeration and heatpump cycles.

## Features

* Model refrigeration and heat pump cycles and evaluate performance (COP, heating/cooling capacity, and power requirements).
* Pick from two compressor models: efficiency based and polynomial based (AHRI 540 / DIN EN 12900).
* Model simple heat exchanger, or plate heat exchangers that can size itself and determine pressure drops
* Optimize for a given heat exchanger geometry COP.
* Simulate how the cycle behaves over time as boundary conditions change.
* Plot pressure-enthalpy and temperature-entropy diagrams.
* Fluid properties come from CoolProp, so many refrigerants are supported.

## Outlook
These are the changes that I have planned:
 * Clear error handling for non physical behavior in the HeatExchanger models
 * Documentation of implemented physical models
 * More extensive examples 

## Getting started

- [Installation](installation.md) - how to install RefrigerationToolbox
- [Usage](usage.md) - how to use RefrigerationToolbox
- [API Reference](api.md) - auto-generated API documentation
