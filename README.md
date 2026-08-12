# RefrigerationToolbox

Toolbox for modelling refrigeration and heatpump cycles.

* [GitHub](https://github.com/KaiLiersch/RefrigerationToolbox/) | [Documentation](https://KaiLiersch.github.io/RefrigerationToolbox/)
* Created by [Kai Liersch](-) | GitHub [@KaiLiersch](https://github.com/KaiLiersch)
* MIT License

## Installation

Currently, this package is not hosted on any python package index and can only
be installed from source.

### Install from source

Simply clone from github:

```sh
git clone git@github.com:KaiLiersch/RefrigerationToolbox.git
```

Enter cloned repository and run:
```sh
cd RefrigerationToolbox
pip install .
```

If you just want to execute the jupyter notebooks I recommend using the approach described in the following section,
as the setup with uv automatically creates a virtual environment and lets you use the virtual environment as a jupyter environment.

### Setup development environment
If you want to make changes to the package itself, follow this workflow:

```sh
git clone git@github.com:KaiLiersch/RefrigerationToolbox.git
```
Enter the cloned repository and set up uv:
```sh
cd RefrigerationToolbox
uv sync
```
To install the virtual environment as a jupyter environment run:
```sh
source .venv/bin/activate
python -m ipykernel install --user --name=refrigerationtoolbox-venv --display-name="Python (RefrigerationToolbox .venv)"
```
VSCode will only find the jupyter environment if the .venv folder is in the root folder of the currently open workspace.

## Getting started

I recommend having a look at the jupyter notebooks in the examples folder.
For a more in-depth look check out the [documentation](https://KaiLiersch.github.io/RefrigerationToolbox/).


## Features

* Model refrigeration and heat pump cycles and evaluate performance (COP, heating/cooling capacity, and power requirements).
* Pick from two compressor models: efficiency-based and polynomial-based (AHRI 540 / DIN EN 12900).
* Model simple heat exchangers, or plate heat exchangers that can size itself and determine pressure drops
* Optimize for a given heat exchanger geometry COP.
* Simulate how the cycle behaves over time as boundary conditions change.
* Plot pressure-enthalpy and temperature-entropy diagrams.
* Fluid properties come from CoolProp, so many refrigerants are supported.

## Documentation

Documentation is built with [Zensical](https://zensical.org/) and deployed to GitHub Pages.

* **Live site:** https://KaiLiersch.github.io/RefrigerationToolbox/
* **Preview locally:** `just docs-serve` (serves at http://localhost:8000)
* **Build:** `just docs-build`

## Author

RefrigerationToolbox was created in 2026 by Kai Liersch.

Built with [Cookiecutter](https://github.com/cookiecutter/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
