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

### Setup development environment
If you want to make changes to the packages itself, follow this workflow:

```sh
git clone git@github.com:KaiLiersch/RefrigerationToolbox.git
```
Enter the cloned repository and setup uv:
```sh
cd RefrigerationToolbox
uv sync
```

## Features

* Build a refrigeration or heat pump cycle and work out how well it performs (COP, heat, and power).
* Pick from two compressor models: one based on efficiencies, one based on manufacturer data.
* Use a simple heat exchanger, or a detailed plate heat exchanger that also sizes itself and checks pressure drop.
* Automatically find the operating point that gives the best COP.
* Simulate how the cycle behaves over time as conditions change.
* Draw the common diagrams (pressure-enthalpy and temperature-entropy) and plots of the results.
* Fluid properties come from CoolProp, so many refrigerants are supported.

## Documentation

Documentation is built with [Zensical](https://zensical.org/) and deployed to GitHub Pages.

* **Live site:** https://KaiLiersch.github.io/RefrigerationToolbox/
* **Preview locally:** `just docs-serve` (serves at http://localhost:8000)
* **Build:** `just docs-build`


## Author

RefrigerationToolbox was created in 2026 by Kai Liersch.

Built with [Cookiecutter](https://github.com/cookiecutter/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
