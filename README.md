# RefrigerationToolbox

![PyPI version](https://img.shields.io/pypi/v/RefrigerationToolbox.svg)

Toolbox for modelling refrigeration and heatpump cycles.

* [GitHub](https://github.com/KaiLiersch/RefrigerationToolbox/) | [PyPI](https://pypi.org/project/RefrigerationToolbox/) | [Documentation](https://KaiLiersch.github.io/RefrigerationToolbox/)
* Created by [Kai Liersch](-) | GitHub [@KaiLiersch](https://github.com/KaiLiersch) | PyPI [@KaiLiersch](https://pypi.org/user/KaiLiersch/)
* MIT License

## Features

* Currently not working as it is very early in development

## Documentation

Documentation is built with [Zensical](https://zensical.org/) and deployed to GitHub Pages.

* **Live site:** https://KaiLiersch.github.io/RefrigerationToolbox/
* **Preview locally:** `just docs-serve` (serves at http://localhost:8000)
* **Build:** `just docs-build`

API documentation is auto-generated from docstrings using [mkdocstrings](https://mkdocstrings.github.io/).

Docs deploy automatically on push to `main` via GitHub Actions. To enable this, go to your repo's Settings > Pages and set the source to **GitHub Actions**.

## Development

To set up for local development:

```bash
# Clone your fork
git clone git@github.com:your_username/RefrigerationToolbox.git
cd RefrigerationToolbox

# Install in editable mode with live updates
uv tool install --editable .
```

This installs the CLI globally but with live updates - any changes you make to the source code are immediately available when you run `refrigerationtoolbox`.

Run tests:

```bash
uv run pytest
```

Run quality checks (format, lint, type check, test):

```bash
just qa
```

## Author

RefrigerationToolbox was created in 2026 by Kai Liersch.

Built with [Cookiecutter](https://github.com/cookiecutter/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
