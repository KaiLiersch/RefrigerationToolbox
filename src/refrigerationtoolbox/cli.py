"""Console script for refrigerationtoolbox."""

import typer
from rich.console import Console

from refrigerationtoolbox import utils

app = typer.Typer()
console = Console()


@app.command()
def main() -> None:
    """Console script for refrigerationtoolbox."""
    console.print("Replace this message by putting your code into refrigerationtoolbox.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    utils.do_something_useful()


if __name__ == "__main__":
    app()
