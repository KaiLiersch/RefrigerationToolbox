"""Tests for `refrigerationtoolbox` package."""

import refrigerationtoolbox
import refrigerationtoolbox.cycle as cycle


def test_import():
    """Verify the package can be imported."""
    assert refrigerationtoolbox

def test_cycle():
    """Quick and dirty testcase"""
    assert cycle.cycle("sdkh")
