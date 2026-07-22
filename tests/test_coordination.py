"""Unit tests for the coordination number module."""

import numpy as np


from electrolyte_toolkit.coordination import coordination_number, first_minimum


def test_first_minimum_simple_peak():
    """first_minimum should locate the dip after a clear peak in g(r)."""
    bins = np.linspace(0, 5, 100)
    rdf = np.ones_like(bins)
    rdf[15:25] = 3.0   # peak
    rdf[35:45] = 0.5   # minimum after the peak

    r_min = first_minimum(bins, rdf)
    assert bins[35] <= r_min <= bins[45]


def test_coordination_number_ideal_gas():
    """For g(r) = 1 everywhere (uncorrelated/ideal gas), the coordination
    number has an exact analytical value: N(r_max) = (4/3) * pi * r_max^3 * density.
    """
    bins = np.linspace(0.001, 5, 2000)
    rdf = np.ones_like(bins)
    density = 0.05
    r_max = 3.0

    cn = coordination_number(bins, rdf, density, r_max=r_max)
    expected = (4 / 3) * np.pi * r_max**3 * density

    assert np.isclose(cn, expected, rtol=0.01)