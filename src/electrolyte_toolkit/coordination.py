"""Coordination number calculations from radial distribution functions."""

import numpy as np


def first_minimum(bins, rdf):
    """Find the first local minimum of g(r) after its first peak.

    This is the standard, physically motivated cutoff radius for defining
    a "first solvation shell" when computing coordination numbers.
    """
    peak_idx = np.argmax(rdf)
    for i in range(peak_idx, len(rdf) - 1):
        if rdf[i] < rdf[i + 1]:
            return bins[i]
    return bins[-1]


def coordination_number(bins, rdf, number_density, r_max=None):
    """Integrate g(r) up to r_max to get the average coordination number.

    Parameters
    ----------
    bins : np.ndarray
        Bin centers (distance r) from an RDF calculation.
    rdf : np.ndarray
        g(r) values at each bin.
    number_density : float
        Number density of the coordinating species (atoms per unit volume).
    r_max : float, optional
        Integration cutoff. If None, uses the first minimum after the
        first peak of g(r).

    Returns
    -------
    float
        The average coordination number within r_max.
    """
    if r_max is None:
        r_max = first_minimum(bins, rdf)

    mask = bins <= r_max
    integrand = 4 * np.pi * bins[mask] ** 2 * rdf[mask] * number_density
    return np.trapz(integrand, bins[mask])