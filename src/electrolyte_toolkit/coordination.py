"""Coordination number calculations from radial distribution functions."""

import numpy as np


def first_minimum(bins, rdf, r_max=None, smooth_window=5):
    """Find the first local minimum of g(r) after its first peak.

    This is the standard, physically motivated cutoff radius for defining
    a "first solvation shell" when computing coordination numbers.
    """
    if r_max is not None:
        mask = bins <= r_max
        bins = bins[mask]
        rdf = rdf[mask]

    if smooth_window > 1:
        kernel = np.ones(smooth_window) / smooth_window
        rdf_smooth = np.convolve(rdf, kernel, mode="same")
    else:
        rdf_smooth = rdf

    peak_idx = np.argmax(rdf_smooth)
    # Global minimum after the peak, within the search window - more robust
    # to noise than searching for the first local uptick.
    tail = rdf_smooth[peak_idx:]
    min_idx = peak_idx + np.argmin(tail)
    return bins[min_idx]


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