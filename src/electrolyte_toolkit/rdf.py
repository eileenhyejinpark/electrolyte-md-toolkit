"""Radial distribution function (RDF) calculations for electrolyte systems."""

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF


def compute_rdf(universe, sel1, sel2, nbins=75, range=(0.0, 10.0)):
    """Compute g(r) between two atom selections in an MDAnalysis Universe.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe containing the trajectory.
    sel1, sel2 : str
        MDAnalysis selection strings (e.g. "type 1", "type 2").
    nbins : int
        Number of bins for the histogram.
    range : tuple of float
        (min, max) distance range in Angstroms (or reduced units).

    Returns
    -------
    bins : np.ndarray
        Bin centers (distance r).
    rdf : np.ndarray
        g(r) values at each bin.
    """
    group1 = universe.select_atoms(sel1)
    group2 = universe.select_atoms(sel2)

    rdf_analysis = InterRDF(group1, group2, nbins=nbins, range=range)
    rdf_analysis.run()

    return rdf_analysis.results.bins, rdf_analysis.results.rdf