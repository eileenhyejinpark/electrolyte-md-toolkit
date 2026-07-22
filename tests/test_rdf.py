"""Unit tests for the RDF module."""

import numpy as np
import MDAnalysis as mda

from electrolyte_toolkit.rdf import compute_rdf


def make_two_atom_universe(distance):
    """Build a minimal 2-atom Universe with a known, fixed separation."""
    u = mda.Universe.empty(n_atoms=2, trajectory=True)
    u.add_TopologyAttr("type", ["1", "2"])
    u.dimensions = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    u.atoms.positions = np.array([
        [0.0, 0.0, 0.0],
        [distance, 0.0, 0.0],
    ])
    return u


def test_rdf_peak_at_known_distance():
    """The RDF between two atoms at a fixed known separation should peak
    at (approximately) that separation distance.
    """
    distance = 2.0
    u = make_two_atom_universe(distance)

    bins, rdf = compute_rdf(u, "type 1", "type 2", nbins=50, range=(0.0, 5.0))

    peak_r = bins[np.argmax(rdf)]
    assert np.isclose(peak_r, distance, atol=0.15)  # within one bin width