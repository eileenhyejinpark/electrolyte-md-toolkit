"""Diffusion coefficient calculations from mean-squared displacement (MSD)."""

import numpy as np
from MDAnalysis.analysis.msd import EinsteinMSD


def compute_msd(universe, select="all", msd_type="xyz", fft=True):
    """Compute the time-averaged MSD for a selection of atoms.

    Requires a Universe built from UNWRAPPED coordinates - wrapped/periodic
    coordinates will produce incorrect (artificially suppressed) MSD values.

    Parameters
    ----------
    universe : MDAnalysis.Universe
    select : str
        MDAnalysis selection string (e.g. "type 2").
    msd_type : str
        Dimensionality to compute MSD over ("xyz", "xy", "x", etc.).
    fft : bool
        Use the (much faster) FFT-based algorithm.

    Returns
    -------
    np.ndarray
        MSD as a function of lag time (in frames).
    """
    msd_analysis = EinsteinMSD(universe, select=select, msd_type=msd_type, fft=fft)
    msd_analysis.run()
    return msd_analysis.results.timeseries


def diffusion_coefficient(msd, timestep, dim=3, fit_start=0.1, fit_end=0.9):
    """Estimate the diffusion coefficient from the linear region of the MSD
    via the Einstein relation: MSD(t) = 2 * dim * D * t

    Parameters
    ----------
    msd : np.ndarray
        MSD values from compute_msd, one per trajectory frame.
    timestep : float
        Time between frames (dump frequency * simulation timestep).
    dim : int
        Number of dimensions the MSD was computed over (3 for "xyz").
    fit_start, fit_end : float
        Fractional bounds (0-1) of the MSD curve to fit linearly, excluding
        the short-time ballistic regime and the noisy long-time tail.

    Returns
    -------
    float
        Estimated diffusion coefficient D.
    """
    n = len(msd)
    lags = np.arange(n) * timestep
    start = int(fit_start * n)
    end = int(fit_end * n)

    slope, intercept = np.polyfit(lags[start:end], msd[start:end], 1)
    return slope / (2 * dim)