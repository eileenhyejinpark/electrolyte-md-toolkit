"""Unit tests for the diffusion coefficient module."""

import numpy as np

from electrolyte_toolkit.diffusion import diffusion_coefficient


def test_diffusion_coefficient_recovers_known_slope():
    """Construct a perfectly linear synthetic MSD via the exact Einstein
    relation MSD(t) = 2*dim*D*t, and check that diffusion_coefficient()
    recovers the true D.
    """
    D_true = 0.25
    dim = 3
    timestep = 0.5
    n_frames = 200

    lags = np.arange(n_frames) * timestep
    msd = 2 * dim * D_true * lags  # exact, noise-free linear relationship

    D_estimated = diffusion_coefficient(msd, timestep=timestep, dim=dim)

    assert np.isclose(D_estimated, D_true, rtol=1e-6)