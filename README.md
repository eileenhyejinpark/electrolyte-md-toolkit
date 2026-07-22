# Electrolyte MD Analysis Toolkit

A Python toolkit for analyzing molecular dynamics trajectories of electrolyte
systems, built on MDAnalysis and LAMMPS. Computes standard structural and
dynamical properties used in battery electrolyte research: radial distribution
functions (RDF), coordination numbers, and diffusion coefficients.

## Status

Core analysis pipeline implemented and validated against a demo molecular
dynamics simulation (see `examples/demo_electrolyte/`): a binary Lennard-Jones
mixture representing a minority "ion" species in a majority "solvent."

## Results (demo system)

**Solvent-ion radial distribution function:**

![RDF plot](examples/demo_electrolyte/rdf_plot.png)

**Ion mean squared displacement:**

![MSD plot](examples/demo_electrolyte/msd_plot.png)


|
 Property 
|
 Value 
|
|
---
|
---
|
|
 RDF peak 
|
 r ≈ 1.1 
|
|
 First solvation shell cutoff 
|
 r = 1.625 
|
|
 Coordination number 
|
 0.965 
|
|
 Diffusion coefficient 
|
 D = 0.170 (reduced LJ units) 
|

Full walkthrough with code: [`examples/demo_electrolyte/analysis_demo.ipynb`](examples/demo_electrolyte/analysis_demo.ipynb)

## Results (TIP3P electrolyte system)

A more realistic system: TIP3P water + Na⁺/Cl⁻ ions (Joung-Cheatham
parameters), PPPM long-range electrostatics, rigid-body water dynamics via
`fix rigid/nvt/small`, thermostatted at 300 K.

![Ion-water RDF](examples/demo_tip3p_electrolyte/ion_water_rdf.png)

| Property | Value |
|---|---|
| Na⁺-O(water) first shell peak | r = 3.08 Å |
| Cl⁻-O(water) first shell peak | r = 5.24 Å |

Coordination numbers were not reliably extractable from this system (only 5
ions of each species, 5 ps of production dynamics — see the analysis
notebook for a full discussion of why). The RDF peaks themselves are robust
across repeated analysis and correctly show Cl⁻'s hydration shell sitting
farther out than Na⁺'s, consistent with expected trends.

Full walkthrough: [`examples/demo_tip3p_electrolyte/analysis_demo.ipynb`](examples/demo_tip3p_electrolyte/analysis_demo.ipynb)

## Setup

```bash
conda env create -f environment.yml
conda activate electrolyte-portfolio
```

## Usage

```python
import MDAnalysis as mda
from electrolyte_toolkit.rdf import compute_rdf
from electrolyte_toolkit.coordination import coordination_number, first_minimum
from electrolyte_toolkit.diffusion import compute_msd, diffusion_coefficient

u = mda.Universe("trajectory.lammpstrj", format="LAMMPSDUMP")

# Radial distribution function
bins, rdf = compute_rdf(u, "type 1", "type 2", nbins=100, range=(0.0, 5.0))

# Coordination number within the first solvation shell
r_cut = first_minimum(bins, rdf)
cn = coordination_number(bins, rdf, number_density, r_max=r_cut)

# Diffusion coefficient (requires an unwrapped-coordinate trajectory)
msd = compute_msd(u, select="type 2")
D = diffusion_coefficient(msd, timestep=frame_dt, dim=3)
```

## Project Structure

- `src/electrolyte_toolkit/` — analysis toolkit source code (`rdf.py`, `coordination.py`, `diffusion.py`)
- `examples/demo_electrolyte/` — demo LAMMPS simulation input, trajectory analysis notebook, and output plots
- `tests/` — unit tests (in progress)

## Roadmap

- [x] Unit tests validating each module against known synthetic data
- [ ] Extend to realistic electrolyte force fields (e.g., TIP3P water + Na⁺/Cl⁻)
      with literature-parameterized pair potentials

## License

MIT