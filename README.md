# Electrolyte MD Analysis Toolkit

A Python toolkit for analyzing molecular dynamics trajectories of electrolyte
systems, built on MDAnalysis and LAMMPS. Computes standard structural and
dynamical properties used in battery electrolyte research: radial distribution
functions (RDF), coordination numbers, and diffusion coefficients.

## Motivation

Built to demonstrate and practice the computational electrolyte analysis
workflows used in molecular dynamics research on battery electrolyte systems
(e.g., ion coordination and solvation structure in liquid electrolytes).

## Status

🚧 In progress — starting with a demo system (NaCl in water) to validate the
analysis pipeline before extending to more complex electrolyte systems.

## Setup

```bash
conda env create -f environment.yml
conda activate electrolyte-portfolio
```

## Project Structure

- `src/` — analysis toolkit source code
- `examples/` — demo LAMMPS simulations and Jupyter notebooks
- `tests/` — unit tests

## License

MIT