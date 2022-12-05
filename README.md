# pydusty
A generic python wrapper for performing MCMC analysis with the SED modeling code `dusty`.

This code is based on Scott Adam's rewrite of Chris Kochanek's MCMC wrapper around `dusty`, and my subsequent rewrite of Scott's wrapper using `emcee'. The original codes are available under `pydusty/legacy`.

## Installation
1. Clone the repo using `git clone`
2. `cd pydusty`
3. `pip install -e .`

The dustyv2 source code files are under `data/dusty_files`. You may need to recreate the executable file with `gfortran dusty.f -o dusty`. Please see the dusty github repository (esp. https://github.com/ivezic/dusty/issues/7 issue) for more details.

## Usage
An example usage file is given under `examples/run_dusty_mcmc.py (--h)`. Detailed documentation of the code will follow soon.

The following functionalities present in Scott's `dusty/' have not yet been implemented in this version
1. Using observed expansion velocities to constrain dust shell radii
2. Using observed dL/dt as an additional constraint.

It should be straightforward to implement these, I just haven;t encountered data where these are necessary. <br>
Happy to add these in for compelling science cases, or even better, accept contributions here.
