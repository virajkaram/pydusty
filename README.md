# pydusty
A collection of two python wrappers around the SED modeling code `dusty`

`dusty/` : Scott Adam's python rewrite of Chris Kochanek's MCMC wrapper around `dusty`. This is mainly focused on modeling SEDs of SN2008S-like transients.<br>
`emcee/` : My rewrite of Scott's code, using `emcee` to improve the MCMC and `multiprocessing` to give a significant speed increase. The downside is not all functionalities in Scott's code are implemented yet.
This mainly focuses on modeling stellar SEDs. Constraints on the following parameters can be derived, using either uniform or gaussian priors
1. Stellar temperature
2. Dust temperature
3. Dust optical depth
4. Dust shell thickness
5. Foreground extinction

This can be done for different grain compositions.

The following functionalities are present in Scott's `dusty/' but not in `emcee/`
1. Using observed expansion velocities to constrain dust shell radii
2. Using observed dL/dt as an additional constraint.

It should be straightforward to implement these, I just haven;t encountered data where these are necessary. <br>
Happy to add these in for compelling science cases, or even better, accept contributions here.

## Instructions
For `emcee`, first edit `emcee/dusty_emcee_parallel_molabs.py` with the required prior types, prior ranges and initial values.<br>
An example file with luminosities is given in `data/lums_WISE_J175749.76-075314.9_ebv1.18_nodisterr.dat`.<br>
Run `python emcee/dusty_emcee_parallel_molabs.py`.<br>
For `nprocesses=4`, this will make 4 directories where the total number of iterations are split. These can be visualised using `emcee/plotting_parallel.py` (which will make posterior corner plots, and best-fits).

## Improvements
Agreed that the code structure and ways to provide user-inputs are very suboptimal, but will be improved soon (ish). 
