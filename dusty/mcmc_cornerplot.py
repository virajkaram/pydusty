#!/usr/bin/env python
# This program was copied from plot2d_hist.py on 9/15/2014
# Changed to make multi-panel figure for paper

import numpy as np
import corner
import argparse
from sys import exit
from astropy.io import ascii

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Generate corner plot")
parser.add_argument("--fields", help="parameters to plot", nargs='+', default=['lstar','tstar','tau','td','r1','thick'])
parser.add_argument("--input", help="input data file", default='results.dat')
parser.add_argument("--output", help="output plot name", default='mcmc.pdf')
parser.add_argument("--burn_in_length", help="number of lines to skip from beginning of input file", type=int, default=500)
#parser.add_argument("--range")
args = parser.parse_args()

try:
  data = ascii.read(args.input, names=['chi2','tstar','tau','td','thick','lstar','r1','v','ebv','amin','amax','rv'])
except:
  try:
    data = ascii.read(args.input, names=['chi2','tstar','tau','td','thick','lstar','r1','v','ebv','amin','amax'])
  except:
    print 'error reading in %s' % (args.input)
    exit(-1)

kappa_V = 100.0 # [cm^2 / g]
v_w = 10.0 # [km/s]

if 'Rstar' in args.fields:
  data['Rstar'] = np.log10( np.sqrt( 10**data['lstar']*3.86e33 / (4.*np.pi*5.67e-5*data['tstar']**4) ) )
if 'Rphot' in args.fields:
  data['Rphot'] = np.log10( ( (10**data['r1'])**-1 - (data['tau']*taulamphot_over_taulamV - 1.0) / (data['tau']*taulamphot_over_taulamV) * ( (10**data['r1'])**-1 - (10**data['r1']*data['thick'])**-1 ) )**-1 )
if 'windparam' in args.fields:
  data['windparam'] = 4.0*np.pi*10**(data['r1'])*data['tau']/kappa_V * (np.pi*1.0e12/1.99e33) # wind density parameter [Msun/yr]
if 'mdot' in args.fields:
  data['mdot'] = 4.0*np.pi*10**(data['r1'])*data['tau']/kappa_V * (np.pi*1.0e12/1.99e33)*v_w # wind density parameter [Msun/yr]
if 'mej' in args.fields:
  data['mej'] = 4*3.14*(10**data['r1'])**2*data['tau']/100.0 / 2e33

# identify best-fit parameters
argmin = data['chi2'].argmin()
best = [data[argmin][field] for field in args.fields]

# trim off burn-in
data = data[args.burn_in_length:]

# select data fields to plot and transpose to format expected by corner module
samples = np.transpose([data[field] for field in args.fields])

# convert fields to pretty names
names = []
for field in args.fields:
  if field == 'lstar': names.append('$\mathrm{log}$ $L/L_{\odot}$')
  elif field == 'tstar': names.append('$T_{*}$ [K]')
  elif field == 'tau': names.append(r'$\tau_{\mathrm{V}}$')
  elif field == 'td': names.append('$T_{\mathrm{dust}}$ [K]')
  elif field == 'r1': names.append('log( Inner Dust Radius / cm )')
  elif field == 'ebv': names.append('E(B-V)')
  elif field == 'thick': names.append('$R_{\mathrm{out}}/R_{\mathrm{in}}$')
  elif field == 'mej': names.append('$M_{\mathrm{ej}}/M_{\odot}$')
  elif field == 'windparam': names.append('$\dot{M}/v_{\mathrm{w}}$')
  elif field == 'mdot': names.append('$\dot{M}/M_{\odot}\mathrm{yr}^{-1}$')
#fig = corner.corner(samples, labels=names, truths=best)
fig = corner.corner(samples, labels=names, truths=best, range=[(0,8),(400,2000),(14,16),(0,0.0001),(0,1e-5)])
fig.savefig(args.output)
