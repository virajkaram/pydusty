#!/usr/bin/env python

# Edits:
#	2015-09-30: added calculation of Rphot,lambda

import argparse
import numpy as np
import matplotlib.pyplot as mpl
import pandas as pd
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
# personal libraries
from mylib import read_header, return_ascii_starting_header

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Return parameter range for a given confidence interval')
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("filename", help="file containing MCMC results")
parser.add_argument("parameter", help="parameter to report", choices=["Rphot","chi2","tstar","tau","td","thick","lstar","rdust","v","ebv","amax","windparam","mej","Rstar"])
parser.add_argument("--confidence", help="confidence interval to report", default=90.0, type=float)
#parser.add_argument("--confidence", help="confidence interval to report", default=68.0, type=float)
parser.add_argument("--burn", help="burn-in length", default=500, type=int)
parser.add_argument("--lamphot", help="wavelength [microns] of photosphere to calculate", type=float, default=4.5)
parser.add_argument("--three_epochs", choices=['1','2','3'])
args = parser.parse_args()

kappa_V = 100.0 # [cm^2 / g]

filename = args.filename
if args.parameter == 'Rphot':   # assumed to mean Rphot_4.5-micron
  # find tau_4.5/tau_V
  file_spectrum = filename.replace('results.dat','best.stb')
  # find the wavelength in spectrum closest to the desired value
  lam, tauT = np.loadtxt(file_spectrum, comments='#', usecols=(0,6), unpack=True)
  lambdaphot_diff = np.abs(lam - args.lamphot)
  indexlamphot = np.argmin(lambdaphot_diff)
  # find the wavelength in spectrum closest to optical
  lambdaV_diff = np.abs(lam - 0.55)
  indexlamV = np.argmin(lambdaV_diff)
  taulamphot_over_taulamV = tauT[indexlamphot]/tauT[indexlamV]
  print 'taulamphot_over_taulamV = %0.4f' % taulamphot_over_taulamV
if args.three_epochs:
  try:
    data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','tau2','tau3','td','td2','td3','thick','lstar','lstar2','lstar3','v','rdust','ebv','amin','amax'], delimiter='\t')
  except:
    try:
      data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','tau2','tau3','td','td2','td3','thick','lstar','lstar2','lstar3','rdust','ebv','amin','amax'], delimiter='\t')
    except:
      print 'Error reading %s' % (filename)
  if args.three_epochs == '2':
    data['tau'] = data['tau2']
    data['td'] = data['td2']
    data['lstar'] = data['lstar2']
  if args.three_epochs == '3':
    data['tau'] = data['tau3']
    data['td'] = data['td3']
    data['lstar'] = data['lstar3']
else:
  try:
    data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','v','ebv','amin','amax','rv','V'], delimiter='\t')
  except:
    try:
      data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','v','ebv','amin','amax'], delimiter='\t')
    except: 
      try:
        data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','v','ebv'], delimiter='\t')
      except:
        print 'Error reading %s' % (filename)
data = data[args.burn::]
if args.parameter == 'Rstar':
  data['Rstar'] = np.log10( np.sqrt( 10**data['lstar']*3.86e33 / (4.*np.pi*5.67e-5*data['tstar']**4) ) )
if args.parameter == 'Rphot':
  data['Rphot'] = np.log10( ( (10**data['rdust'])**-1 - (data['tau']*taulamphot_over_taulamV - 1.0) / (data['tau']*taulamphot_over_taulamV) * ( (10**data['rdust'])**-1 - (10**data['rdust']*data['thick'])**-1 ) )**-1 )
data['windparam'] = 4.0*np.pi*10**(data['rdust'])*data['tau']/kappa_V * (np.pi*1.0e12/1.99e33) # wind density parameter [Msun/yr]
if args.parameter == 'mej':
  data['mej'] = 4*3.14*(10**data['rdust'])**2*data['tau']/100.0 / 2e33
header = return_ascii_starting_header(filename,'#')
tstart = float(read_header(header,'tstart'))
tnow = float(read_header(header,'tnow'))
telap = tnow - tstart
y = np.percentile(data[args.parameter],50)
ymax = np.percentile(data[args.parameter],50+args.confidence/2.0)
ymin = np.percentile(data[args.parameter],50-args.confidence/2.0)
print '  %s = %0.4g < %0.4g < %0.4g \n    %0.4g (+%0.4g -%0.4g)' % (args.parameter, ymin, y, ymax, y, ymax-y, y-ymin)
print '    log %s = %0.4g < %0.4g < %0.4g \n    %0.4g (+%0.4g -%0.4g)' % (args.parameter, np.log10(ymin), np.log10(y), np.log10(ymax), np.log10(y), np.log10(ymax)-np.log10(y), np.log10(y)-np.log10(ymin))
