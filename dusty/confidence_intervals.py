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
#parser.add_argument("parameters", help="parameters to report (Rphot,chi2,tstar,tau,td,thick,lstar,rdust,v,ebv,amax,windparam,mej)")
#, choices=["Rphot","chi2","tstar","tau","td","thick","lstar","rdust","v","ebv","amax","windparam","mej"])
parser.add_argument("--confidence", help="confidence interval to report", default=90.0, type=float)
#parser.add_argument("--confidence", help="confidence interval to report", default=68.0, type=float)
parser.add_argument("--latex", help="output with latex formatting", action="store_true")
parser.add_argument("--burn", help="burn-in length", default=500, type=int)
parser.add_argument("--lamphot", help="wavelength [microns] of photosphere to calculate", type=float, default=4.5)
args = parser.parse_args()

kappa_V = 100.0 # [cm^2 / g]

filename = args.filename
#parameters = args.parameters.split(',')
parameters = ['lstar','tstar','tau','td','thick','v','mej','windparam','ebv','chi2']
#  if param == 'Rphot':   # assumed to mean Rphot_4.5-micron
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
try:
  data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','v','ebv','amin','amax','Rv'], delimiter='\t')
except:
  try:
    data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','v','ebv','amin','amax'], delimiter='\t')
  except:
    try:
      data = pd.read_csv(filename, comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','v','ebv'], delimiter='\t')
    except:
      print 'Error reading %s' % (filename)
data = data[args.burn::]
#if param == 'Rphot':
data['Rphot'] = np.log10( ( (10**data['rdust'])**-1 - (data['tau']*taulamphot_over_taulamV - 1.0) / (data['tau']*taulamphot_over_taulamV) * ( (10**data['rdust'])**-1 - (10**data['rdust']*data['thick'])**-1 ) )**-1 )
data['windparam'] = np.log10(4.0*np.pi*10**(data['rdust'])*data['tau']/kappa_V * (np.pi*1.0e12/1.99e33)) # wind density parameter log [Msun/yr]
#if args.parameter == 'mej':
data['mej'] = np.log10(4*3.14*(10**data['rdust'])**2*data['tau']/100.0 / 2e33) # log M/Msun
  
header = return_ascii_starting_header(filename,'#')
tstart = float(read_header(header,'tstart'))
tnow = float(read_header(header,'tnow'))
telap = tnow - tstart

results = {}
for param in parameters:
  y = np.percentile(data[param],50)
  ymax = np.percentile(data[param],50+args.confidence/2.0)
  ymin = np.percentile(data[param],50-args.confidence/2.0)
  print '  %s = %0.4g < %0.4g < %0.4g \n    %0.4g (+%0.4g -%0.4g)' % (param, ymin, y, ymax, y, ymax-y, y-ymin)
  print '    log %s = %0.4g < %0.4g < %0.4g \n    %0.4g (+%0.4g -%0.4g)' % (param, np.log10(ymin), np.log10(y), np.log10(ymax), np.log10(y), np.log10(ymax)-np.log10(y), np.log10(y)-np.log10(ymin))
  results[param] = {}
  results[param]['perr'] = ymax-y
  results[param]['merr'] = y-ymin
  results[param]['med'] = y

bestdata = pd.read_csv(filename.replace('results.dat','best.dat'), comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','ebv','amin','amax','Rv'], delimiter='\t', skiprows=1)
#bestdata = pd.read_csv(filename.replace('results.dat','best.dat'), comment='#', names=['chi2','tstar','tau','td','thick','lstar','rdust','ebv','amin','amax'], delimiter='\t')

#for i in xrange(len(parameters)):
#parameters = ['lstar','tstar','tau','td','thick','ebv']
parameters = ['lstar','tstar','tau','td','thick','v','mej','ebv']
print '%', args.filename
if len(data)<10000:
  print '(unfinished)'
else:
  print ''
for i in xrange(len(parameters)):
  if parameters[i] != 'v' and parameters[i] != 'mej' and parameters[i] != 'windparam':
    if bestdata[parameters[i]].values[0] > results[parameters[i]]['med']+results[parameters[i]]['perr'] or (bestdata[parameters[i]].values[0] < results[parameters[i]]['med']-results[parameters[i]]['merr']):
      #print '%% WARNING: best %s' % (parameters[i])
      print '%%\t WARNING: best %s of %0.2f is outside of confidence interval of %0.2f < %0.2f < %0.2f' % (parameters[i], bestdata[parameters[i]], results[parameters[i]]['med']-results[parameters[i]]['merr'], results[parameters[i]]['med'], results[parameters[i]]['med']+results[parameters[i]]['perr'])
  if results[parameters[i]]['perr'] == 0:
    #if parameters[i] == 'lstar':
    print '%0.2f (fixed) &\t' % (results[parameters[i]]['med'])
    #else:
    #  print '%d (fixed) &\t' % (results[parameters[i]]['med'])
  elif parameters[i] == 'lstar':
    print '$%0.2f^{+%0.2g}_{-%0.2f}$ &\t' % (results['lstar']['med'], results['lstar']['perr'], results['lstar']['merr'])
  elif parameters[i] == 'tstar':
    print '$%d^{+%d}_{-%d}$ &\t' % (results['tstar']['med'], results['tstar']['perr'], results['tstar']['merr'])
  elif parameters[i] == 'tau':
    print '$%0.1f^{+%0.1f}_{-%0.1f}$ &\t' % (results['tau']['med'], results['tau']['perr'], results['tau']['merr'])
  elif parameters[i] == 'td':
    print '$%d^{+%d}_{-%d}$ &\t' % (results['td']['med'], results['td']['perr'], results['td']['merr'])
  elif parameters[i] == 'thick':
    print '$%0.1f^{+%0.1f}_{-%0.1f}$ &\t' % (results['thick']['med'], results['thick']['perr'], results['thick']['merr'])
  elif parameters[i] == 'v':
    print '$%0.1f^{+%0.1f}_{-%0.1f}$ &\t' % (results['v']['med'], results['v']['perr'], results['v']['merr'])
  elif parameters[i] == 'mej':
    if results['thick']['med'] != 2:
      print '$%0.1f^{+%0.1f}_{-%0.1f}$ &\t' % (results['windparam']['med'], results['windparam']['perr'], results['windparam']['merr'])
    else:
      print '$%0.1f^{+%0.1f}_{-%0.1f}$ &\t' % (results['mej']['med'], results['mej']['perr'], results['mej']['merr'])
  elif parameters[i] == 'ebv':
    print '$%0.2f^{+%0.2f}_{-%0.2f}$ &\t' % (results['ebv']['med'], results['ebv']['perr'], results['ebv']['merr'])
  #elif parameters[i] == 'chi2':
    #  print '$%0.1f^{+%0.1f}_{-%0.1f}$ \t\\\\' % (results['chi2']['med'], results['chi2']['perr'], results['chi2']['merr'])
  else:
    print 'unrecognized parameter: %s' % (parameters[i])
print '$%0.1f$ \t\\\\' % (bestdata['chi2'].values[0])

