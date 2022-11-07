#!/usr/bin/env python

import numpy as np
import argparse
from astropy.io import ascii
from sys import exit

def filter_effective_wavelength(filter):
  "return effective wavelength for given filter"
  lambdas = {}
  lambdas['U'] = 0.36   # microns
  lambdas['B'] = 0.44
  lambdas['V'] = 0.55
  lambdas['R'] = 0.658
  lambdas['sdssr'] = 0.612
  lambdas['sdssg'] = 0.464
  lambdas['WFC3uvF555W'] = 0.51874
  lambdas['WFC3uvF814W'] = 0.79011
  lambdas['WFC3irF110W'] = 1.15340
  lambdas['WFC3irF160W'] = 1.53690
  lambdas['lbcbBesselU'] = 0.36
  lambdas['lbcbBesselB'] = 0.44
  lambdas['lbcrBesselV'] = 0.55
  lambdas['lbcrBesselR'] = 0.658
  lambdas['lbcrBesselI'] = 0.806
  lambdas['lbcrBesselY'] = 1.03
  try:
    lambda_eff = lambdas[filter]
  except:
    print 'ERROR: %s-band effective wavelength not found.  Supported filters are:\n%s' % (filter, ''.join("%s\n" % x for x in lambdas.keys()))
    exit(1)
  return lambda_eff

def filter_zeropoint(filter):
  zps = {}
  zps['U'] = 1884.
  zps['B'] = 4620.
  zps['V'] = 3590.
  zps['R'] = 3009.9
  zps['sdssu'] = 3631.
  zps['sdssg'] = 3631.
  zps['sdssr'] = 3631.
  zps['sdssi'] = 3631.
  zps['sdssz'] = 3631.
  try:
    zp = zps[filter]
  except:
    print 'ERROR: %s-band zeropint not found.  Supported filters are:\n%s' % (filter, ''.join("%s\n" % x for x in zps.keys()))
    exit(1)
  return zp

def rl(rv,x):
  "rv is typically 3.1, x is 1/lambda [microns^-1]"
  if x < 1.1:
     a =  0.574*x**1.61
     b = -0.527*x**1.61
  else:
    if x < 3.3:
      y = x-1.82
      a = (1.0+0.17699*y-0.50447*y*y-0.02427*y**3+0.72085*y**4 +
          0.01979*y**5-0.77530*y**6+0.32999*y**7)
      b = (1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-
          0.62251*y**5+5.30260*y**6-2.09002*y**7)
    else:
      if x < 5.9:
        fa = 0.0
        fb = 0.0
      else:
        fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
        fb =  0.21300*(x-5.9)**2 + 0.120700*(x-5.9)**3
      a =  1.752 - 0.316*x - 0.104/((x-4.67)**2+0.341) + fa
      b = -3.090 + 1.825*x + 1.206/((x-4.67)**2+0.263) + fb
  rlval = rv*(a+b/rv)
  return rlval

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("filter", help="name of filter")
parser.add_argument("lum", help="luminosity [Lsun]", type=float)
parser.add_argument("distance", help="Distance [Mpc]", type=float)
parser.add_argument("ebv", help="E(B-V)", type=float)
parser.add_argument("--rv", help="R_V law", type=float, default=3.1)
args = parser.parse_args()

pc = 3.09e18 # [cm]
mpc = 1.e6*pc 
clight = 2.99e10 # cm/s
lsun = 3.83e33 # erg/s
jansky = 1.e-23 # erg
micron = 1.e-4

lam = filter_effective_wavelength(args.filter)
zp = filter_zeropoint(args.filter)

rval = rl(args.rv, 1./lam)
observed_lum = args.lum * 10**(-0.4*rval*args.ebv)
dcm = args.distance*mpc
parser.add_argument("lum", help="luminosity [Lsun]", type=float)
jy = observed_lum * lam*micron*lsun / (jansky*4*np.pi*dcm**2*clight)
mag = -2.5*np.log10(jy/zp)
print 'Observed Luminosity: %0.2f' % (observed_lum)
print 'Apparent mag: %0.2f' % (mag)
