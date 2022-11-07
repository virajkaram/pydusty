#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse

# getplot.pl best
def getplot(fileroot):
  # find luminosity
  llstar, junk = np.loadtxt(fileroot+'.dat', usecols=(5,6), skiprows=1, unpack=True)
  lstar = 10**llstar
  return lstar

def dospec1(fileroot, color='black'):
  lstar = getplot(fileroot)
  # normlum = emerging spectral flux
  # xatt    = fractional contribution of attenuated input radiation
  # xds     = fractional contribution of scattered radiation
  # xde     = fractional contribution of dust emission
  # finp    = input spectral shape
  # taut    = overal optical depth at lambda
  # albedo = albedo at wavelength lambda
  lam, normlum, xatt, xds, xde, finp, taut, albedo = np.loadtxt(fileroot+'.stb', usecols=(0,1,2,3,4,5,6,7), unpack=True)
  lum = lstar*normlum
  llam = np.log10(lam)
  llum = np.log10(lum+1.e-20)
  plt.plot(10**llam, 10**llum,'-', color=color)
  #return llam, llum

def plotdata(filename, color='black'):
  olam, olum, oerr = np.loadtxt(filename, usecols=(0,1,2), unpack=True)
  measured = (olum>0)
  llamuse = np.log10(olam[measured])
  llumuse = np.log10(olum[measured])
  lerruse = oerr[measured]/olum[measured]/np.log(10) # double check that this should be a natural log...
  plt.errorbar(10**llamuse, 10**llumuse, yerr=oerr[measured], fmt='.', color=color)
  limit = (olum<0)
  llamuse = np.log10(olam[limit])
  llumuse = np.log10(oerr[limit])
  plt.errorbar(10**llamuse, 10**llumuse, yerr=0.3*10**(llumuse), uplims=True, fmt='.', color=color)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Plot best-fit SED")
parser.add_argument("--magfile", help="name of file containing observed flux constraints", default="premags.dat")
parser.add_argument("--minlam", help="minimum wavelength to plot", type=float, default=10**-0.5)
parser.add_argument("--maxlam", help="maximum wavelength to plot", type=float, default=10**1)
parser.add_argument("--minlum", help="minimum luminosity to plot", type=float, default=10**1)
parser.add_argument("--maxlum", help="maximum luminosity to plot", type=float, default=10**5.5)
parser.add_argument("--output", help="filename of output plot", default="tmp.pdf")
args = parser.parse_args()

plt.figure()
#plt.rc('font', family='serif')
#plt.rc('text', usetex=True)
#dospec1('best')
plotdata(args.magfile)
plt.xlim(args.minlam,args.maxlam)
plt.ylim(args.minlum,args.maxlum)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\lambda L_{\lambda}$ [$L_{\odot}$]')
plt.xlabel(r'Wavelength [$\mu\mathrm{m}$]')
plt.savefig(args.output)
