#!/usr/bin/env python
# This program was copied from plot2d_hist.py on 9/15/2014
# Changed to make multi-panel figure for paper

import numpy as np
import matplotlib.pyplot as mpl
import argparse
from astropy.io import ascii
from sys import argv, exit
from matplotlib.ticker import AutoMinorLocator

def getdata(mcmc_datafile):
  data = ascii.read(mcmc_datafile)
  chi2 = data['col1'][burn_in_length:args.max_length]
  tstar = data['col2'][burn_in_length:args.max_length]
  tau = data['col3'][burn_in_length:args.max_length]
  tdust = data['col4'][burn_in_length:args.max_length]
  thick = data['col5'][burn_in_length:args.max_length]
  sluml = data['col6'][burn_in_length:args.max_length]
  r1 = data['col7'][burn_in_length:args.max_length]
  return chi2,tstar,tau,tdust,thick,sluml,r1

parser = argparse.ArgumentParser()
parser.add_argument("x", help='parameter for x-axis', choices=['chi2','tstar','tau','tdust','thick','sluml','r1'])
parser.add_argument("y", help='parameter for y-axis',choices=['chi2','tstar','tau','tdust','thick','sluml','r1'])
parser.add_argument("z", help='parameter for color-axis',choices=['chi2','tstar','tau','tdust','thick','sluml','r1'])
parser.add_argument("--input", help="input data file", default='results.dat')
parser.add_argument("--output", help="output plot name", default='mcmc.pdf')
parser.add_argument("--burn_in_length", help="number of lines to skip from beginning of input file", type=int, default=500)
parser.add_argument("--max_length", help="total number of lines to include", type=int, default=100000)
#parser.add_argument("--progenitor_file", help="name of progenitor data file", default="run1/progenitor3_tau0.dat")
parser.add_argument("--cbar_min", help="Minimum value for colorbar", type=float, default='-10')
parser.add_argument("--cbar_max", help="Maximum value for colorbar", type=float, default='-10')
args = parser.parse_args()

#progenitor_file = args.progenitor_file
burn_in_length = args.burn_in_length

#progenitor_data = ascii.read(progenitor_file)

fig = mpl.figure(figsize=(7,5))
#fig, p1 = mpl.subplots(1)
fig.subplots_adjust(wspace=0, hspace=0) # no space between subplots (no pad)
mymap = mpl.get_cmap('jet')

bigAxes = fig.add_axes([0.1,0.1,0.8,0.8])
bigAxes.set_frame_on(False)
#bigAxes.spines['top'].set_color('none')
#bigAxes.spines['bottom'].set_color('none')
#bigAxes.spines['left'].set_color('none')
#bigAxes.spines['right'].set_color('none')
bigAxes.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
#bigAxes.axes.get_yaxis().set_visible(False)
#bigAxes.axes.get_xaxis().set_visible(False)
#bigAxes.set_xlabel(r'T$_{*}$ [K]')
#bigAxes.set_ylabel(r'log(L/L$_{\odot}$')
bigAxes.set_xlabel(args.x)
bigAxes.set_ylabel(args.y)

# Top left panel
data = {}
data['chi2'],data['tstar'],data['tau'],data['tdust'],data['thick'],data['sluml'],data['r1'] = getdata(args.input)
ax1 = fig.add_axes([0.13,0.11,0.7,0.8])
s = ax1.scatter(data[args.x],data[args.y],s=3,c=data[args.z],edgecolors='None',cmap=mymap, vmin=data[args.z].min(), vmax=data[args.z].max())
#ax1.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
#ax1.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
ax1.set_xlim(data[args.x].min(),data[args.x].max())
ax1.set_ylim(data[args.y].min(),data[args.y].max())
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
#mpl.setp(ax1.get_yticklabels(), visible=False)
#mpl.setp(ax1.get_xticklabels(), visible=False)

#ax1.legend(prop={'size':12},numpoints=1,loc=4,markerscale=4)

cax = fig.add_axes([0.83,0.11,0.02,0.8])
cbar = fig.colorbar(s,cax)
cbar.set_label(args.z)

mpl.savefig(args.output)
