#!/usr/bin/env python
# This program was copied from mcmc_plot_paper.py on 11/14/2014
# Make multi-panel figures showing the MCMC shell and wind models

import numpy as np
import matplotlib.pyplot as mpl
import matplotlib
from matplotlib.ticker import AutoMinorLocator
# non-standard libraries:
from datetime import datetime
from astropy.io import ascii
# personal libraries:
from mylib import read_header, return_ascii_starting_header

def getmcmcdata(mcmc_datafile):
  data = ascii.read(mcmc_datafile)
  rdr = ascii.get_reader()
  header = return_ascii_starting_header(mcmc_datafile,'#')
  chi2 = data['col1'][burn_in_length:max_length]
  tstar = data['col2'][burn_in_length:max_length]
  tau = data['col3'][burn_in_length:max_length]
  tdust = data['col4'][burn_in_length:max_length]
  thick = data['col5'][burn_in_length:max_length]
  sluml = data['col6'][burn_in_length:max_length]
  r1 = data['col7'][burn_in_length:max_length]
  return chi2,tstar,tau,tdust,thick,sluml,r1, header

def plot_wind_models():
  fig = mpl.figure(figsize=(7,6))
  fig.subplots_adjust(wspace=0, hspace=0) # no space between subplots (no pad)
  mymap = mpl.get_cmap('jet')
  #mymap = mpl.get_cmap('cubehelix')
  bigAxes = fig.add_axes([0.12,0.1,0.8,0.8])
  bigAxes.set_frame_on(False)
  bigAxes.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
  bigAxes.set_xlabel(r'T$_{*}$ [K]')
  bigAxes.set_ylabel(r'log L/L$_{\odot}$')

  #progenitor_data = ascii.read('run1/progenitor3_tau0.dat')
  #colorcycle = mpl.rcParams['axes.color_cycle']
  #colors = [colorcycle[1],colorcycle[2],colorcycle[3],colorcycle[4]]

  ax1 = fig.add_axes([0.13,0.11,0.8,0.8])
  #ax1.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  #ax1.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  data = {}
  minchilum = 1e10
  for tau in taus:
    data['tau'+str(tau)] = ascii.read('limits_grid_tau'+str(tau)+'.dat')
    fileminchilum = data['tau'+str(tau)]['chi2'].min()
    minchilum = min(minchilum,fileminchilum)
  for tau in taus:
    x = np.reshape(data['tau'+str(tau)]['tstar'], (-1,ncol))
    y = np.reshape(data['tau'+str(tau)]['lum'], (-1,ncol))
    zlum = np.reshape(data['tau'+str(tau)]['chi2']-minchilum, (-1,ncol))
    linecolor = mymap(normcolor(tau))
    ax1.contourf(x,y,zlum, levels=sigmas, colors=(linecolor[0:3],(0,0,0)), alpha=myalpha)
  #ax1.text(7000,6.0,'Limits Only',size=15)
  ax1.set_xlim(3500,30000)
  ax1.set_ylim(4.0,6.5)
  ax1.set_xticks([5000,15000,25000])
  ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
  #mpl.setp(ax1.get_yticklabels(), visible=False)

  # color bar
  #cax = fig.add_axes([0.83,0.11,0.02,0.8])
  #cbar = fig.colorbar(s,cax)
  #if plot_tau_eff:
  #  cbar.set_label(r'$\tau_{V,\mathrm{eff}}$')
  #else:
  #  cbar.set_label(r'$\tau_{V,\mathrm{tot}}$')
  mpl.savefig('fig6.pdf',bbox_inches='tight')
  mpl.savefig('fig6.png',dpi=100,bbox_inches='tight')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# set output verbosity
verbose = 1

# use tau_eff (rather than tau_tot) for color-coding points
plot_tau_eff = 0

sigmas = np.array([6.25, 21.1])	# delta chi^2 contours to draw 
#	6.25 and 21.1 correspond to 90 and 99.99% confidence intervals (respectively) for 3 parameters (L,T,tau)
ncol = 20	# number of luminosity columns in limits plots
taus = [0,1]	# tau lines to draw in limits plots

pointsize = 1
myalpha = 0.8
cbar_min = 0
cbar_max = 17
normcolor = matplotlib.colors.Normalize(vmin=cbar_min,vmax=cbar_max)
burn_in_length = 500
#burn_in_length = 75
max_length = 10000 + burn_in_length

#plot_shell_models()
plot_wind_models()

