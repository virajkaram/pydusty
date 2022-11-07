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

# Shell models
def plot_shell_models():
  # 2x4 panels
  # left - no dL/dt prior
  # right - with dL/dt prior
  # row 1 - V + I
  # row 2 - V only
  # row 3 - V + 4.5-micron
  # row 4 - only upper limits
  fig = mpl.figure(figsize=(7,10))
  fig.subplots_adjust(wspace=0, hspace=0) # no space between subplots (no pad)
  mymap = mpl.get_cmap('jet')
  #mymap = mpl.get_cmap('cubehelix')
  bigAxes = fig.add_axes([0.12,0.1,0.8,0.8])
  bigAxes.set_frame_on(False)
  bigAxes.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
  bigAxes.set_xlabel(r'T$_{*}$ [K]')
  bigAxes.set_ylabel(r'log L/L$_{\odot}$')

  progenitor_data = ascii.read('run1/progenitor3_tau0.dat')
  #colorcycle = mpl.rcParams['axes.color_cycle']
  #colors = [colorcycle[1],colorcycle[2],colorcycle[3],colorcycle[4]]

  # row 1, left panel - V + I, no dL/dt
  ax1 = fig.add_axes([0.13,0.71,0.35,0.2])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_shell_f814f555_mcmc/results.dat')
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated7_aMRN_2014_silicate_shell_f814f555/results.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated7_aMRN_2014_silicate_shell_f814f555/results.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax1.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax1.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax1.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax1.text(7000,6.0,'V+I',size=15)
  ax1.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax1.set_xlim(3500,29999.99)
  ax1.set_ylim(3.0,6.5)
  ax1.set_xticks([5000,15000,25000])
  ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax1.get_xticklabels(), visible=False)

  # row 1, right panel - V + I, with dL/dt
  ax2 = fig.add_axes([0.48,0.71,0.35,0.2])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_shell_f814f555_dL_mcmc/results.dat')
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated7_aMRN_2014_silicate_shell_dL_f814f555/results.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated7_aMRN_2014_silicate_shell_dL_f814f555/results.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax2.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax2.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax2.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax2.text(7000,6.0,'V+I+LC',size=15)
  ax2.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax2.set_xlim(3500,30000)
  ax2.set_ylim(3.0,6.5)
  ax2.set_xticks([5000,15000,25000])
  ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax2.get_xticklabels(), visible=False)
  mpl.setp(ax2.get_yticklabels(), visible=False)

  # row 2, left panel - V only, no dL/dt
  ax3 = fig.add_axes([0.13,0.51,0.35,0.2])
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated7_aMRN_2014_silicate_shell_f814only/results.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated7_aMRN_2014_silicate_shell_f814only/results.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax3.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax3.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax3.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax3.text(7000,6.0,'I',size=15)
  ax3.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax3.set_xlim(3500,29999.99)
  ax3.set_ylim(3.0,6.4999)
  ax3.set_xticks([5000,15000,25000])
  ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax3.get_xticklabels(), visible=False)

  # row 2, right panel - V only, with dL/dt
  ax4 = fig.add_axes([0.48,0.51,0.35,0.2])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_shell_f814only_dL_mcmc/results.dat')
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated7_aMRN_2014_silicate_shell_dL_f814only/results.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated7_aMRN_2014_silicate_shell_dL_f814only/results.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax4.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax4.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax4.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax4.text(7000,6.0,'I+LC',size=15)
  ax4.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax4.set_xlim(3500,30000)
  ax4.set_ylim(3.0,6.5)
  ax4.set_xticks([5000,15000,25000])
  ax4.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax4.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax4.get_xticklabels(), visible=False)
  mpl.setp(ax4.get_yticklabels(), visible=False)

  # row 3, left panel - V + 4.5-micron, no dL/dt
  ax5 = fig.add_axes([0.13,0.31,0.35,0.2])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_shell_f814ch2_mcmc/results.dat')
  #if verbose: print '%s -- tau_max = %0.1f' % ('2014_silicate_shell_f814ch2_mcmc/results.dat', np.max(tau))
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated9_aMRN_2014_silicate_shell_f814ch2/results_manual.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated9_aMRN_2014_silicate_shell_f814ch2/results_manual.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax5.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax5.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax5.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2) 
  ax5.text(7000,6.0,r'I+4.5$\mu\mathrm{m}$',size=15)
  ax5.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax5.set_xlim(3500,29999.99)
  ax5.set_ylim(3.0,6.4999)
  ax5.set_xticks([5000,15000,25000])
  ax5.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax5.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax5.get_xticklabels(), visible=False)

  # row 3, right panel - V + 4.5-micron, with dL/dt
  ax6 = fig.add_axes([0.48,0.31,0.35,0.2])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_shell_f814ch2_dL_mcmc/results.dat')
  try:
    chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated9_aMRN_2014_silicate_shell_dL_f814ch2/results_manual.dat')
    if verbose: print '%s -- tau_max = %0.1f' % ('updated9_aMRN_2014_silicate_shell_dL_f814ch2/results_manual.dat', np.max(tau))
    if plot_tau_eff:
      tau_coeff = float(read_header(header, 'tau_eff_coeff'))
      tau = tau*tau_coeff
    s = ax6.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
    ax6.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
    ax6.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
    ax6.text(7000,6.0,r'I+4.5$\mu\mathrm{m}$+LC',size=15)
    ax6.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  except:
    print 'error plotting updated9_aMRN_2014_silicate_shell_dL_f814ch2/results_manual.dat -- skipping'
  ax6.set_xlim(3500,30000)
  ax6.set_ylim(3.0,6.5)
  ax6.set_xticks([5000,15000,25000])
  ax6.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax6.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax6.get_xticklabels(), visible=False)
  mpl.setp(ax6.get_yticklabels(), visible=False)

  # row 4, left panel - only upper limits, no dL/dt
  # row 4, right panel - limits only, with dL/dt
  ax7 = fig.add_axes([0.13,0.11,0.35,0.2])
  ax8 = fig.add_axes([0.48,0.11,0.35,0.2])
  ax7.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax8.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax7.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax8.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  data = {}
  minchilum = 1e10
  minchitot = 1e10
  for tau in taus:
    data['tau'+str(tau)] = ascii.read('grid'+str(ncol)+'_aMRN_2014_silicate_shell_dL_limits/limits_grid_tau'+str(tau)+'.dat')
    data['tau'+str(tau)]['totalchi'] = data['tau'+str(tau)]['chi2'] + data['tau'+str(tau)]['chi_dLdt']    
    fileminchilum = data['tau'+str(tau)]['chi2'].min()
    fileminchitot = data['tau'+str(tau)]['totalchi'].min()
    minchilum = min(minchilum,fileminchilum)
    minchitot = min(minchitot,fileminchitot)
  for tau in taus:
    x = np.reshape(data['tau'+str(tau)]['tstar'], (-1,ncol))
    y = np.reshape(data['tau'+str(tau)]['lum'], (-1,ncol))
    zlum = np.reshape(data['tau'+str(tau)]['chi2']-minchilum, (-1,ncol))
    ztot = np.reshape(data['tau'+str(tau)]['totalchi']-minchitot, (-1,ncol))
    linecolor = mymap(normcolor(tau))
    ax7.contourf(x,y,zlum, levels=sigmas, colors=(linecolor[0:3],(0,0,0)), alpha=myalpha)
    ax8.contourf(x,y,ztot, levels=sigmas, colors=(linecolor[0:3],(0,0,0)), alpha=myalpha)
  ax7.text(7000,6.0,'Limits Only',size=15)
  ax7.set_xlim(3500,29999.99)
  ax7.set_ylim(3.0,6.4999)
  ax7.set_xticks([5000,15000,25000])
  ax7.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax7.yaxis.set_minor_locator(AutoMinorLocator(5))
  ax8.text(7000,6.0,'Limits+LC',size=15)
  ax8.set_xlim(3500,29999.99)
  ax8.set_ylim(3.0,6.4999)
  ax8.set_xticks([5000,15000,25000])
  ax8.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax8.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax8.get_yticklabels(), visible=False)

  cax = fig.add_axes([0.83,0.11,0.02,0.8])
  cbar = fig.colorbar(s,cax)
  if plot_tau_eff:
    cbar.set_label(r'$\tau_{V,\mathrm{eff}}$')
  else:
    cbar.set_label(r'$\tau_{V,\mathrm{tot}}$')
  mpl.savefig('fig4.pdf',bbox_inches='tight')
  mpl.savefig('fig4.png',dpi=100,bbox_inches='tight')

# Shell models
def plot_wind_models():
  # 2x2 panels
  # row 1 left - V+I
  # row 1 right - V
  # row 2 left - V+4.5-micron
  # row 2 right - limits only
  # fake figure to get solid colorbar
  fig = mpl.figure(figsize=(7,6))
  fig.subplots_adjust(wspace=0, hspace=0) # no space between subplots (no pad)
  mymap = mpl.get_cmap('jet')
  #mymap = mpl.get_cmap('cubehelix')
  bigAxes = fig.add_axes([0.12,0.1,0.8,0.8])
  bigAxes.set_frame_on(False)
  bigAxes.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
  bigAxes.set_xlabel(r'T$_{*}$ [K]')
  bigAxes.set_ylabel(r'log L/L$_{\odot}$')

  progenitor_data = ascii.read('run1/progenitor3_tau0.dat')
  #colorcycle = mpl.rcParams['axes.color_cycle']
  #colors = [colorcycle[1],colorcycle[2],colorcycle[3],colorcycle[4]]

  # row 1, left panel - V + I
  ax1 = fig.add_axes([0.13,0.51,0.35,0.4])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_wind_f814f555_mcmc/results.dat')
  #if verbose: print '%s -- tau_max = %0.1f' % ('2014_silicate_wind_f814f555_mcmc/results.dat', np.max(tau))
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated7_aMRN_2014_silicate_wind_f814f555/results.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated7_aMRN_2014_silicate_wind_f814f555/results.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax1.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax1.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax1.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax1.text(7000,6.0,'V+I',size=15)
  ax1.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax1.set_xlim(3500,29999.99)
  ax1.set_ylim(3.0,6.5)
  ax1.set_xticks([5000,15000,25000])
  ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax1.get_xticklabels(), visible=False)

  # row 1, right panel - V
  ax2 = fig.add_axes([0.48,0.51,0.35,0.4])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_wind_f814only_mcmc/results.dat')
  #if verbose: print '%s -- tau_max = %0.1f' % ('2014_silicate_wind_f814only_mcmc/results.dat', np.max(tau))
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated7_aMRN_2014_silicate_wind_f814only/results.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated7_aMRN_2014_silicate_wind_f814only/results.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax2.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax2.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax2.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax2.text(7000,6.0,'I',size=15)
  ax2.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax2.set_xlim(3500,30000)
  ax2.set_ylim(3.0,6.5)
  ax2.set_xticks([5000,15000,25000])
  ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax2.get_xticklabels(), visible=False)
  mpl.setp(ax2.get_yticklabels(), visible=False)

  # row 2, left panel - V+4.5-micron
  ax3 = fig.add_axes([0.13,0.11,0.35,0.4])
  #chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('2014_silicate_wind_f814ch2_mcmc/results.dat')
  #if verbose: print '%s -- tau_max = %0.1f' % ('2014_silicate_wind_f814ch2_mcmc/results.dat', np.max(tau))
  chi2,tstar,tau,tdust,thick,sluml,r1,header = getmcmcdata('updated8_aMRN_2014_silicate_wind_f814ch2/results.dat')
  if verbose: print '%s -- tau_max = %0.1f' % ('updated8_aMRN_2014_silicate_wind_f814ch2/results.dat', np.max(tau))
  if plot_tau_eff:
    tau_coeff = float(read_header(header, 'tau_eff_coeff'))
    tau = tau*tau_coeff
  s = ax3.scatter(tstar,sluml,s=pointsize,c=tau,edgecolors='None',cmap=mymap, vmin=cbar_min, vmax=cbar_max)
  ax3.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax3.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  ax3.text(7000,6.0,r'I+4.5$\mu\mathrm{m}$',size=15)
  ax3.text(18000,3.3,r'$\chi^2_{\mathrm{min}}=%0.1f$' % (np.min(chi2)), size=14)
  ax3.set_xlim(3500,29999.99)
  ax3.set_ylim(3.0,6.4999)
  ax3.set_xticks([5000,15000,25000])
  ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax3.yaxis.set_minor_locator(AutoMinorLocator(5))

  # row 2, right panel - limits only
  ax4 = fig.add_axes([0.48,0.11,0.35,0.4])
  ax4.plot(progenitor_data['tstar'],progenitor_data['lum'],'-', color='black', label='Prog')
  ax4.fill_between(progenitor_data['tstar'],progenitor_data['lum_max'],progenitor_data['lum_min'], color='gray', alpha=0.2)
  data = {}
  minchilum = 1e10
  for tau in taus:
    data['tau'+str(tau)] = ascii.read('grid'+str(ncol)+'_aMRN_2014_silicate_wind_limits/limits_grid_tau'+str(tau)+'.dat')
    fileminchilum = data['tau'+str(tau)]['chi2'].min()
    minchilum = min(minchilum,fileminchilum)
  for tau in taus:
    x = np.reshape(data['tau'+str(tau)]['tstar'], (-1,ncol))
    y = np.reshape(data['tau'+str(tau)]['lum'], (-1,ncol))
    zlum = np.reshape(data['tau'+str(tau)]['chi2']-minchilum, (-1,ncol))
    linecolor = mymap(normcolor(tau))
    ax4.contourf(x,y,zlum, levels=sigmas, colors=(linecolor[0:3],(0,0,0)), alpha=myalpha)
  ax4.text(7000,6.0,'Limits Only',size=15)
  ax4.set_xlim(3500,29999.99)
  ax4.set_ylim(3.0,6.4999)
  ax4.set_xticks([5000,15000,25000])
  ax4.xaxis.set_minor_locator(AutoMinorLocator(5))
  ax4.yaxis.set_minor_locator(AutoMinorLocator(5))
  mpl.setp(ax4.get_yticklabels(), visible=False)

  # color bar
  cax = fig.add_axes([0.83,0.11,0.02,0.8])
  cbar = fig.colorbar(s,cax)
  if plot_tau_eff:
    cbar.set_label(r'$\tau_{V,\mathrm{eff}}$')
  else:
    cbar.set_label(r'$\tau_{V,\mathrm{tot}}$')
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
ncol = 50	# number of luminosity columns in limits plots
taus = [0,1,3,10]	# tau lines to draw in limits plots

pointsize = 1
myalpha = 0.8
cbar_min = 0
cbar_max = 17
normcolor = matplotlib.colors.normalize(vmin=cbar_min,vmax=cbar_max)
burn_in_length = 500
#burn_in_length = 75
max_length = 10000 + burn_in_length

plot_shell_models()
plot_wind_models()

