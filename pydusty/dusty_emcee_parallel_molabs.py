#!/usr/bin/env python
# coding: utf-8
# Author : Viraj Karambelkar
# My rewrite of Scott's rewrite of Kochanek's MCMC wrapper
# This code uses emcee instead of Scott's manual mcmc
# This code can run in parallel, significantly speeding up the process
# Not all functionalities of Scott's code are implemented yet, the parameters that can be varied are described
# in my 2021 R-Cor-Bor paper (Karambelkar et al. 2021)
# At some point, this will be refactored so that it can be better used by everyone, stay tuned..

import numpy as np
import emcee
from astropy.io import ascii
from os import system
from astropy.time import Time
import argparse
from multiprocessing import Pool
import os
from datetime import datetime
import sys
import corner


# Initialize parameters

# In[2]:


def geninput(x):
    # x = {}
    # for n in range(len(varxnames)):
    #    x[varxnames[n]] = varx[n]
    # for n in range(len(constxnames)):
    #    x[constxnames[n]] = constx[n]
    tstar = x['tstar']
    td = x['td']
    tau = x['tau']
    thick = x['thick']
    ebv = x['ebv']
    idtype = x['idtype']
    fileuse = x['fileuse']
    output = open('foo1.inp', 'w')
    if (tstar < tstarmin or tstar > tstarmax) or (blackbody):
        output.write('Spectrum = 1\n')
        output.write('Number of BB = 1\n')
        output.write('Temperature = %0.2f\n' % (tstar))
    else:
        output.write('Spectrum = 5   \n')
        output.write('%s\n' % (fileuse))
    output.write('   optical properties index = 1 \n')
    output.write('   #   Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg \n')
    if idtype == 0:
        output.write('    x = 0.00    0.00   0.00    1.00    0.00    0.00 \n')
    else:
        output.write('    x = 0.00    0.00   1.00    0.00    0.00    0.00 \n')
    if custom_grain_distribution:
        output.write('- size distribution = 2  % custom       \n')
        output.write('  q = 3.5, a(min) = %s micron, a(max) = %s micron\n' % (amin, amax))
    else:
        output.write('- size distribution = 1  % standard MRN    \n')
    output.write('- temperature = %s K \n' % (td))
    output.write('- density type = 1                   \n')
    output.write('- number of powers = 1              \n')
    output.write('- shells relative thickness = %s\n' % (thick))
    output.write('- power = 2 \n')
    output.write('- grid type = 1                  % linear grid \n')
    output.write('- lambda0 = 0.55 micron          % optical depth specified  \n')
    output.write('- tau(min) = ' + str(tau) + ' ; tau(max) = 1000.0   % for the visual wavelength \n')
    output.write('- number of models = 1           \n')
    output.write('- accuracy for flux conservation = 0.05             \n')
    output.write('- verbosity flag;                              verbose = 1  \n')
    output.write('- properties of emerging spectra;            fname.spp = 1  \n')
    output.write('- detailed spectra for each model;          fname.s### = 1  \n')
    output.write('- images at specified wavelengths;          fname.i### = 1  \n')
    output.write('     number of wavelengths = 5  \n')
    output.write('     wavelengths = 3.5, 4.5, 6.0, 8.0, 24.0 micron  \n')
    output.write('- radial profiles for each model;           fname.r### = 1  \n')
    output.write('- detailed run-time messages;               fname.m### = 1  \n')
    output.write('- visibility function at spec. wavelengths; fname.v### = 0  \n')
    output.close()


# In[3]:


def log_gaussian_prior(x, pars):
    mu, sigma = pars
    return -np.log(sigma) - ((x - mu) ** 2) / (2 * sigma ** 2)


def log_linear_prior(x, pars):
    xmin, xmax = pars
    if xmin < x < xmax:
        return np.log(1 / (xmax - xmin))

    else:
        return -np.inf


def get_log_prior(pvals, pnames, pri_types, pri_params):
    '''
    pvals : array of vlues at which priors have to be evaluated
    pnames : names of parameters, not necessary?
    pri_types : array, with values gaussian or linear for each parameter
    pri_params : array of tuples, relevant parameters for each type of prior
    '''
    lp = 0
    for i in range(len(pnames)):
        if pri_types[i] == 'gaussian':
            lp = lp + log_gaussian_prior(pvals[i], pri_params[i])

        elif pri_types[i] == 'linear':
            lp = lp + log_linear_prior(pvals[i], pri_params[i])
    return lp


def run_dusty(varx):
    if verbose:
        pd = {}
        for n in range(len(varxnames)):
            pd[varxnames[n]] = varx[n]
        for n in range(len(constxnames)):
            pd[constxnames[n]] = constx[n]
        print('calling dusty with tstar = %0.1f; tau = %0.2f; td = %0.1f; thick = %0.2f; E(B-V) = %0.4f' % (
        pd['tstar'], pd['tau'], pd['td'], pd['thick'], pd['ebv']))
    system('./dusty')

    ierror = 0
    # try:
    data = ascii.read('foo1.stb')
    lam = data['col1']
    flx = data['col2']
    npt = len(lam)

    i = 0
    with open('foo1.out', 'r') as f:
        for line in f:
            if i == 42:
                # print 'line read in from foo1.out:', line
                line_s = line.split()
                id = int(line_s[0])
                tau0 = float(line_s[1])
                f1 = float(line_s[2])
                r1 = float(line_s[3])
                r1torstar = float(line_s[4])
                theta1 = float(line_s[5])
                tdout = float(line_s[6])
                break
            i += 1
            # except:
    #  ierror = 1
    #  lam = np.nan
    #  flx = np.nan
    #  npt = np.nan
    #  r1 = np.nan
    return lam, flx, npt, r1, ierror


def rl(rv, x):
    "rv is typically 3.1, x is 1/lambda [microns^-1]"
    if x < 1.1:
        a = 0.574 * x ** 1.61
        b = -0.527 * x ** 1.61
    else:
        if x < 3.3:
            y = x - 1.82
            a = (1.0 + 0.17699 * y - 0.50447 * y * y - 0.02427 * y ** 3 + 0.72085 * y ** 4 +
                 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7)
            b = (1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 -
                 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7)
        else:
            if x < 5.9:
                fa = 0.0
                fb = 0.0
            else:
                fa = -0.04473 * (x - 5.9) ** 2 - 0.009779 * (x - 5.9) ** 3
                fb = 0.21300 * (x - 5.9) ** 2 + 0.120700 * (x - 5.9) ** 3
            a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) + fa
            b = -3.090 + 1.825 * x + 1.206 / ((x - 4.67) ** 2 + 0.263) + fb
    rlval = rv * (a + b / rv)
    return rlval


def ext_corr(obsdat, ebv, verbose=True, veryverbose=False, debug=False):
    mlam = obsdat['mlam']
    mlumobs = obsdat['mlumobs']
    merrobs = obsdat['merrobs']
    filters = obsdat['filters']
    mlum = np.zeros(len(mlumobs))
    merr = np.zeros(len(mlumobs))
    if veryverbose: print('Extinction-corrected luminosities:')
    for i in range(len(mlam)):
        rlval = rl(rv, 1.0 / mlam[i])
        ecor = 10.0 ** (-0.4 * rlval * ebv)
        mlum[i] = mlumobs[i] / ecor
        merr[i] = merrobs[i] / ecor
    if veryverbose: print('    %s %s %s %s (%s)' % (mlam[i], mlum[i], merr[i], filters[i], ecor))
    gtzero = (mlum > 0)
    ltzero = np.invert(gtzero)
    mluml = np.zeros(len(mlumobs))
    merrl = np.zeros(len(mlumobs))
    mluml[gtzero] = np.log10(mlum[gtzero])
    merrl[gtzero] = merr[gtzero] / mlum[gtzero] / np.log(10)
    merrl[ltzero] = np.log10(merr[ltzero])
    return {'mlam': mlam, 'mluml': mluml, 'merrl': merrl, 'mlum': mlum, 'merr': merr, 'filters': filters}


def mol_abs_fractions(T, filters, abs_path='/scr2/viraj/Gattini/dusty_modeling/mol_abs'):
    print('Interpolating to calculate molecular absorption fraction')
    abs_fracs = []
    for f in filters:
        if os.path.exists('%s/%s_molecular_absorption_fractions.dat' % (abs_path, f)):
            fracs = ascii.read('%s/%s_molecular_absorption_fractions.dat' % (abs_path, f))
            abs_fracs.append(np.interp(T, fracs['col1'], fracs['col2']))
        elif os.path.exists('%s/%s_molecular_absorption_fractions.dat' % (abs_path, f.lower())):
            fracs = ascii.read('%s/%s_molecular_absorption_fractions.dat' % (abs_path, f.lower()))
            abs_fracs.append(np.interp(T, fracs['col1'], fracs['col2']))
        else:
            print('Could not find file for filter %s. Assuming no molecular absorption.' % (f))
            abs_fracs.append(1)
    return np.array(abs_fracs)


def calc_chi2(lam, flx, npt, r1, par_dict, data):
    '''
    data: extinction corrected data
    '''
    if 'ebv' in varxnames:
        corr_data = ext_corr(obsdata, par_dict['ebv'])
    else:
        corr_data = obsdata

    mlam = corr_data['mlam']
    mlum = corr_data['mlum']
    merr = corr_data['merr']
    mluml = corr_data['mluml']
    merrl = corr_data['merrl']
    filters_names = corr_data['filters']
    nm = len(mlum)
    thick = par_dict['thick']

    if corr_mol_absorption:
        abs_fracs = mol_abs_fractions(par_dict['tstar'], filters_names, path_to_abs_luts)

    '''
    if corr_mol_absorption:
        gtzero = (mlum > 0)
        ltzero = np.invert(gtzero)
        abs_fracs = mol_abs_fractions(par_dict['tstar'],filters_data[gtzero])
        mlum[gtzero] = mlum/abs_fracs
        merr[gtzero] = merr[gtzero]/abs_fracs
        mluml[gtzero] = mluml[gtzero] - np.log10(abs_fracs)
        merrl[gtzero] = merrl[gtzero] - np.log10(abs_fracs)
    '''

    aa = 0.0
    bb = 0.0
    val = np.zeros(100)
    vall = np.zeros(100)
    for i in range(nm):
        j = np.argmin(np.abs(lam - mlam[i]))
        if lam[j] > mlam[i]:
            j = j - 1
        # j = locate(lam,npt,mlam[i]) - 1    # -1 because indexed to 0
        val[i] = flx[j] + (flx[j + 1] - flx[j]) * (mlam[i] - lam[j]) / (lam[j + 1] - lam[j])
        if corr_mol_absorption and not lim_only:
            val[i] = val[i] * abs_fracs[i]
        if debug: print(flx[j], flx[j + 1], lam[j], lam[j + 1], mlam[i], val[i], mlum[i])
        if lim_only:
            aa = aa + (val[i] / merr[i]) ** 2
            if mlum[i] > 0:
                print('ERROR: I am expecting only upper limits')
                exit(1)
        else:
            vall[i] = np.log10(val[i] + 1.e-32)
            if debug: print('vall[%s] = %s' % (i, vall[i]))
            if mlum[i] > 0:
                aa = aa + (mluml[i] - vall[i]) / merrl[i] ** 2
                bb = bb + 1.0 / merrl[i] ** 2
    if lim_only:
        slum = np.sqrt(chilim / aa)
        sluml = np.log10(slum)
    else:
        sluml = aa / bb
    #if extrapolation:
    #    sluml = llum  # if extrapolating tau from a best-fit model then use the same luminosity
    if fixLstar:
        if verbose: print("forcing Lstar to 10^%0.2f" % (fixLstar))
        sluml = fixLstar
    # r1 output from DUSTY is the distance at which a point source w/luminosity 10^4 L_sun produces the bolometric flux F_e1
    r1 = np.log10(r1) + 0.5 * (sluml - 4.0)  # scale radius for luminosity
    r2 = r1 + np.log10(thick)
    # dusty reports a radius scaled to a luminosity of 10^4
    # so this scales it to the luminosity you just worked out
    chi = 0.0
    chis = {}
    if not extrapolation and not lim_only:
        for i in range(nm):
            if mlum[i] > 0:
                chis[filters[i]] = ((mluml[i] - sluml - vall[i]) / merrl[i]) ** 2
                chi = chi + chis[filters[i]]
        # chi = chi + ((mluml[i]-sluml-vall[i])/merrl[i])**2
        if debug: print('%s = chi + ((%s-%s-%s)/%s)**2' % (chi, mluml[i], sluml, vall[i], merrl[i]))
        if debug: print('sluml = %s' % (sluml))
        if mlum[i] < 0:
            rat = 10.0 ** (sluml + vall[i] - merrl[i])
            if debug: print('rat = %s = 10**(%s+%s-%s)' % (rat, sluml, vall[i], merrl[i]))
            chis[filters[i]] = rat * rat
            chi = chi + chis[filters[i]]
        # chi = chi + rat*rat
    chi_Lobs = chi
    if verbose: print('chi^2 from L_obs = %0.1f' % (chi_Lobs))
    eweight = 1.0
    chi = chi / eweight ** 2

    # chi from velocity
    if velocity:
        if dust_inout == 1:
            dustdist = r1
        else:
            dustdist = r2
        vlog = dustdist - dlog
        # chi = chi + ((vlog-vlog0)/evlog)**2
        if shell:
            chis['vlog'] = ((vlog - vlog0) / evlog) ** 2
            if verbose: print('chi^2 from velocity = %0.5f for Tdust = %0.3f K' % (chis['vlog'], par_dict['Td']))
            r_modeled = 10 ** dustdist
            r_actual = 10 ** vlog0 * 8.640e9 * telapsed
            if verbose: print('modeled r%d = %0.3g, actual r%d = %0.3g, diff = %0.3g for Tdust = %0.3f K' % (
            dust_inout, r_modeled, dust_inout, r_actual, r_actual - r_modeled, par_dict['Td']))
            if verbose:
                if r_actual - r_modeled > 0:
                    print('  try decreasing t_dust')
                else:
                    print('  try increasing t_dust')
        else:
            chis['vlog'] = 0
        chi = chi + chis['vlog']
        if verbose: print('sluml = %0.1f; newchi = %0.1f' % (sluml, chi))

    return -0.5 * chi, sluml


def generate_spectrum():
    return ''


def log_posterior(varx):
    par_dict = {}
    for i in range(len(varxnames)):
        par_dict[varxnames[i]] = varx[i]
    for i in range(len(constxnames)):
        par_dict[constxnames[i]] = constx[i]
    fileuse = generate_spectrum()
    par_dict['fileuse'] = fileuse

    lgpri = get_log_prior(varx, varxnames, varpritypes, varpripars)
    if lgpri == -np.inf:
        return lgpri, lgpri, lgpri

    geninput(par_dict)

    lam, flx, npt, r1, ierror = run_dusty(varx)

    chi2, sluml = calc_chi2(lam, flx, npt, r1, par_dict, obsdata)
    log_post = chi2 + lgpri
    return log_post, sluml, r1


def plot(sampler):
    samples = sampler.get_chain()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["Tstar", 'Td', 'tau']
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    plt.savefig('%s_chains.pdf' % (object_photometry_file))

    flat_samples = sampler.get_chain(flat=True)

    fig = corner.corner(
        flat_samples);
    plt.savefig('%s_corner.pdf' % (object_photometry_file))


def generate_data():
    # Generate simulated data
    wavs = [0.3, 0.5, 0.8, 3.6, 4.5, 5.8, 8.0]
    filts = ['g', 'r', 'i', 'F36', 'F45', 'F58', 'F80']
    lstar = 1e5
    pdict = {'tstar': 6000, 'td': 500, 'idtype': idtype, 'tau': 3.0, 'thick': thick, 'ebv': 0.2, 'fileuse': ''}
    geninput(pdict)
    lam, flx, npt, r1, ierror = run_dusty(pdict)
    simdatafile = 'sim_data_7pts_3pars_scat.dat'
    with open(simdatafile, 'w') as f:
        f.write('#Generated data with %s\n' % (pdict))
        for i in range(len(wavs)):
            j = np.argmin(np.abs(lam - wavs[i]))
            if lam[j] > wavs[i]:
                j = j - 1
            # j = locate(lam,npt,mlam[i]) - 1    # -1 because indexed to 0
            val = flx[j] + (flx[j + 1] - flx[j]) * (wavs[i] - lam[j]) / (lam[j + 1] - lam[j])
            val = flx[j] + (flx[j + 1] - flx[j]) * (wavs[i] - lam[j]) / (lam[j + 1] - lam[j])
            print('%s %s %.2f %s\n' % (
            wavs[i], np.random.normal(val * lstar, np.sqrt(val * lstar)), np.sqrt(val * lstar), filts[i]))
            f.write('%s %s %.2f %s\n' % (
            wavs[i], np.random.normal(val * lstar, np.sqrt(val * lstar)), np.sqrt(val * lstar), filts[i]))

    print('Generated data with', pdict)


def run_mcmc(dusty_inpdir):
    print('Changing directory to %s' % (dusty_inpdir))
    os.chdir(dusty_inpdir)

    seedind = int(dusty_inpdir.split('_')[-1]) ** 2
    np.random.seed(seedind)
    print('Seed for np is', seedind)
    state = np.random.RandomState(seedind)
    # print(state.get_state())
    outfilename = '%s_results.dat' % (object_photometry_file)
    if not continue_from_file:
        f = open(outfilename, "w")
        f.write('#Initial points %s\n' % (pos))
        f.write('#Priors on %s %s %s\n' % (varxnames, varpritypes, varpripars))
        f.close()

    print('Running emcee by varyng parameters %s with initial values %s and priors %s' % (varxnames, pos, varpripars))
    dtype = [("sluml", float), ("r1", float)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, blobs_dtype=dtype)
    print(sampler.random_state)
    # sampler.random_state.setter(state.get_state)
    ntrial = ntrials / 4
    for result in sampler.sample(pos, iterations=ntrial, progress=True, store=True):
        print(os.getcwd())
        position = result.coords
        lp = result.log_prob
        blobs = result.blobs
        sluml = blobs['sluml']
        r1s = blobs['r1']
        f = open(outfilename, "a")
        for k in range(len(position)):
            for j in position[k]:
                f.write('%.4f\t' % (j))
            f.write('%.4f\t' % (sluml[k]))
            f.write('%.4f\t' % (r1s[k]))
            f.write('%.4f\n' % (lp[k]))
        f.close()
    system('cd ..')


############################# Initialise parameters #########################################
#### Edit below ########
object_photometry_file = sys.argv[1]  # name of file giving photometry of the object
corr_mol_absorption = True  # Correct for molecular absorption? You need to provide lookup tables (absorption fraction vs temperature) in format xx_molecular_absorption_fractions.dat with xx being the filter
path_to_abs_luts = '/scr2/viraj/Gattini/dusty_modeling/mol_abs'  # Path to directory where luts are located

tstar = float(sys.argv[2])  # stellar temperature
ivarytstar = 1  # vary stellar temperature?
tprimin, tprimax = float(sys.argv[3]), float(sys.argv[4])
tstarrange = (tprimin, tprimax)
tstarmin, tstarmax = tstarrange
tstarprior = str(sys.argv[5])  # options : linear, gaussian

tau = 1  # tau
ivarytau = 1  # vary tau?
tauprior = 'linear'
taurange = (0.001, 7)

td = 600  # dust temperature
ivarytd = 1  # vary dust temperature?
tdprior = 'linear'
tdrange = (100, 1500)

thick = 10  # thickness
ivarythick = 0  # vary thickness?
thickprior = 'linear'
thickrange = (2, 100)

ebv = 0.0  # galaxy extinction E(B-V)
ivaryext = 0  # vary extinction
ebvprior = 'linear'
ebvrange = (0, 2.17)

idtype = 0  # dust type: 0=graphite 1=silicate

shell = 1  # shell/wind: 0=wind 1=shell

custom_grain_distribution = 0  # use a custom grain distribution? If 0 then use the standard MRN distribution
amin = 0.005  # micron, if providing custom grain distribution
amax = 0.25  # micron, if providing custom grain distribution

continue_from_file = 0  # continue interupted MCMC chain from the last line of outputfile
outfilename = '' # output file name to continue from
verbose = True
veryverbose = False
debug = False # Set to True if debugging
path_to_filter_files = '../dusty/MARCS/'
fixLstar = None
lim_only = False
chilim = 4.0 #chi2 for lim_only case
blackbody = True

#emcee parameters
nwalkers = 10
ntrials = 1000
nprocesses = 4

#####################################################################################################################
#####################################################################################################################
#### WARNING!!! The parameters below have not been implemented as of 2022-11-06, do not edit below #######
#####################################################################################################################
#####################################################################################################################
progenitor = 0  # if progenitor=1 use progenitor magnitude and write separate output
progenitor_photometry_file = 'premags.dat'  # name of file giving progenitor photometry

v0 = 765.0  # velocity [km/s]   # SN 1997bs according to Smith et al. 2011
vlog0 = np.log10(v0)
v0prior = ''
v0range = None

evlog = 0.3  # log error
evlogprior = ''
evlogrange = None
rv = 3.1
dust_inout = 0
velocity = 0

extrapolation = 0  # extrapolate model from a previously calculated date?
extrapolation_date = '2009-11-29'
extrapolation_name = '2009'

dLdt_limit = 0  # enforce a limit on dL/dt?
dLdt_filter = ['WFC3uvF555W']  # filter to apply dL/dt limit to
dLdt_obs = [320]  # [L_sun/yr]
dLdt_obs_err = [394]  # [L_sun/yr] (1-sigma uncertainty)

tstart_date = '1997-4-15'  # SN 1997bs
tnow_date = '2014-01-01'  # SN 1997bs observation
tnow = Time(tnow_date).jd
tstart = Time(tstart_date).jd
telapsed = tnow - tstart
telapsed_years = telapsed / 365.25
dlog = np.log10(8.640e9 * telapsed)  # coefficient 8.64e09 is (1 km/s)(day) so dlog = log(distance moved given time baseline at 1km/s)

##### Blank default Prior ranges and types for parameters where priors don't make sense, no need to edit ########
ivarydtype = 0
idtypeprior = ''
idtyperange = None
shellprior = ''
shellrange = None
customprior = ''
customrange = None
aminprior = ''
aminrange = None
amaxprior = ''
amaxrange = None
############################# End initialise parameters ##########################################

############################# Pre-process parameters #############################################
parnames = np.array(['tstar', 'td', 'tau', 'thick', 'v0', 'evlog', 'ebv', 'idtype', 'shell', 'custom_grain_distribution', 'amin',
     'amax'])
par_init = np.array([tstar, td, tau, thick, v0, evlog, ebv, idtype, shell, custom_grain_distribution, amin, amax])
vary = np.array([ivarytstar, ivarytd, ivarytau, ivarythick, 0, 0, ivaryext, ivarydtype, 0, 0, 0, 0])
pritypes = np.array(
    [tstarprior, tdprior, tauprior, thickprior, v0prior, evlogprior, ebvprior, idtypeprior, shellprior, customprior,
     aminprior, amaxprior])
pripars = np.array(
    [tstarrange, tdrange, taurange, thickrange, v0range, evlogrange, ebvrange, idtyperange, shellrange, customrange,
     aminrange, amaxrange])

varx_init = []
varxnames = []
varpripars = []
varpritypes = []

constx_init = []
constxnames = []
for i in range(len(vary)):
    if vary[i]:
        varx_init.append(par_init[i])
        varxnames.append(parnames[i])
        varpripars.append(pripars[i])
        varpritypes.append(pritypes[i])
    else:
        constx_init.append(par_init[i])
        constxnames.append(parnames[i])

varx_init = np.array(varx_init)
varxnames = np.array(varxnames)
varpripars = np.array(varpripars)
varpritypes = np.array(varpritypes)

constx = np.array(constx_init)
constxnames = np.array(constxnames)
print("Varying %s with initial value %s and prior range %s" % (varxnames, varx_init, varpripars))

data = ascii.read(object_photometry_file)  # ,delimiter=' ')
mlam = data['col1']
mlumobs = data['col2']
merrobs = data['col3']
filters = data['col4']
obsdat = {'mlam': mlam, 'mlumobs': mlumobs, 'merrobs': merrobs, 'filters': filters}
nm = np.shape(mlumobs)[0]
if verbose: print('read %d data points' % (nm))
if ivaryext == 0:
    obsdata = ext_corr(obsdat, ebv)
    ebv = 0.0
else:
    print('Will vary extinction\n')
    obsdata = obsdat
print(obsdata)
############################# End pre-process parameters #############################################


if __name__ == '__main__':
    ############################## Run emcee #############################################################
    ndim = len(varxnames)

    pos_c = []
    pos_size = []
    for i in range(len(varxnames)):
        pos_c.append(varx_init[i])
        pos_size.append(min(varx_init[i] - varpripars[i][0], varpripars[i][1] - varx_init[i]))
    pos_c = np.array(pos_c)
    pos_size = np.array(pos_size)
    psize = pos_size / 2
    pos = [pos_c + psize * (np.random.rand(ndim) - 0.5) for i in range(nwalkers)]
    print('pos_size is ', pos_size)

    if continue_from_file:
        pos = []
        rfile = ascii.read(outfilename)
        rfile = rfile[-nwalkers:]
        rcols = rfile.colnames[:ndim]
        for r in rfile:
            arr = []
            for colnm in rcols:
                arr.append(r[colnm])
            arr = np.array(arr)
            pos.append(arr)

    pool = Pool(nprocesses)
    print('Generating %s processes' % (nprocesses))
    tstart = datetime.utcnow()
    dusty_inpdir_list = ['proc_%s' % (i) for i in range(nprocesses)]
    for i in dusty_inpdir_list:
        if not os.path.exists(i):
            print('Making directory %s' % (i))
            os.mkdir(i)
            os.system('cp dusty %s/' % (i))
            os.system('cp dusty.inp %s/' % (i))
            os.system('cp dusty.f %s/' % (i))
            os.system('cp lambda_grid.dat %s/' % (i))

    pool.map(run_mcmc, dusty_inpdir_list)
    tend = datetime.utcnow()
    print(tend)
    print('Start:%s, End%s' % (tstart, tend))
