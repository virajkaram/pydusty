import numpy as np
import logging
import os
from astropy.io import ascii
import argparse
import sys


logger = logging.getLogger(__name__)

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


def get_extinction_corrected_fluxes(mlam, mlumobs, merrobs, ebv, rv=3.1):
    mlum = np.zeros(len(mlumobs))
    merr = np.zeros(len(mlumobs))
    logger.debug('Extinction-corrected luminosities:')
    for i in range(len(mlam)):
        rlval = rl(rv, 1.0 / mlam[i])
        ecor = 10.0 ** (-0.4 * rlval * ebv)
        mlum[i] = mlumobs[i] / ecor
        merr[i] = merrobs[i] / ecor
        logger.debug('    %s %s %s (%s)' % (mlam[i], mlum[i], merr[i], ecor))

    return mlum, merr


def extinction_correct_obsdata(obsdat, ebv, rv=3.1):
    mlam = obsdat['mlam']
    mlumobs = obsdat['mlumobs_uncor']
    merrobs = obsdat['merrobs_uncor']

    mlum, merr = get_extinction_corrected_fluxes(mlam, mlumobs, merrobs, ebv, rv)
    gtzero = (mlum > 0)
    ltzero = np.invert(gtzero)
    mluml = np.zeros(len(mlumobs))
    merrl = np.zeros(len(mlumobs))
    mluml[gtzero] = np.log10(mlum[gtzero])
    merrl[gtzero] = merr[gtzero] / mlum[gtzero] / np.log(10)
    merrl[ltzero] = np.log10(merr[ltzero])

    obsdat['mluml'] = mluml
    obsdat['merrl'] = merrl
    obsdat['mlum'] = mlum
    obsdat['merr'] = merr
    return obsdat


def calculate_molecular_absorption_fractions(T, filters, abs_path='/scr2/viraj/Gattini/dusty_modeling/mol_abs'):
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


def load_and_extcor_data(object_photometry_file: str, ebv=0):
    data = ascii.read(object_photometry_file)
    mlam = data['col1']
    mlumobs = data['col2']
    merrobs = data['col3']
    filters = data['col4']
    obsdata = {'mlam': mlam, 'mlumobs_uncor': mlumobs, 'merrobs_uncor': merrobs, 'filters': filters}

    extinction_corrected_data = extinction_correct_obsdata(obsdata, ebv=ebv)
    return extinction_corrected_data


def get_default_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('object_photometry_file', type=str, help='photometry file')
    # parser.add_argument('--tstar', type=float, default=6000, help='Initial stellar temperature value')
    # parser.add_argument('--tstarmin', type=float, default=6000, help='Initial stellar temperature value')
    # parser.add_argument('--tstarmax', type=float, default=6000, help='Initial stellar temperature value')
    parser.add_argument('--workdir', type=str, default=None, help='dusty workdir name')
    parser.add_argument('--nwalkers', type=int, default=10, help='Number of emcee walkers')
    parser.add_argument('--ntrials', type=int, default=1000, help='Total number of emcee trials per walker')
    parser.add_argument('--nprocesses', type=int, default=1, help='Number of processors to use')
    parser.add_argument('--ebv', type=float, default=0, help='E_B-V to use')
    parser.add_argument('--molecular_absorption', action="store_true", help='Correct for molecular absorption?')
    parser.add_argument('--molecular_table_path', type=str, default=None, help='Path to molecular lookup tables')
    parser.add_argument('--limits_only',action="store_true", help='Only upper limits?')
    parser.add_argument('--fixLstar', action="store_true", help='fix luminosity?')
    parser.add_argument('--extrapolation', action="store_true", help='extrapolation?')
    parser.add_argument('--continue_from_file', action="store_true", help='Continue from an existing file?')
    parser.add_argument('--chi_square_limits_only', type=float, default=4.0, help='chisq. for limits only case')
    parser.add_argument('--loglevel', type=str, default='DEBUG', help='logging level')
    parser.add_argument('--logfile', type=str, default=None, help='log file')

    return parser


def getLogger(level, logfile=None):
    logger = logging.getLogger("pydusty")
    formatter = logging.Formatter('%(name)s [l %(lineno)d] - %(levelname)s - %(message)s')
    if logfile is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(logfile)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(level)
    return logger