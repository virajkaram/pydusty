import os

from pydusty.dusty import DustyParameters, DustyCustomInputSpectrum
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger
from pathlib import Path
import numpy as np
from glob import glob
from astropy.table import Table
from astropy.convolution import Box1DKernel, convolve
from multiprocessing import Pool


def convert_phoenix_file_to_dusty_format(phoenix_filename, outfilename, log_gkey='g00'):
    phnx = Table.read(phoenix_filename)
    # Restrict to wavelengths between 0.1 and 20 microns
    phnx['wav_um'] = phnx['WAVELENGTH'] / 1e4
    phnx = phnx[(phnx['wav_um']>0.1) & (phnx['wav_um']<50)]
    phnx_flx = np.array(phnx[log_gkey])
    phnx_wavs = phnx['wav_um']
    phnx_flx /= np.nanmax(phnx_flx)
    # convolve to fewer than 10000 points
    kernel = Box1DKernel(20)
    phnx_flx = convolve(phnx_flx, kernel)
    phnx_wavs_short = phnx_wavs[::10]
    phnx_flx_short = phnx_flx[::10]
    assert len(phnx_wavs_short) < 10000
    with open(outfilename, 'w') as f:
        f.write(f'# Custom input spectrum file for dusty from {phoenix_filename}\n')
        f.write('# lambda F_lambda\n')
        f.write('# (micron) (arbitrary)\n')
        for wav, flux in zip(phnx_wavs_short, phnx_flx_short):
            f.write(f'{wav} {flux}\n')


def run_dusty_for_params(params):
    (custom_input_spectrum_file,tstarval, tdust, tau, blackbody, shell_thickness, dust_type, tstarmin,
     tstarmax, custom_grain_distribution, tau_wav_micron, dusty_file_dir, work_subdir) = params

    dusty_parameters = DustyParameters(
        custom_input_spectrum_file=custom_input_spectrum_file,
        tdust=tdust,
        tau=tau,
        blackbody=blackbody,
        shell_thickness=shell_thickness,
        dust_type=dust_type,
        tstarmin=tstarmin,
        tstarmax=tstarmax,
        custom_grain_distribution=custom_grain_distribution,
        tau_wavelength_microns=tau_wav_micron,
    )

    dusty_runner = DustyCustomInputSpectrum(parameters=dusty_parameters,
                                            dusty_working_directory=work_subdir,
                                            dusty_file_directory=dusty_file_dir
                                            )

    base_filename = (f'sed_{tstarval}_{tdust.value}_{tau.value}_'
                     f'{dust_type.value}_{shell_thickness.value}_{tau_wav_micron.value}um.dat')
    filename = (f'{work_subdir}/{base_filename}')
    if len(glob(f'{working_dir}/*/{base_filename}')) > 0:
        return

    os.chdir(work_subdir)
    dusty_runner.generate_input()
    dusty_runner.run()

    lam, flx, npt, r1, ierror = dusty_runner.get_results()
    with open(filename, 'w') as f:
        f.write(f"# {r1}\n")
        f.write("lam, flux\n")
        for ind in range(len(lam)):
            f.write(f"{lam[ind]}, {flx[ind]}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--phoenix_directory', type=str, default='/scr2/viraj/phoenix_files',)
    parser.add_argument("--log_gkey", type=str,
                        default='g00')
    parser.add_argument('--skip_remaking_dusty_format', action='store_true')
    parser.add_argument("--tau_wav_micron", type=float, default=100.0,
                        help="wavelength in um at which tau is specified")
    parser.add_argument("--thick", type=float, default=2.0)
    parser.add_argument("--dtype", choices=['graphite', 'silicate',
                                            'amorphous_carbon', 'silicate_carbide',
                                            'silow'],
                        default='silow')
    parser.add_argument('workdir', type=str, default=None, help='dusty workdir name')
    parser.add_argument('--dusty_file_dir', type=str, default='data/dusty_files',
                        help='Directory with dusty code files')
    parser.add_argument('--loglevel', type=str, default='DEBUG', help='logging level')
    parser.add_argument('--logfile', type=str, default=None, help='log file')
    parser.add_argument('--ncpus', type=int, default=12, )

    args = parser.parse_args()

    logger = getLogger(args.loglevel, args.logfile)

    tstar_values = [2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
                    3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600,
                    4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700,
                    5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700,]
    tdust_values = [200, 300, 400, 500,  600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]

    tau_values = 10**np.linspace(-0.5, 2.0, 15)
    log_gkey = args.log_gkey

    phoenix_dusty_format_dirname = f"{args.phoenix_directory}/phoenix_{log_gkey}_dusty_format"
    Path(phoenix_dusty_format_dirname).mkdir(parents=True, exist_ok=True)
    if not args.skip_remaking_dusty_format:
        phoenix_filelist = glob(f"{args.phoenix_directory}/phoenix*.fits")
        for phoenix_filename in phoenix_filelist:
            convert_phoenix_file_to_dusty_format(phoenix_filename,
                                                 f"{phoenix_dusty_format_dirname}/dusty_{Path(phoenix_filename).name.replace('.fits', '.dat')}",
                                                 log_gkey=log_gkey,
                                                 )

    ncpus = args.ncpus

    blackbody = Parameter(name='blackbody',
                          value=False)
    shell_thickness = Parameter(name='shell_thickness',
                                value=args.thick)
    dust_type = Parameter(name='dust_type',
                          value=args.dtype)
    tstarmin = Parameter(name='tstarmin',
                         value=3500)
    tstarmax = Parameter(name='tstarmin',
                         value=48999)
    custom_grain_distribution = Parameter(name='custom_grain_distribution',
                                          value=False)
    tau_wav_micron = Parameter(name='tau_wav', value=args.tau_wav_micron,
                               is_variable=False)

    workdir = args.workdir + f'/{dust_type.value}_thick_{shell_thickness.value}'
    Path(workdir).mkdir(parents=True, exist_ok=True)

    # Check that we have phoenix files for each Tstar, otherwise raise an error
    for tstarval in tstar_values:
        phoenix_file = glob(f'{phoenix_dusty_format_dirname}/dusty_phoenix*_{tstarval}.dat')
        if len(phoenix_file) == 0:
            raise FileNotFoundError(f"Could not find phoenix file for Tstar = {tstarval}"
                                    f"in directory {phoenix_dusty_format_dirname}")

    working_dir = args.workdir + '/phoenix_input_spectra_grid_g00'
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    params_list = []
    i = 0
    for tstarval in tstar_values:
        for tdustval in tdust_values:
            for tauval in tau_values:
                phoenix_file = glob(
                    f'{phoenix_dusty_format_dirname}/dusty_phoenix*_{tstarval}.dat')[0]
                custom_input_spectrum_file = Parameter(
                    name='custom_input_spectrum_file',
                    value=phoenix_file,
                    is_variable=False)

                tdust = Parameter(name='tdust',
                                  value=tdustval,
                                  is_variable=True)

                tau = Parameter(name='tau',
                                value=tauval,
                                is_variable=False)

                i += 1
                work_subdir = f'{working_dir}/{i % ncpus}'
                params_list.append([custom_input_spectrum_file, tstarval, tdust, tau,
                                    blackbody, shell_thickness, dust_type, tstarmin,
                                    tstarmax, custom_grain_distribution, tau_wav_micron,
                                    args.dusty_file_dir,
                                    work_subdir,])
    pool =  Pool(processes=ncpus)
    pool.map(run_dusty_for_params, params_list)
