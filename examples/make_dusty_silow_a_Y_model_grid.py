import os
from pydusty.dusty import DustyParameters, Dusty
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger
from pathlib import Path
import numpy as np
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tau_wav_micron", type=float, default=0.55,
                        help="wavelength in um at which tau is specified")
    parser.add_argument('workdir', type=str,
                        default=None, help='dusty workdir name')
    parser.add_argument('--dusty_file_dir',
                        type=str, default='data/dusty_files',
                        help='Directory with dusty code files')
    parser.add_argument('--loglevel', type=str, default='DEBUG', help='logging level')
    parser.add_argument('--logfile', type=str, default=None, help='log file')

    args = parser.parse_args()

    logger = getLogger(args.loglevel, args.logfile)

    tstar_values = [2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000,
                    6500, 7000, 7500, 8000]

    tdust_values = [300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
                    800, 850, 900, 950, 1000, 1050, 1100]

    tau_values = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24,
                  0.26, 0.28, 0.3]

    shell_thickness_values = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

    a_values = 10**(np.linspace(-2, 0, 10))

    blackbody = Parameter(name='blackbody',
                          value=True)

    dust_type = Parameter(name='dust_type',
                          value='silow')
    tstarmin = Parameter(name='tstarmin',
                         value=3500)
    tstarmax = Parameter(name='tstarmin',
                         value=48999)
    custom_grain_distribution = Parameter(name='custom_grain_distribution',
                                          value=True)
    tau_wav_micron = Parameter(name='tau_wav', value=args.tau_wav_micron,
                               is_variable=False)

    workdir = args.workdir + f'/{dust_type.value}_grain_size_Y_grid'
    Path(workdir).mkdir(parents=True, exist_ok=True)
    params_list = []
    for tstarval in tstar_values:
        for tdustval in tdust_values:
            for tauval in tau_values:
                for a_value in a_values:
                    for shell_thickness_value in shell_thickness_values:
                        tstar = Parameter(name='tstar',
                                          value=tstarval,
                                          is_variable=False)

                        tdust = Parameter(name='tdust',
                                          value=tdustval,
                                          is_variable=True)

                        tau = Parameter(name=f'tau',
                                        value=tauval,
                                        is_variable=False)

                        grain_size = Parameter(name='min_grain_size',
                                                   value=a_value, is_variable=False)

                        shell_thickness = Parameter(name='shell_thickness',
                                                    value=shell_thickness_value)
                        params_list.append([tstar, tdust, tau, grain_size, shell_thickness])

    for params in tqdm(params_list, total=len(params_list)):
        tstar, tdust, tau, grain_size, shell_thickness = params
        dusty_parameters = DustyParameters(
            tstar=tstar,
            tdust=tdust,
            tau=tau,
            blackbody=blackbody,
            shell_thickness=shell_thickness,
            dust_type=dust_type,
            tstarmin=tstarmin,
            tstarmax=tstarmax,
            custom_grain_distribution=custom_grain_distribution,
            tau_wavelength_microns=tau_wav_micron,
            min_grain_size=grain_size,
            max_grain_size=grain_size,
        )

        dusty_runner = Dusty(parameters=dusty_parameters,
                             dusty_working_directory=workdir,
                             dusty_file_directory=args.dusty_file_dir
                             )

        os.chdir(workdir)
        dusty_runner.generate_input()
        dusty_runner.run()

        lam, flx, npt, r1, ierror = dusty_runner.get_results()
        with open(
                f'{workdir}/sed_{tstar.value}_{tdust.value}_{tau.value}_'
                f'{grain_size.value}_{shell_thickness.value}_{dust_type.value}.dat',
                'w') as f:
            f.write(f"# {r1}\n")
            f.write("lam, flux\n")
            for ind in range(len(lam)):
                f.write(f"{lam[ind]}, {flx[ind]}\n")
