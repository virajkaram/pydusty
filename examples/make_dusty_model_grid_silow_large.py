import os
from pydusty.dusty import DustyParameters, Dusty
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger
from pathlib import Path
from glob import glob
from multiprocessing import Pool
import numpy as np


def run_dusty_for_params(params):
    (tstar, tdust, tau,
    shell_thickness, dust_type,
    work_subdir,
    dusty_file_dir,
    working_dir) = params

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
    )

    dusty_runner = Dusty(parameters=dusty_parameters,
                         dusty_working_directory=work_subdir,
                         dusty_file_directory=dusty_file_dir
                        )

    base_filename = (f'sed_{tstar.value}_{tdust.value}_{tau.value}_'
                     f'{dust_type.value}_{shell_thickness.value}_'
                     f'{tau_wav_micron.value}um.dat')
    filename = (f'{work_subdir}/{base_filename}')

    os.chdir(work_subdir)
    dusty_runner.generate_input()
    dusty_runner.run()

    lam, flx, npt, r1, ierror = dusty_runner.get_results()
    with open(filename, 'w') as f:
        f.write(f"# {r1}\n")
        f.write("lam, flux\n")
        for ind in range(len(lam)):
            f.write(f"{lam[ind]}, {flx[ind]}\n")

    # Cleanup directory - remove all files except the sed file
    for file in glob(f'{work_subdir}/*'):
        if file != filename:
            os.remove(file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('workdir', type=str,
                        default=None, help='dusty workdir name')
    parser.add_argument('--dusty_file_dir',
                        type=str, default='data/dusty_files',
                        help='Directory with dusty code files')
    parser.add_argument('--loglevel', type=str, default='DEBUG',
                        help='logging level')
    parser.add_argument('--logfile', type=str, default=None,
                        help='log file')
    parser.add_argument('--ncpus', type=int, default=12,)

    args = parser.parse_args()

    logger = getLogger(args.loglevel, args.logfile)

    # tstar values 2000 to 7500 in steps of 100
    tstar_values = np.arange(2000, 7501, 100)

    # tdust values 100 to 1500 in steps of 100
    tdust_values = np.arange(100, 1501, 100)

    #tdust values 0.01 to 100 in log steps
    tau_values = np.logspace(-2, 2, 20)

    shell_thickness_values = np.array([2.0, 5.0, 10.0, 20.0, 40.0, 50.0])

    # Total parameters :
    n_cpus = args.ncpus

    blackbody = Parameter(name='blackbody',
                          value=True)

    tstarmin = Parameter(name='tstarmin',
                         value=2000)
    tstarmax = Parameter(name='tstarmin',
                         value=48999)
    custom_grain_distribution = Parameter(name='custom_grain_distribution',
                                          value=False)
    tau_wav_micron = Parameter(name='tau_wav', value=0.55,
                               is_variable=False)

    dust_type = Parameter(name='dust_type',
                          value='silow')

    working_dir = args.workdir + f'/silow_grid_large'
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    params_list = []
    i = 0
    for tstarval in tstar_values:
        for tdustval in tdust_values:
            for tauval in tau_values:
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

                        shell_thickness = Parameter(name='shell_thickness',
                                                    value=shell_thickness_value)
                        base_filename = (f'sed_{tstar.value}_{tdust.value}_{tau.value}_'
                                         f'{dust_type.value}_{shell_thickness.value}_'
                                         f'{tau_wav_micron.value}um.dat')
                        if len(glob(f'{working_dir}/*/{base_filename}')) > 0:
                            continue
                        i+=1
                        work_subdir = f'{working_dir}/run_{i}'

                        params_list.append([tstar, tdust, tau,
                                            shell_thickness, dust_type,
                                            work_subdir,
                                            args.dusty_file_dir,
                                            working_dir])

    print("Total number of models to run: ", len(params_list))
    print("OK? (y/n)")
    if input().lower() != 'y':
        print("Aborting.")
        exit(0)
    pool = Pool(processes=n_cpus)
    pool.map(run_dusty_for_params, params_list)
