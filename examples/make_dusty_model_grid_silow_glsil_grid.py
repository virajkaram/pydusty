import os
from pydusty.dusty import DustyParameters, Dusty_Multi_Composition
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger
from pathlib import Path
from glob import glob
from tqdm import tqdm
from multiprocessing import Pool


def run_dusty_for_params(params):
    (tstar, tdust, tau, grain_size, shell_thickness, dust_abundances, dust_types, dust_type_param,
     work_subdir, dusty_file_dir, working_dir) = params

    dusty_parameters = DustyParameters(
        tstar=tstar,
        tdust=tdust,
        tau=tau,
        blackbody=blackbody,
        shell_thickness=shell_thickness,
        dust_type=dust_type_param,
        tstarmin=tstarmin,
        tstarmax=tstarmax,
        custom_grain_distribution=custom_grain_distribution,
        tau_wavelength_microns=tau_wav_micron,
        dust_composition_elements=dust_types,
        dust_composition_abundances=dust_abundances,
        min_grain_size=grain_size,
        max_grain_size=grain_size,
    )

    dusty_runner = Dusty_Multi_Composition(parameters=dusty_parameters,
                                           dusty_working_directory=work_subdir,
                                           dusty_file_directory=dusty_file_dir
                                           )

    base_filename = (f'sed_{tstar.value}_{tdust.value}_{tau.value}_'
                     f'{grain_size.value}_{shell_thickness.value}_'
                     f'{dust_type_param.value}.dat')
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

    tstar_values = [4000, 4500, 3500, 5000, 3000]

    tdust_values = [800, 900, 700, 1000, 600]

    tau_values = [0.1, 0.14, 0.06, 0.18, 0.02]

    shell_thickness_values = [2.0, 5.0, 10.0, 20.0, 40.0]

    a_values = [0.01, 0.05, 0.1, 0.15, 0.20]

    si_ratios = [0.6, 0.7, 0.8, 0.9, 1.0]

    n_cpus = args.ncpus

    blackbody = Parameter(name='blackbody',
                          value=True)

    tstarmin = Parameter(name='tstarmin',
                         value=3500)
    tstarmax = Parameter(name='tstarmin',
                         value=48999)
    custom_grain_distribution = Parameter(name='custom_grain_distribution',
                                          value=True)
    tau_wav_micron = Parameter(name='tau_wav', value=0.55,
                               is_variable=False)

    working_dir = args.workdir + f'/silow_glsil_grain_size_Y_grid'
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    params_list = []
    i = 0
    for tstarval in tstar_values:
        for tdustval in tdust_values:
            for tauval in tau_values:
                for a_value in a_values:
                    for shell_thickness_value in shell_thickness_values:
                        for si_ratio in si_ratios:
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

                            dust_type_param = Parameter(name='dust_type',
                                                  value=f'silow_{si_ratio}_glsil_tau_'
                                                        f'{1-si_ratio}_{tau_wav_micron.value}um')

                            dust_abundances = [si_ratio, 1 - si_ratio]
                            dust_types = ['silow', 'glassy_silicate']

                            i+=1
                            work_subdir = f'{working_dir}/{i % n_cpus}'

                            Path(work_subdir).mkdir(parents=True, exist_ok=True)
                            params_list.append([tstar, tdust, tau, grain_size,
                                                shell_thickness, dust_abundances,
                                                dust_types, dust_type_param,
                                                work_subdir, args.dusty_file_dir,
                                                working_dir])

    pool = Pool(processes=n_cpus)
    pool.map(run_dusty_for_params, params_list)
