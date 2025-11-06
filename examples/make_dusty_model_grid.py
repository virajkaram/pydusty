import os

from pydusty.dusty import DustyParameters, Dusty_Alumina_SilDL
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger
from pathlib import Path
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tau_wav_micron", type=float, default=100.0,
                        help="wavelength in um at which tau is specified")
    parser.add_argument("--thick", type=float, default=2.0)
    parser.add_argument("--al", type=float, default=0.5,
                        help="Aluminum abundance, Silicon abundance is 1-al")
    parser.add_argument("--al_type", type=str, default="compact",
                        choices=['compact', 'porous'])
    parser.add_argument('workdir', type=str, default=None, help='dusty workdir name')
    parser.add_argument('--dusty_file_dir', type=str, default='data/dusty_files',
                        help='Directory with dusty code files')
    parser.add_argument('--loglevel', type=str, default='DEBUG', help='logging level')
    parser.add_argument('--logfile', type=str, default=None, help='log file')

    args = parser.parse_args()

    logger = getLogger(args.loglevel, args.logfile)

    # tstar_values = [1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750,
    #  4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500,]
    tstar_values = [1500, 2000, 2500, 3000, 3500, 4000]
    # tdust_values = [500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100,
    #                 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500]
    tdust_values = [100, 200, 300, 400, 500, 600, 700]
    tau_values =10**np.linspace(1, np.log10(500), 10)
    blackbody = Parameter(name='blackbody',
                          value=True)
    shell_thickness = Parameter(name='shell_thickness',
                                value=args.thick)
    dust_type = Parameter(name='dust_type',
                          value=f'si_{(1 - args.al)}_al_{args.al}_'
                                f'{args.al_type}_tau_{args.tau_wav_micron}um')
    tstarmin = Parameter(name='tstarmin',
                         value=3500)
    tstarmax = Parameter(name='tstarmin',
                         value=48999)
    custom_grain_distribution = Parameter(name='custom_grain_distribution',
                                          value=False)
    tau_wav_micron = Parameter(name='tau_wav', value=args.tau_wav_micron,
                               is_variable=False)
    al_abundance = Parameter(name='al', value=args.al, is_variable=False)

    workdir = args.workdir + f'/{dust_type.value}_thick_{shell_thickness.value}'
    Path(workdir).mkdir(parents=True, exist_ok=True)
    for tstarval in tstar_values:
        for tdustval in tdust_values:
            for tauval in tau_values:

                tstar = Parameter(name='tstar',
                                  value=tstarval,
                                  is_variable=False)

                tdust = Parameter(name='tdust',
                                  value=tdustval,
                                  is_variable=True)

                tau = Parameter(name=f'tau',
                                value=tauval,
                                is_variable=False)

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
                    al_com_abundance=al_abundance,
                )

                dusty_runner = Dusty_Alumina_SilDL(parameters=dusty_parameters,
                                                   dusty_working_directory=workdir,
                                                   dusty_file_directory=args.dusty_file_dir
                                                   )

                os.chdir(workdir)
                dusty_runner.generate_input()
                dusty_runner.run()

                lam, flx, npt, r1, ierror = dusty_runner.get_results()
                with open(
                        f'{workdir}/sed_{tstar.value}_{tdust.value}_{tau.value}_{dust_type.value}_{shell_thickness.value}.dat',
                        'w') as f:
                    f.write(f"# {r1}\n")
                    f.write("lam, flux\n")
                    for ind in range(len(lam)):
                        f.write(f"{lam[ind]}, {flx[ind]}\n")
