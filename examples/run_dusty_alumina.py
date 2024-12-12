import os

from pydusty.dusty import DustyParameters, Dusty_Alumina_SilDL
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tstar", type=float, default=5000)
    parser.add_argument("--tau", type=float, default=50)
    parser.add_argument("--tau_wav_micron", type=float, default=100.0,
                        help="wavelength in um at which tau is specified")
    parser.add_argument("--tdust", type=float, default=1000)
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

    tstar = Parameter(name='tstar',
                      value=args.tstar,
                      is_variable=False)

    tdust = Parameter(name='tdust',
                      value=args.tdust,
                      is_variable=True)

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

    tau = Parameter(name=f'tau',
                    value=args.tau,
                    is_variable=False)

    al_abundance = Parameter(name='al', value=args.al, is_variable=False)
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
                                       dusty_working_directory=args.workdir,
                                       dusty_file_directory=args.dusty_file_dir
                                       )

    os.chdir(args.workdir)
    dusty_runner.generate_input()
    dusty_runner.run()

    lam, flx, npt, r1, ierror = dusty_runner.get_results()
    with open(
            f'{args.workdir}/sed_{tstar.value}_{tdust.value}_{tau.value}_{dust_type.value}_{shell_thickness.value}.dat',
            'w') as f:
        f.write(f"# {r1}\n")
        f.write("lam, flux\n")
        for ind in range(len(lam)):
            f.write(f"{lam[ind]}, {flx[ind]}\n")
