import os

from pydusty.dusty import Dusty, DustyParameters
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tstar", type=float, default=5000)
    parser.add_argument("--tau", type=float, default=50)
    parser.add_argument("--tdust", type=float, default=1000)
    parser.add_argument("--thick", type=float, default=2.0)
    parser.add_argument("--dtype", choices=['graphite', 'silicate'], default='graphite')
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

    tau = Parameter(name='tau',
                    value=args.tau,
                    is_variable=False)

    blackbody = Parameter(name='blackbody',
                          value=True)
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

    dusty_parameters = DustyParameters(
        tstar=tstar,
        tdust=tdust,
        tau=tau,
        blackbody=blackbody,
        shell_thickness=shell_thickness,
        dust_type=dust_type,
        tstarmin=tstarmin,
        tstarmax=tstarmax,
        custom_grain_distribution=custom_grain_distribution
    )

    dusty_runner = Dusty(parameters=dusty_parameters,
                         dusty_working_directory=args.workdir,
                         dusty_file_directory=args.dusty_file_dir
                         )

    os.chdir(args.workdir)
    dusty_runner.generate_input()
    dusty_runner.run()

    lam, flx, npt, r1, ierror = dusty_runner.get_results()
    with open(f'{args.workdir}/sed_{tstar.value}_{tdust.value}_{tau.value}_{dust_type.value}.dat','w') as f:
        f.write("lam, flux\n")
        for ind in range(len(lam)):
            f.write(f"{lam[ind]}, {flx[ind]}\n")

