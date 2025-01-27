import os

from pydusty.dusty import DustyCustomInputSpectrum, DustyParameters
from pydusty.parameters import Parameter
import argparse
from pydusty.utils import getLogger

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--tau", type=float, default=50)
    parser.add_argument("--tau_wav_micron", type=float, default=100.0,
                        help="wavelength in um at which tau is specified")
    parser.add_argument("--tdust", type=float, default=1000)
    parser.add_argument("--thick", type=float, default=2.0)
    parser.add_argument("--dtype", choices=['graphite', 'silicate',
                                            'amorphous_carbon', 'silicate_carbide'],
                        default='graphite')
    parser.add_argument('workdir', type=str, default=None, help='dusty workdir name')
    parser.add_argument('--dusty_file_dir', type=str, default='data/dusty_files',
                        help='Directory with dusty code files')
    parser.add_argument('--loglevel', type=str, default='DEBUG', help='logging level')
    parser.add_argument('--logfile', type=str, default=None, help='log file')

    args = parser.parse_args()

    logger = getLogger(args.loglevel, args.logfile)


    custom_input_spectrum_file = Parameter(name='custom_input_spectrum_file',
                                            value=args.input_file,
                                            is_variable=False)
    tdust = Parameter(name='tdust',
                      value=args.tdust,
                      is_variable=True)

    tau = Parameter(name='tau',
                    value=args.tau,
                    is_variable=False)

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
                               dusty_working_directory=args.workdir,
                               dusty_file_directory=args.dusty_file_dir
                               )

    os.chdir(args.workdir)
    dusty_runner.generate_input()
    dusty_runner.run()

    lam, flx, npt, r1, ierror = dusty_runner.get_results()
    with open(f"{args.workdir}/sed_{args.input_file.split('/')[-1].split('.dat')[0]}_"
              f"{tdust.value}_{tau.value}_{dust_type.value}_{shell_thickness.value}_"
              f"{tau_wav_micron.value}um.dat", 'w') as f:
        f.write(f"# {r1}\n")
        f.write("lam, flux\n")
        for ind in range(len(lam)):
            f.write(f"{lam[ind]}, {flx[ind]}\n")
