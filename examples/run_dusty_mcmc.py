from pydusty.dusty import Dusty, DustyParameters
from pydusty.mcmc import Emcee
from pydusty.parameters import Parameter
from pydusty.priors import GaussianPrior, UniformPrior
from pydusty.utils import get_default_argparser, load_and_extcor_data, getLogger
import logging
import sys


if __name__ == '__main__':

    argparser = get_default_argparser()
    args = argparser.parse_args()
    logger = getLogger(args.loglevel, args.logfile)

    # Set initial parameters
    tstar = Parameter(name='tstar',
                      value=4000,
                      prior=UniformPrior(prior_parameters=(3000, 6000)),
                      is_variable=True)

    tdust = Parameter(name='tdust',
                      value=1000,
                      prior=UniformPrior(prior_parameters=(500, 2000)),
                      is_variable=True)

    tau = Parameter(name='tau',
                    value=1,
                    prior=UniformPrior(prior_parameters=(0, 5)),
                    is_variable=True)

    blackbody = Parameter(name='blackbody',
                          value=True)
    shell_thickness = Parameter(name='shell_thickness',
                                value=2.0)
    dust_type = Parameter(name='dust_type',
                          value=1)
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

    working_dir = args.workdir
    if working_dir is None:
        working_dir = './run'

    basic_dusty = Dusty(parameters=dusty_parameters,
                        dusty_working_directory=working_dir)

    ext_corrected_obsdata = load_and_extcor_data(args.object_photometry_file)

    emcee_runner = Emcee(obsdata=ext_corrected_obsdata,
                         object_photometry_file=args.object_photometry_file,
                         nwalkers=args.nwalkers,
                         ntrials=args.ntrials,
                         dusty=basic_dusty,
                         correct_for_molecular_absorption=args.molecular_absorption,
                         molecular_absorption_lookup_table_path=args.molecular_table_path,
                         limits_only=args.limits_only,
                         fixLstar=args.fixLstar,
                         chi_square_limits_only=args.chi_square_limits_only,
                         extrapolation=args.extrapolation,
                         continue_from_file=args.continue_from_file,
                         )

    emcee_runner.run_mcmc()
