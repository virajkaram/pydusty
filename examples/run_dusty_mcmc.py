from pydusty.dusty import Dusty, DustyParameters
from pydusty.mcmc import Emcee
from pydusty.parameters import Parameter
from pydusty.priors import UniformPrior
from pydusty.utils import get_default_argparser, load_and_extcor_data, getLogger
from multiprocessing import Pool
from datetime import datetime


def run_mcmc_process(emcee_runner: Emcee):
    emcee_runner.run_mcmc()


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

    working_dir = args.workdir
    if working_dir is None:
        working_dir = './run'

    nprocesses = args.nprocesses
    logger.info(f"Creating {nprocesses} processes.")
    pool = Pool(nprocesses)
    tstart = datetime.utcnow()
    emcee_runners = []

    for process_num in range(nprocesses):
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

        dusty_process_working_dir = working_dir + f'_{process_num}'
        basic_dusty = Dusty(parameters=dusty_parameters,
                            dusty_working_directory=dusty_process_working_dir,
                            )

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
                             random_seed=process_num
                             )

        emcee_runners.append(emcee_runner)

    pool.map(run_mcmc_process, emcee_runners)
    tend = datetime.utcnow()
    logger.info(f'Start {tstart}, End {tend}')