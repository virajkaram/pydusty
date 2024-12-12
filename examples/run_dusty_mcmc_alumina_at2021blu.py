from pydusty.dusty import DustyParameters, Dusty_Alumina_SilOW
from pydusty.mcmc import Emcee
from pydusty.parameters import Parameter
from pydusty.priors import UniformPrior, GaussianPrior
from pydusty.utils import get_default_argparser, load_and_extcor_data, getLogger
from pydusty.parallel import ParallelEmceeRunner


def run_mcmc_process(emcee_runner: Emcee):
    emcee_runner.run_mcmc()


if __name__ == '__main__':

    argparser = get_default_argparser()
    args = argparser.parse_args()
    logger = getLogger(args.loglevel, args.logfile)

    # Set initial parameters
    tstar = Parameter(name='tstar',
                      value=3000,
                      is_variable=True,
                      prior=UniformPrior(prior_parameters=(1000, 6000)),
                      )

    tdust = Parameter(name='tdust',
                      value=1100,
                      is_variable=True,
                      prior=UniformPrior(prior_parameters=(700, 1500))
                      )

    tau = Parameter(name='tau',
                    value=0.05,
                    prior=UniformPrior(prior_parameters=(1e-3, 0.1)),
                    is_variable=True)

    tau_wav_micron = Parameter(name='tau_wav', value=100,
                               is_variable=False)

    blackbody = Parameter(name='blackbody',
                          value=True)

    shell_thickness = Parameter(name='shell_thickness',
                                value=2.0, is_variable=False)

    dust_type = Parameter(name='dust_type',
                          value=f'si_al_tau_{tau_wav_micron.value}um')

    tstarmin = Parameter(name='tstarmin',
                         value=3500)

    tstarmax = Parameter(name='tstarmin',
                         value=48999)

    custom_grain_distribution = Parameter(name='custom_grain_distribution',
                                          value=False)

    al_abundance = Parameter(name='al', value=0.5,
                             prior=UniformPrior(prior_parameters=(0, 1)),
                             is_variable=True)

    # error_underestimate_scale_factor = Parameter(name='error_underestimate_scale_factor',
    #                                              value=0.15,
    #                                              prior=UniformPrior(
    #                                                  prior_parameters=(0.05, 0.3)
    #                                              ),
    #                                              is_variable=True
    #                                              )

    error_underestimate_scale_factor = Parameter(
        name='error_underestimate_scale_factor',
        value=0,
        is_variable=False
        )

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
        error_underestimate_factor=error_underestimate_scale_factor,
    )

    working_dir = args.workdir

    ext_corrected_obsdata = load_and_extcor_data(args.object_photometry_file)

    parallel_emcee_runner = ParallelEmceeRunner(dusty=Dusty_Alumina_SilOW,
                                                dusty_parameters=dusty_parameters,
                                                nwalkers=args.nwalkers,
                                                working_dir=working_dir,
                                                obsdata=ext_corrected_obsdata,
                                                object_photometry_file=args.object_photometry_file,
                                                nprocesses=args.nprocesses,
                                                ntrials=args.ntrials,
                                                correct_for_molecular_absorption=args.molecular_absorption,
                                                molecular_absorption_lookup_table_path=args.molecular_table_path,
                                                limits_only=args.limits_only,
                                                fixLstar=args.fixLstar,
                                                chi_square_limits_only=args.chi_square_limits_only,
                                                extrapolation=args.extrapolation,
                                                continue_from_file=args.continue_from_file,
                                                dusty_file_dir=args.dusty_file_dir,
                                                )
    parallel_emcee_runner.run()
    parallel_emcee_runner.make_plots()
    emcee_runners = []
