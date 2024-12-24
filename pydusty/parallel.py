import os.path

from pydusty.mcmc import Emcee
from pydusty.dusty import Dusty, BaseDusty
from multiprocessing import Pool
from pydusty.parameters import DustyParameters
from datetime import datetime
import logging
from astropy.table import Table, vstack
from astropy.io import ascii
from pydusty.utils import make_corner_plot

logger = logging.getLogger(__name__)


class ParallelEmceeRunner:
    def __init__(self,
                 nwalkers: int,
                 dusty_parameters: DustyParameters,
                 working_dir: str,
                 obsdata: dict,
                 object_photometry_file: str,
                 nprocesses: int,
                 ntrials: int = 250,
                 correct_for_molecular_absorption: bool =False,
                 molecular_absorption_lookup_table_path: str = None,
                 limits_only: bool = False,
                 fixLstar: bool = False,
                 chi_square_limits_only: float = 4,
                 extrapolation: bool = False,
                 continue_from_file: bool = False,
                 dusty_file_dir: str = 'data/dusty_files',
                 dusty: BaseDusty = Dusty,
                 ):
        self.nwalkers = nwalkers
        self.dusty_parameters = dusty_parameters
        self.working_dir = working_dir
        self.obsdata = obsdata
        self.correct_for_molecular_absorption = correct_for_molecular_absorption
        self.limits_only = limits_only
        self.fixLstar = fixLstar
        self.chi_square_limits_only = chi_square_limits_only
        self.extrapolation = extrapolation
        self.molecular_absorption_lookup_table_path = molecular_absorption_lookup_table_path
        self.continue_from_file = continue_from_file
        self.object_photometry_file = object_photometry_file
        self.ntrials = ntrials
        self.nprocesses = nprocesses
        self.dusty_file_dir = dusty_file_dir
        self.emcee_runners = []
        self.dusty = dusty
        self.full_results_filename = ''

    @staticmethod
    def run_mcmc_process(emcee_runner: Emcee):
        emcee_runner.run_mcmc()

    def run(self):
        logger.info(f"Creating {self.nprocesses} processes.")
        pool = Pool(self.nprocesses)
        tstart = datetime.utcnow()
        for process_num in range(self.nprocesses):
            dusty_process_working_dir = self.working_dir + f'_{process_num}'
            basic_dusty = self.dusty(parameters=self.dusty_parameters,
                                dusty_working_directory=dusty_process_working_dir,
                                dusty_file_directory=self.dusty_file_dir
                                )
            emcee_runner = Emcee(dusty=basic_dusty,
                                 obsdata=self.obsdata,
                                 object_photometry_file=self.object_photometry_file,
                                 nwalkers=self.nwalkers,
                                 ntrials=int(self.ntrials / self.nprocesses),
                                 correct_for_molecular_absorption=self.correct_for_molecular_absorption,
                                 molecular_absorption_lookup_table_path=self.molecular_absorption_lookup_table_path,
                                 limits_only=self.limits_only,
                                 fixLstar=self.fixLstar,
                                 chi_square_limits_only=self.chi_square_limits_only,
                                 extrapolation=self.extrapolation,
                                 continue_from_file=self.continue_from_file,
                                 random_seed=process_num
                                 )

            self.emcee_runners.append(emcee_runner)

        pool.map(self.run_mcmc_process, self.emcee_runners)
        tend = datetime.utcnow()
        logger.info(f'Start {tstart}, End {tend}')

    def write_results_file(self):
        result = Table()
        for emcee_runner in self.emcee_runners:
            tmp = ascii.read(emcee_runner.outfilename)
            tmp['run'] = emcee_runner.random_seed
            result = vstack([result, tmp])

        full_results_basename = os.path.basename(self.object_photometry_file)+'_full_results.dat'
        self.full_results_filename = os.path.join(self.working_dir, full_results_basename)
        result.write(self.full_results_filename, format='csv', overwrite=True)

    def make_plots(self):
        self.write_results_file()
        chains = ascii.read(self.full_results_filename)
        #
        plotfile_basename = os.path.basename(self.object_photometry_file).replace('.dat','.pdf')
        plotfilename = os.path.join(self.working_dir, plotfile_basename)
        chains.remove_columns(['log_scaling', 'r1', 'log_posterior'])
        make_corner_plot(chains, savefilename=plotfilename)

