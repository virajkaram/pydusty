from pydusty.utils import extinction_correct_obsdata
import numpy as np
import logging
from pydusty.utils import calculate_molecular_absorption_fractions
import emcee
import os
from pydusty.dusty import Dusty
from tqdm import tqdm

logger = logging.getLogger(__name__)


class Emcee:
    def __init__(self,
                 nwalkers: int,
                 dusty: Dusty,
                 obsdata: dict,
                 object_photometry_file: str,
                 ntrials: int = 250,
                 correct_for_molecular_absorption: bool =False,
                 molecular_absorption_lookup_table_path: str = None,
                 limits_only: bool = False,
                 fixLstar: bool = False,
                 chi_square_limits_only: float = 4,
                 extrapolation: bool = False,
                 continue_from_file: bool = False,
                 random_seed:int =0,
                 ):
        self.nwalkers = nwalkers
        self.dusty = dusty
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
        self.random_seed = random_seed

        output_file_basename = os.path.basename(self.object_photometry_file) + '_results.dat'
        self.outfilename = f'{self.dusty.dusty_working_directory}/{output_file_basename}'
        self.variable_parameters = [x for x in self.dusty.parameters if x.is_variable]
        self.constant_parameters = [x for x in self.dusty.parameters if not x.is_variable]
        self.ndim = len(self.variable_parameters)

    @staticmethod
    def get_luminosity_log_chisq(data, abs_fracs=None):
        pass

    @staticmethod
    def get_interpolated_model_fluxes(observed_wav, model_wav, model_flux):
        interpolated_model_fluxes = np.interp(observed_wav, model_wav, model_flux)
        return interpolated_model_fluxes

    @staticmethod
    def correct_fluxes_for_molecular_absorption(fluxes, filter_names, tstar, absorption_lookup_table_path):
        absorption_fractions = calculate_molecular_absorption_fractions(tstar, filter_names,
                                                                        absorption_lookup_table_path)
        assert len(absorption_fractions) == len(fluxes)
        fluxes *= absorption_fractions
        return fluxes

    @staticmethod
    def get_log_luminosity_scaling(observed_fluxes,
                                   observed_fluxerrs,
                                   model_fluxes,
                                   limits_only=False,
                                   chi_square_limits_only: float = 4,
                                   fixLstar: float = None
                                   ):
        '''
        calculate a weighted scaling
        '''
        if fixLstar:
            log_weighted_luminosity_scaling = fixLstar
        elif limits_only:
            if np.any(observed_fluxes>0):
                err = 'limits_only option is set, expecting only upper limits (i.e. negative flux entries)'
                raise ValueError(err)
            scalings_array = (model_fluxes/observed_fluxerrs)**2
            weighted_scaling = chi_square_limits_only/(np.sum(scalings_array))
            log_weighted_luminosity_scaling = np.log10(weighted_scaling)
        else:
            scale = (np.sum(observed_fluxes*model_fluxes/(observed_fluxerrs**2))
                     /np.sum(model_fluxes**2/(observed_fluxerrs**2)))

            log_weighted_luminosity_scaling = np.log10(scale)

        return log_weighted_luminosity_scaling

    @staticmethod
    def get_log_scaled_inner_radius(r1, log_luminosity_scaling):
        log_scaled_radius = np.log10(r1) + 0.5 * (log_luminosity_scaling - 4.0)
        return log_scaled_radius

    @staticmethod
    def get_log_outer_radius(log_inner_radius, thickness):
        log_outer_radius = log_inner_radius + np.log10(thickness)
        return log_outer_radius

    @staticmethod
    def calculate_luminosity_chi2(observed_fluxes, observed_fluxerrs,
                                  scaled_model_fluxes,  limits_only=False,
                                  extrapolation=False,
                                  error_underestimation_scaling_factor=0,):
        chi_Lobs = 0
        if not extrapolation and not limits_only:
            chi_array = np.zeros(len(observed_fluxes))

            limit_mask = (observed_fluxes < 0)
            detection_mask = np.invert(limit_mask)

            # log_observed_fluxes = np.log10(observed_fluxes)
            # log_observed_fluxerrs = observed_fluxerrs / (observed_fluxes * np.log(10))
            # log_model_fluxes = np.log10(scaled_model_fluxes)

            chi_array[detection_mask] = (((observed_fluxes[detection_mask]
                                          - scaled_model_fluxes[detection_mask]))**2
                                         /(observed_fluxerrs[detection_mask]
                                           + error_underestimation_scaling_factor
                                           * observed_fluxes[detection_mask])**2)
            chi_array[limit_mask] = (scaled_model_fluxes[limit_mask]
                                     /observed_fluxerrs[limit_mask])**2

            chi_Lobs = np.sum(chi_array)
            eweight = 1.0

            chi_Lobs /= eweight**2
        return -0.5*chi_Lobs

    def calculate_chi2(self, lam, flx, parameter_dict):
        if self.dusty.parameters.ebv.is_variable:
            corr_data = extinction_correct_obsdata(self.obsdata, parameter_dict['ebv'])
        else:
            corr_data = self.obsdata

        mlam = corr_data['mlam']
        mlum = corr_data['mlum']
        merr = corr_data['merr']
        filter_names = corr_data['filters']

        interpolated_model_fluxes = self.get_interpolated_model_fluxes(mlam, lam, flx)

        if self.correct_for_molecular_absorption:
            interpolated_model_fluxes = self.correct_fluxes_for_molecular_absorption(interpolated_model_fluxes,
                                                                                     filter_names,
                                                                                     parameter_dict['tstar'],
                                                                                     self.molecular_absorption_lookup_table_path)

        log_luminosity_scaling = self.get_log_luminosity_scaling(observed_fluxes=mlum,
                                                                 observed_fluxerrs=merr,
                                                                 model_fluxes=interpolated_model_fluxes,
                                                                 limits_only=self.limits_only,
                                                                 chi_square_limits_only=self.chi_square_limits_only,
                                                                 fixLstar=self.fixLstar)

        scaled_model_fluxes = interpolated_model_fluxes * 10**log_luminosity_scaling
        chi2_lum = self.calculate_luminosity_chi2(
            observed_fluxes=mlum,
            observed_fluxerrs=merr,
            scaled_model_fluxes=scaled_model_fluxes,
            limits_only=self.limits_only,
            extrapolation=self.extrapolation,
            error_underestimation_scaling_factor=
            self.dusty.parameters.error_underestimate_factor.value
        )

        return chi2_lum, log_luminosity_scaling

    def calculate_chi2_manually(self, lam, flx, r1, par_dict):
        if self.dusty.parameters.ebv.is_variable:
            corr_data = extinction_correct_obsdata(self.obsdata, par_dict['ebv'])
        else:
            corr_data = self.obsdata

        mlam = corr_data['mlam']
        mlum = corr_data['mlum']
        merr = corr_data['merr']
        mluml = corr_data['mluml']
        merrl = corr_data['merrl']
        filters_names = corr_data['filters']
        nm = len(mlum)
        thick = par_dict['thick']

        abs_fracs = []
        if self.correct_for_molecular_absorption:
            abs_fracs = calculate_molecular_absorption_fractions(par_dict['tstar'], filters_names,
                                                                 self.molecular_absorption_lookup_table_path)

        aa = 0.0
        bb = 0.0
        values = [] # np.zeros(100)
        log_values = [] # np.zeros(100)

        for i in range(nm):
            # j = np.argmin(np.abs(lam - mlam[i]))
            # if lam[j] > mlam[i]:
            #     j = j - 1
            # j = locate(lam,npt,mlam[i]) - 1    # -1 because indexed to 0
            # val = flx[j] + (flx[j + 1] - flx[j]) * (mlam[i] - lam[j]) / (lam[j + 1] - lam[j])
            val = np.interp(x=mlam[i], xp=lam, fp=flx)
            if self.correct_for_molecular_absorption and not self.limits_only:
                val = val * abs_fracs[i]

            # logger.debug(flx[j], flx[j + 1], lam[j], lam[j + 1], mlam[i], val[i], mlum[i])

            if self.limits_only:
                aa = aa + (val / merr[i]) ** 2
                if mlum[i] > 0:
                    err = 'ERROR: I am expecting only upper limits'
                    raise ValueError(err)

            else:
                logval = np.log10(val + 1.e-32)
                logger.debug(f'{i} logval = {logval}')
                if mlum[i] > 0:
                    aa = aa + (mluml[i] - logval) / merrl[i] ** 2
                    bb = bb + 1.0 / merrl[i] ** 2

            values.append(val)
            log_values.append(logval)

        if self.limits_only:
            slum = np.sqrt(self.chi_square_limits_only / aa)
            sluml = np.log10(slum)
        else:
            sluml = aa / bb
        # if extrapolation:
        #    sluml = llum  # if extrapolating tau from a best-fit model then use the same luminosity
        if self.fixLstar:
            logger.debug(f"forcing Lstar to 10^{self.fixLstar}")
            sluml = self.fixLstar
        # r1 output from DUSTY is the distance at which a point source w/luminosity 10^4 L_sun produces the
        # bolometric flux F_e1
        r1 = np.log10(r1) + 0.5 * (sluml - 4.0)  # scale radius for luminosity
        r2 = r1 + np.log10(thick)
        # dusty reports a radius scaled to a luminosity of 10^4
        # so this scales it to the luminosity you just worked out
        chi = 0.0
        chis = {}
        if not self.extrapolation and not self.limits_only:
            for i in range(nm):
                if mlum[i] > 0:
                    chis[filters_names[i]] = ((mluml[i] - sluml - logval[i]) / merrl[i]) ** 2
                    chi = chi + chis[filters_names[i]]
            # chi = chi + ((mluml[i]-sluml-vall[i])/merrl[i])**2
                logger.debug('%s = chi + ((%s-%s-%s)/%s)**2' % (chi, mluml[i], sluml, logval[i], merrl[i]))
                logger.debug('sluml = %s' % (sluml))
                if mlum[i] < 0:
                    rat = 10.0 ** (sluml + logval[i] - merrl[i])
                    logger.debug('rat = %s = 10**(%s+%s-%s)' % (rat, sluml, logval[i], merrl[i]))
                    chis[filters_names[i]] = rat * rat
                    chi = chi + chis[filters_names[i]]
            # chi = chi + rat*rat
        chi_Lobs = chi
        logger.info('chi^2 from L_obs = %0.1f' % (chi_Lobs))
        eweight = 1.0
        chi = chi / eweight ** 2

        return -0.5 * chi, sluml

    def log_posterior(self, varx):
        parameter_dict = {}
        for parameter in self.dusty.parameters.parameter_list:
            parameter_dict[parameter.name] = parameter.value

        fileuse = ''
        if not self.dusty.parameters.blackbody:
            err = "The non-blackbody mode has not been implemented yet. Please set blackbody to True."
            # Need to add stellar spectra here, and write update fileuse
            raise AttributeError(err)
        parameter_dict['fileuse'] = fileuse

        lgpri = 0
        for idx in range(len(varx)):
            self.variable_parameters[idx].update_value(varx[idx])
            lgpri += self.variable_parameters[idx].prior.get_log_value(varx[idx])

        if lgpri == -np.inf:
            return lgpri, lgpri, lgpri

        self.dusty.generate_input()

        self.dusty.run()

        lam, flx, npt, r1, ierror = self.dusty.get_results()

        chi2, sluml = self.calculate_chi2(lam, flx, parameter_dict)
        log_r1 = self.get_log_scaled_inner_radius(r1=r1, log_luminosity_scaling=sluml)
        log_post = chi2 + lgpri

        logger.info(f'log_chi2: {chi2}')
        logger.info(f'log_post: {log_post}')
        return log_post, sluml, log_r1

    def run_mcmc(self):
        curdir = os.getcwd()
        dusty_inpdir = self.dusty.dusty_working_directory
        logger.info(f'Changing directory to {dusty_inpdir}')
        os.chdir(dusty_inpdir)

        seedind = self.random_seed ** 2
        np.random.seed(seedind)
        logger.info(f'Seed for np is {seedind}')
        # state = np.random.RandomState(seedind)

        variable_parameter_names = [x.name for x in self.variable_parameters]
        variable_parameter_pritypes = [x.prior.name for x in self.variable_parameters]
        variable_parameter_priparams = [x.prior.prior_parameters for x in self.variable_parameters]
        variable_parameter_initpos = [x.value for x in self.variable_parameters]
        variable_parameter_initpos_nwalkers = np.array([x.prior.get_samples(nsamples=self.nwalkers) for x in self.variable_parameters]).T

        if not self.continue_from_file:
            f = open(self.outfilename, "w")
            f.write(f'#Initial points\n')
            # for initpos in variable_parameter_initpos_nwalkers:
            #     f.write(f'# {initpos}\n')
            f.write(f'#Variables: {variable_parameter_names}\n')
            f.write(f'#Priors: {variable_parameter_pritypes}\n')
            f.write(f'#Prior params: {variable_parameter_priparams}\n')
            for varname in variable_parameter_names:
                f.write(f"{varname}\t")
            f.write(f"log_scaling\t")
            f.write(f"r1\t")
            f.write(f"log_posterior\n")
            f.close()

        logger.info(
            f'Running emcee by varying parameters {variable_parameter_names} with initial values '
            f'{variable_parameter_initpos} and priors {variable_parameter_priparams}')

        dtype = [("sluml", float), ("r1", float)]

        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.log_posterior, blobs_dtype=dtype)
        # logger.info(sampler.random_state)
        # sampler.random_state.setter(state.get_state)

        for result in tqdm(sampler.sample(variable_parameter_initpos_nwalkers, iterations=self.ntrials, progress=True, store=True)):
            print(os.getcwd())
            position = result.coords
            lp = result.log_prob
            blobs = result.blobs
            sluml = blobs['sluml']
            r1s = blobs['r1']
            f = open(self.outfilename, "a")
            for k in range(len(position)):
                for j in position[k]:
                    f.write('%.4f\t' % (j))
                f.write('%.4f\t' % (sluml[k]))
                f.write('%.4f\t' % (r1s[k]))
                f.write('%.4f\n' % (lp[k]))
            f.close()
        os.chdir(curdir)