from __future__ import annotations

from pydusty.priors import Prior
import numpy as np
import logging

logger = logging.getLogger(__name__)

class Parameter:
    def __init__(self,
                 name: str,
                 value: float | bool | str,
                 is_variable: bool=False,
                 prior: Prior = None):
        self.name = name
        self.value = value
        self.is_variable = is_variable
        self.prior = prior

    def update_value(self, new_value):
        self.value = new_value

    def set_as_variable(self):
        self.is_variable = True

    def validate_parameter(self):
        if self.is_variable:
            if self.prior is None:
                err = f'Please specify a prior for variable parameter {self.name}'
                raise ValueError(err)

    def get_log_prior_value(self):
        if not self.is_variable:
            err = 'Querying prior value for a parameter that is not being varied is not permitted.'
            raise ValueError(err)

        else:
            log_prior_value = self.prior.get_log_value(self.value)
            return log_prior_value


false_boolean_parameter = Parameter(name='generic_false_param', value=False, is_variable=False)
true_boolean_parameter = Parameter(name='generic_true_param', value=True, is_variable=False)
default_tstarmin_parameter = Parameter(name='tstarmin', value=3500, is_variable=False)
default_tstarmax_parameter = Parameter(name='tstarmax', value=48999, is_variable=False)
default_tau_wavelength_micron_parameter = Parameter(name='tau_wav', value=0.55,
                                                    is_variable=False)
default_error_underestimate_factor = Parameter(name='error_underestimate_factor', value=0.0,
                                               is_variable=False)
class DustyParameters:
    def __init__(self,
                 tdust: Parameter,
                 tau: Parameter,
                 shell_thickness: Parameter,
                 dust_type: Parameter,
                 blackbody: Parameter,
                 tstar: Parameter = None,
                 tstarmin: Parameter = default_tstarmin_parameter,
                 tstarmax: Parameter = default_tstarmax_parameter,
                 custom_grain_distribution: Parameter = false_boolean_parameter,
                 tau_wavelength_microns: Parameter = default_tau_wavelength_micron_parameter,
                 min_grain_size: Parameter = None,
                 max_grain_size: Parameter = None,
                 ebv: Parameter = None,
                 si_dl_abundance: Parameter = None,
                 al_com_abundance: Parameter = None,
                 error_underestimate_factor: Parameter = default_error_underestimate_factor,
                 dust_composition_elements: list[str] = None,
                 dust_composition_abundances: list[float] = None,
                 custom_input_spectrum_file: Parameter = None,
                 ):
        self.tstar = tstar
        self.tdust = tdust
        self.tau = tau
        self.shell_thickness = shell_thickness
        self.dust_type = dust_type
        self.blackbody = blackbody
        self.tstarmin = tstarmin
        self.tstarmax = tstarmax
        self.custom_grain_distribution = custom_grain_distribution
        self.min_grain_size = min_grain_size
        self.max_grain_size = max_grain_size
        self.ebv = ebv
        self.tau_wavelength_microns = tau_wavelength_microns
        self.si_dl_abundance = si_dl_abundance
        self.al_com_abundance = al_com_abundance
        self.error_underestimate_factor = error_underestimate_factor
        self.custom_input_spectrum_file = custom_input_spectrum_file
        if self.ebv is None:
            self.ebv = Parameter(name='ebv', value=0, is_variable=False)

        if not self.custom_grain_distribution.value:
            self.min_grain_size = Parameter(name='min_grain_size', value=0, is_variable=False)
            self.max_grain_size = Parameter(name='max_grain_size', value=0, is_variable=False)

        else:
            if np.logical_or(self.min_grain_size is None, self.max_grain_size is None):
                err = 'custom grain size requested, but no min_grain_size or max_grain_size specified.'
                raise AttributeError(err)
        self.parameter_dictionary = vars(self)
        self.parameter_list = [self.parameter_dictionary[x] for x in self.parameter_dictionary.keys() if isinstance(self.parameter_dictionary[x], Parameter)]

        self.abundances_dict = {}
        if dust_composition_abundances is not None:
            if dust_composition_elements is None:
                err = 'dust_composition_abundances specified without dust_composition_elements'
                raise ValueError(err)
            if len(dust_composition_elements) != len(dust_composition_abundances):
                err = 'dust_composition_elements and dust_composition_abundances must be the same length'
                raise ValueError(err)
            for i in range(len(dust_composition_elements)):
                self.abundances_dict[dust_composition_elements[i]] = dust_composition_abundances[i]
        # logger.debug(self.parameter_dictionary)
        # logger.debug(self.parameter_list)

    def get_printable_string(self):
        info_str = ''
        for x in self.parameter_dictionary:
            if isinstance(self.parameter_dictionary[x], Parameter):
                info_str += f'{x}:{self.parameter_dictionary[x].value}\n'
        return info_str

    def __iter__(self):
        for x in self.parameter_list:
            yield x
