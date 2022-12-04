from pydusty.priors import Prior


class Parameter:
    def __init__(self,
                 name: str,
                 value: float,
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


class DustyParameters:
    def __init__(self,
                 tstar: Parameter =None,
                 tdust: Parameter =None,
                 tau: Parameter =None,
                 shell_thickness: Parameter =None,
                 dust_type: Parameter =None,
                 blackbody: Parameter = None,
                 tstarmin: Parameter = None,
                 tstarmax: Parameter = None,
                 custom_grain_distribution: Parameter = None,
                 min_grain_size: Parameter = None,
                 max_grain_size: Parameter = None,
                 ebv: Parameter = None
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
        if self.ebv is None:
            self.ebv = Parameter(name='ebv', value=0, is_variable=False)

        self.parameter_dictionary = vars(self)
        self.parameter_list = [self.parameter_dictionary[x] for x in self.parameter_dictionary]

    def get_printable_string(self):
        return ''

    def __iter__(self):
        for x in self.parameter_list:
            yield x
