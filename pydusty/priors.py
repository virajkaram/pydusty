import numpy as np
import logging

logger = logging.getLogger()

class Prior:

    @property
    def required_n_parameters(self):
        raise NotImplementedError

    @property
    def name(self):
        raise NotImplementedError

    def prior_limits(self) -> tuple:
        raise NotImplementedError

    def __init__(self,
                 prior_parameters):
        self.prior_parameters = prior_parameters

    def get_log_value(self, value):
        raise NotImplementedError

    def validate_prior_parameters(self):
        if not len(self.prior_parameters) == self.required_n_parameters:
            err = f'Invalid number of parameters for {self.name} prior'
            raise ValueError(err)

    def get_samples(self, nsamples):
        distsamples = self.return_distsamples()
        randindex = np.random.randint(len(distsamples), size=nsamples)
        random_samples = distsamples[randindex]
        return random_samples

    def return_distsamples(self):
        distsamples = []
        prior_extent = self.prior_limits()

        xranges = np.linspace(prior_extent[0], prior_extent[1], 100)
        prob_values = np.array([np.exp(self.get_log_value(x)) for x in xranges])
        min_prob = np.min(prob_values[prob_values > 0])
        power = np.floor(np.log10(min_prob))
        prob_values = prob_values/(10**power)

        for ind, i in enumerate(xranges):
            prob_value = prob_values[ind]
            distsamples.extend([i for j in range(int(1e3 * prob_value))])

        distsamples = np.array(distsamples)
        np.random.shuffle(distsamples)

        return distsamples


class GaussianPrior(Prior):
    required_n_parameters = 2
    name = 'gaussian'

    def __init__(self, *args, **kwargs):
        super(GaussianPrior, self).__init__(*args, **kwargs)
        self.validate_prior_parameters()
        self.mu, self.sigma = self.prior_parameters

    def prior_limits(self):
        return self.sigma - 5 * self.mu, self.sigma + 5 * self.mu

    def get_log_value(self, value):
        return -np.log(self.sigma) - ((value - self.mu) ** 2) / (2 * self.sigma ** 2)


class UniformPrior(Prior):
    required_n_parameters = 2
    name = 'uniform'

    def __init__(self, *args, **kwargs):
        super(UniformPrior, self).__init__(*args, **kwargs)
        self.validate_prior_parameters()
        self.minimum, self.maximum = self.prior_parameters

    def prior_limits(self):
        return self.minimum, self.maximum

    def get_log_value(self, value):
        if np.logical_and((self.minimum < value), (self.maximum > value)):
            return np.log(1 / (self.maximum - self.minimum))
        else:
            return -np.inf
