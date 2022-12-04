import numpy as np


class Prior:
    @property
    def prior_extent(self) -> tuple:
        raise NotImplementedError

    @property
    def required_n_parameters(self):
        raise NotImplementedError

    @property
    def name(self):
        raise NotImplementedError

    def __init__(self,
                 name,
                 prior_parameters,
                 required_n_parameters):
        self.prior_parameters = prior_parameters

    def get_log_value(self, value):
        raise NotImplementedError

    def validate_prior_parameters(self):
        if not len(self.prior_parameters) == self.required_n_parameters:
            err = f'Invalid number of parameters for {self.name} prior'
            raise ValueError(err)

    def get_samples(self, nsamples):
        distsamples = self.return_distsamples()
        randindex = np.random.randint(len(distsamples), nsamples)
        random_samples = distsamples[randindex]
        return random_samples

    def return_distsamples(self):
        distsamples = []
        for i in np.linspace(self.prior_extent[0], self.prior_extent[1], 100):
            prob_value = 10**self.get_log_value(i)
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

    def prior_extent(self):
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

    def prior_extent(self):
        return self.minimum, self.maximum

    def get_log_value(self, value):
        if self.minimum < value < self.maximum:
            return np.log(1 / (self.maximum - self.minimum))
        else:
            return -np.inf
