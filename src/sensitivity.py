from SALib.sample.morris import sample
from SALib.analyze.morris import analyze


class MorrisBem():

    def __init__(self, param, N=12, num_levels=4):
        """

        :param dict param: Uncertain parameters for the sensitivity analysis
        :param int N: The number of trajectories to generate (recommended 10-20)
        :param int num_levels: The number of grid levels (default 4)
        """

        self.num_levels = num_levels
        self.problem ={'num_vars': len(param), 'names': [], 'bounds': []}
        for k, v in param.items():
            self.problem['names'].append(k)
            self.problem['bounds'].append(v)

        self.X = sample(self.problem, N, num_levels)
        self.si = None

    def get_samples(self):
        """

        :return: a numpy.ndarray containing the model inputs
        required for Method of Morris. The resulting matrix has
        (G/D+1)*N/T rows and D columns, where D is the number of
        parameters.
        """
        return self.X

    def get_names(self):
        """

        :return: a list containing the names of parameters
        that will be screened for their sensitivity
        """

        return self.problem['names']

    def evaluate(self, y):
        """

        :param y: A Numpy array containing the model outputs of dtype=float
        :return: A dictionary of sensitivity indices containing the following entries.
            - `mu` - the mean elementary effect
            - `mu_star` - the absolute of the mean elementary effect
            - `sigma` - the standard deviation of the elementary effect
            - `mu_star_conf` - the bootstrapped confidence interval
            - `names` - the names of the parameters
        """

        self.si = analyze(self.problem, self.X, y,
                          num_levels=self.num_levels)
        return self.si
