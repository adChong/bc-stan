from SALib.sample.morris import sample
from SALib.analyze.morris import analyze
from abc import abstractmethod
from idf_functions import EppyUtilIdf

class MorrisBem:

    def __init__(self, idf, param, N=12, num_levels=4):
        """
        :param idf: energyplus input data file
        :param dict param: Uncertain parameters for the sensitivity analysis
        :param int N: The number of trajectories to generate (recommended 10-20)
        :param int num_levels: The number of grid levels (default 4)
        """
        self.idf = idf
        self.num_levels = num_levels
        self.problem ={'num_vars': len(param), 'names': [], 'bounds': []}
        for k, v in param.items():
            self.problem['names'].append(k)
            self.problem['bounds'].append(v)

        self.X = sample(self.problem, N, num_levels)
        self.si = None
        self.names = [name.split(',') for name in self.problem['names']]

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

    @abstractmethod
    def run_models(self):
        pass





class MorrisEppy(MorrisBem):

    def run_models(self):
        util = EppyUtilIdf(self.idf)

        for i, values in enumerate(self.X):
            for j, value in enumerate(values):
                obj_id = self.names[j][0]
                obj_name = self.names[j][1]
                field = self.names[j][2]
                idf = util.mod(obj_id, obj_name, field, value)

        pass

