import os
import numpy as np
from SALib.sample.morris import sample
from SALib.analyze.morris import analyze
from abc import abstractmethod
from idf_functions import EppyUtilIdf
from eppy.runner import run_functions
from eppy.modeleditor import IDF
from eppy.results import readhtml


class MorrisBem:

    def __init__(self, param, N=12, num_levels=4):
        """

        :param dict param: Uncertain parameters for the sensitivity analysis
            {"parameter name": [lower bound, upper bound]}
        :param int N: The number of trajectories to generate (recommended 10-20)
        :param int num_levels: The number of grid levels (default 4)
        """
        self.idf = None;
        self.num_levels = num_levels
        self.problem ={'num_vars': len(param), 'names': [], 'bounds': []}
        for k, v in param.items():
            self.problem['names'].append(k)
            self.problem['bounds'].append(v)

        self.X = sample(self.problem, N, num_levels)
        self.y = None
        self.si = None
        self.objects = [name.split(',') for name in self.problem['names']]

    @abstractmethod
    def set_idf(self, idf_path, epw_path, idd_path):
        pass

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

    @abstractmethod
    def run_models(self, idf, processors):
        """
        run E+ models using variations based on self.X values and
        set self.y values based on output of the simulations
        :param idf:  energyplus idf object
        :param processors: number of processors
        :return:
        """
        pass

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
        if self.y is None:
            raise NameError('y is not defined. run e+ models first')
        self.si = analyze(self.problem, self.X, self.y,
                          num_levels=self.num_levels)
        return self.si


class MorrisEppy(MorrisBem):

    def set_idf(self, idf_path, epw_path, idd_path):
        IDF.setiddname(idd_path)
        self.idf = IDF(idf_path, epw_path)

    def run_models(self, idf, processors):
        if self.idf is None:
            raise TypeError('idf not set')
        util = EppyUtilIdf(self.idf)
        idf_lst = []
        file_dir = os.path.dirname(__file__)
        output_folder = os.path.join(file_dir, 'results')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        for i, values in enumerate(self.X):
            for j, value in enumerate(values):
                obj_id = self.objects[j][0]
                obj_name = self.objects[j][1]
                field = self.objects[j][2]
                idf = util.mod(obj_id, obj_name, field, value)
            idf.idfname = os.path.join(output_folder, 'run-{}.idf'.format(i))
            idf.save()

            parameters = {'ep_version': '8-6-0',
                          'verbose': 'q',
                          'output_directory': output_folder,
                          'readvars': True,
                          'output_prefix': "run-{}-".format(i)}
            idf_lst.append([idf, parameters])

        run_functions.runIDFs(idf_lst, processors=processors)
        # retrieve E+ outputs after simulations are run
        num_samples = self.X.shape[0]
        self.y = np.zeros(num_samples)

        for k in range(num_samples):
            result_file = os.path.join(output_folder, 'run-{}-tbl.htm'.format(k))

            with open(result_file, "r", encoding="ISO-8859-1") as f:
                results_table = readhtml.titletable(f.read())

            total_site_energy = EppyUtilIdf.get_output(results_table, ['Site and Source Energy', 'Total Site Energy'])
            total_site_energy_per_area = total_site_energy[1]
            self.y[k] = total_site_energy_per_area

