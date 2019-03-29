import os
import shutil
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
        self.problem['bounds'] = param['bounds']
        self.problem['num_vars'] = param['num_vars']
        self.objects = param['obj_id']
        for obj in self.objects:
            # TODO find longest string match in list instead of using obj[0]
            self.problem['names'].append(obj[0])
        '''
        for k, v in param.items():
            self.problem['names'].append(k)
            range = v['range']
            self.problem['bounds'].append(range)
        '''

        self.X = sample(self.problem, N, num_levels)
        self.y = None
        self.si = None

    @abstractmethod
    def set_idf(self, idf_path, epw_path, idd_version):
        pass

    def get_samples(self):
        """

        :return: a numpy.ndarray containing the model inputs
        required for Method of Morris. The resulting matrix has
        (G/D+1)*N/T rows and D columns, where D is the number of
        parameters.
        """
        return self.X

    def get_idf(self):
        return self.idf

    def get_names(self):
        """

        :return: a list containing the names of parameters
        that will be screened for their sensitivity
        """

        return self.problem['names']

    @abstractmethod
    def run_models(self, processors):
        """
        run E+ models using variations based on self.X values and
        set self.y values based on output of the simulations
        :param idf:  energyplus idf object
        :param processors: number of processors
        :return:
        """
        pass

    def evaluate(self):
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

    def run_models(self, processors):
        if self.idf is None:
            raise TypeError('idf not set')
        util = EppyUtilIdf()
        idf_lst = []
        file_dir = os.path.dirname(__file__)
        output_folder = os.path.join(file_dir, 'sensitivity_results')
        try:
            shutil.rmtree(output_folder)
        except FileNotFoundError as e:
            print(e)
        os.mkdir(output_folder)

        for i, values in enumerate(self.X):
            idf_temp = util.copy(self.idf)
            for j, value in enumerate(values):
                for obj in self.objects[j]:
                    obj_id, obj_name, field = obj.split(',')
                    util.mod(idf_temp, obj_id, obj_name, field, value)
            idf_temp.idfname = os.path.join(output_folder, 'run-{}.idf'.format(i))
            idf_temp.save()

            sim_settings = {'ep_version': '8-7-0',
                            'verbose': 'q',
                            'output_directory': output_folder,
                            'readvars': True,
                            'output_prefix': "run-{}-".format(i)}
            idf_lst.append([idf_temp, sim_settings])

        run_functions.runIDFs(idf_lst, processors=processors)
        # retrieve E+ outputs after simulations are run
        num_samples = self.X.shape[0]
        self.y = np.zeros(num_samples)

        for k in range(num_samples):
            result_file = os.path.join(output_folder, 'run-{}-tbl.htm'.format(k))
            with open(result_file, "r", encoding="ISO-8859-1") as f:
                results_table = readhtml.titletable(f.read())

            total_site_energy = util.get_output(results_table, ['Site and Source Energy', 'Total Site Energy'])
            total_site_energy_per_area = total_site_energy[1]
            self.y[k] = total_site_energy_per_area


folder = '../calibration_example'
idf_name = 'eplus_model/sensitivity_model.idf'
idd_name = 'idd/Energy+V8_7_0.idd'
epw_name = 'amy_weatherfile/2013_AMY.epw'

idf_path = os.path.abspath(os.path.join(folder, idf_name))
idd_path = os.path.abspath(os.path.join(folder, idd_name))
epw_path = os.path.abspath(os.path.join(folder, epw_name))
parameters = {'bounds': [[21, 25],
                         [15, 20.5],
                         [2, 20],
                         [3, 30],
                         [0.0001, 0.001524],
                         [1200, 1800],
                         [0.1, 0.9],
                         [0.86, 1.0],
                         [1, 6],
                         [0.8, 1.0]],
              'obj_id': [['Schedule:Day:Interval,All Weekends 12,Value_Until_Time_1',
                          'Schedule:Day:Interval,Space Thermostat Cooling Setpoint Default Weekday Schedule,Value_Until_Time_1'],
                         ['Schedule:Day:Interval,All Weekends 13,Value_Until_Time_1',
                          'Schedule:Day:Interval,Space Thermostat Heating Setpoint Default Weekday Schedule,Value_Until_Time_1'],
                         ['Lights,AllLights,Watts_per_Zone_Floor_Area'],
                         ['ELECTRICEQUIPMENT,AllPlugLoads,Watts_per_Zone_Floor_Area'],
                         ['ZoneInfiltration:DesignFlowRate,Utility Leaky Bldg Infiltration 1,Flow_per_Exterior_Surface_Area'],
                         ['ZoneHVAC:Baseboard:Convective:Electric,Elec Baseboard,Heating_Design_Capacity'],
                         ['Fan:ConstantVolume,User defined fan 1,Fan_Total_Efficiency'],
                         ['Fan:ConstantVolume,User defined fan 1,Motor_Efficiency'],
                         ['Coil:Cooling:DX:SingleSpeed,Coil Cooling DX Single Speed 1,Gross_Rated_Cooling_COP'],
                         ['Coil:Heating:Fuel,OS:Coil:Heating:Gas 1,Burner_Efficiency']],
              'num_vars': 10}

IDF.setiddname(idd_path)
base_idf = IDF(idf_path, epw_path)

morris = MorrisEppy(param=parameters, N=4)
morris.set_idf(idf_path, epw_path, idd_path)

morris.run_models(processors=4)
si = morris.evaluate()
