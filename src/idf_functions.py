from abc import abstractmethod

class UtilIdf():

    def __init__(self, idf):
        self.idf = idf


    @abstractmethod
    def copy(self):
        pass


    @abstractmethod
    def mod(self, obj_id, obj_name, field, value):
        """

        :param obj_id: object ID used to identify EnergyPlus object
        :param obj_name: name of idf object
        :param field: name of field to be modified
        :param value: value of idf obj to modify to
        :return:
        """
        pass


    @abstractmethod
    def get_output(self, arg):
        """

        :param arg: name of argument to extract model output
        :return: model output of dtype=float
        """
        pass
