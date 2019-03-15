from abc import abstractmethod
from eppy.modeleditor import IDF


def equal(a, b):
    try:
        return a.upper() == b.upper()
    except AttributeError:
        return a == b


class UtilIdf:

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
    def get_output(self):
        """

        :return: model output of dtype=float
        """
        pass


class EppyUtilIdf(UtilIdf):

    def copy(self):

        idf_txt = self.idf.idfstr()
        try:
            idf_copy = IDF(epw=getattr(self.idf, 'epw'))
        except AttributeError:
            idf_copy = IDF()
        idf_copy.initreadtxt(idf_txt)
        return idf_copy

    def mod(self, obj_id, obj_name, field, value):

        idf = self.copy(self.idf) # copy file to prevent modifying same idf
        obj_lst = idf.idfobjects[obj_id]

        for obj in obj_lst:
            if not obj_name:
                obj[field] = value
                break
            elif obj_name == obj.Name or obj_name == 'ALL':
                obj[field] = value

        return idf

    def get_output(self, html_tbls, arg_lst):

        if not arg_lst:
            return None
        elif len(arg_lst) == 1:
            key = arg_lst.pop(0)
            for tbl in html_tbls:
                heading = tbl[0]
                if equal(heading, key):
                    return tbl[1:]
        else:
            key = arg_lst.pop(0)
            for tbl in html_tbls:
                heading = tbl[0]
                if equal(heading, key):
                    return self.get_output(tbl[1], arg_lst)
        return None
