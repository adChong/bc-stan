from abc import abstractmethod
from eppy.modeleditor import IDF


def equal(a, b):
    """

    :param a: string a
    :param b: string b
    :return: check if string a equals string b ignoring character case
    """
    try:
        return a.upper() == b.upper()
    except AttributeError:
        return a == b


class UtilIdf:

    @abstractmethod
    def copy(self, idf):
        """

        :param idf: EnergyPlus idf to be copied
        :return: a copy of the EnergyPlus idf
        """
        pass

    @abstractmethod
    def mod(self, idf, obj_id, obj_name, field, value):
        """
        modifies an EnergyPlus IDF in place
        :param obj_id: object ID used to identify EnergyPlus object
        :param obj_name: name of idf object
        :param field: field to be modified
        :param value: value of idf obj to modify to
        """
        pass

    @abstractmethod
    def get_output(self):
        """

        :return: model output of dtype=float that would be used
                for the sensitivity analysis
        """
        pass


class EppyUtilIdf(UtilIdf):

    def copy(self, idf):

        idf_txt = idf.idfstr()
        try:
            idf_copy = IDF(epw=getattr(idf, 'epw'))
        except AttributeError:
            idf_copy = IDF()
        idf_copy.initreadtxt(idf_txt)
        return idf_copy

    def mod(self, idf, obj_id, obj_name, field, value):

        # modification done in place
        obj_lst = idf.idfobjects[obj_id.upper()]

        for obj in obj_lst:
            if not obj_name:
                obj[field] = value
                break
            elif obj_name == obj.Name:
                obj[field] = value

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
