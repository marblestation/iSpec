import numpy as np
from CustomDialog import *

class EstimateErrorsDialog(CustomDialog):
    def __init__(self, parent, title, snr):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "SNR (Signal-to-Noise Ratio)"
        component["text-type"] = "float" # float, int or str
        component["default"] = snr
        component["minvalue"] = 1.0e-10
        component["maxvalue"] = np.inf
        self.__components.append(component)

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


