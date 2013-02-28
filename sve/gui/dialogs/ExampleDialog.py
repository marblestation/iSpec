import numpy as np
from CustomDialog import *


class ExampleDialog(CustomDialog):
    def plot(self, axes, component):
        axes.plot(np.arange(10), np.arange(10), lw=1, color='b', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)

    def __init__(self, parent, title):
        self.__parent = parent
        self.__title = title
        self.__components = []
        self.__components = []
        component = {}
        component["type"] = "Label"
        component["text"] = "Suggested number of splines: 1"
        self.__components.append(component)
        component = {}
        component["type"] = "Plot"
        component["function"] = self.plot
        self.__components.append(component)
        component = {}
        component["type"] = "Listbox"
        component["options"] = ["Testing", "Hi!"]
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "First"
        component["text-type"] = "float" # float, int or str
        component["default"] = "1.0"
        component["minvalue"] = 1.0
        component["maxvalue"] = 2.0
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Second"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Third"
        component["options"] = ["a", "b", "c", "d"]
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Fourth"
        component["options"] = ["a", "b", "c"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show():
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


