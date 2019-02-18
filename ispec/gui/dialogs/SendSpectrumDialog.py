#
#    This file is part of iSpec.
#    Copyright Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
from CustomDialog import *

class SendSpectrumDialog(CustomDialog):

    def __init__(self, parent, title, applications):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Application"
        component["options"] = applications
        component["default"] = applications[0]
        self.__components.append(component)

    def show(self, updated_applications=None):
        self.results = None
        if updated_applications is not None:
            self.__components[0]["options"] = updated_applications
            # Validate that the default value (previous user selected value) exists in the new application list
            default_ok = False
            for application in updated_applications:
                if application == self.__components[5]["default"]:
                    default_ok = True
            if not default_ok:
                self.__components[5]["default"] = updated_applications[0]
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)

