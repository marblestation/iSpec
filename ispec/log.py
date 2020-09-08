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
import logging
import logging.handlers

#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
LOG_FILE = "ispec.log"
CONSOLE = True

logger = logging.getLogger() # root logger, common for all
#logger = logging.getLogger(name)
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))

#formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(name)s [%(funcName)s:%(lineno)d]: %(message)s')
formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] [%(module)s:%(funcName)s:%(lineno)d]: %(message)s')
#formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

if CONSOLE:
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

megabyte = 1048576
try:
    handler = logging.handlers.RotatingFileHandler(LOG_FILE, 'a', maxBytes=50*megabyte, backupCount=5)
except IOError as e:
    logging.error("Logging information will not be stored in a file ({})".format(str(e)))
else:
    handler.setFormatter(formatter)
    logger.addHandler(handler)
# It is accessible via import logging; logging.warning("x")
