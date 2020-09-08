#!/usr/bin/env python
#
#    This file is part of the Integrated Spectroscopic Framework (iSpec).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
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
from builtins import str
import os
import numpy as np
from subprocess import call
import tempfile
import shutil
from astropy.io import ascii

def interpolate_isochrone(dirname, age, feh, alpha_over_iron=0.0):
    YY_dir = dirname + "/isochrones/"
    YY_grid = dirname + "/input/isochrones/"

    if not os.path.exists(YY_dir + "/YYmix2"):
        raise Exception("YONSEI-YALE ISOCHRONES programmes not found in '%s'" % YY_dir)
    if age < 0.001 or age > 20.:
        raise Exception("Age out of range!")
    if feh < -3.0 or feh > 0.75:
        raise Exception("Metallicity out of range!")
    if alpha_over_iron < 0 or alpha_over_iron > 0.6:
        raise Exception("Alpha enhancement out of range!")
    tmp_dir = tempfile.mkdtemp()
    os.symlink(YY_grid, tmp_dir+"/i")

    #### YY.age file can contain one age per line
    age_file = open(tmp_dir + "/YY.age", "w")
    age_file.write(str(age))
    age_file.close()

    targetZ = -0.01 # When it is negative it is ignored and FeH is considered
    content = """ $INPUT
 AFe=%.2f
 targetZ=%.2f
 FeH=%.2f
 Age= %i
 YYout='output.txt'
 $END
 $ISET1
 NYYiso=11
 YYiso(1)='i/yy00l.x76997z00001a0o2v2'
 YYiso(2)='i/yy00l.x7697z0001a0o2v2'
 YYiso(3)='i/yy00l.x7688z0004a0o2v2'
 YYiso(4)='i/yy00l.x767z001a0o2v2'
 YYiso(5)='i/yy00l.x758z004a0o2v2'
 YYiso(6)='i/yy00l.x749z007a0o2v2'
 YYiso(7)='i/yy00l.x74z01a0o2v2'
 YYiso(8)='i/yy00l.x71z02a0o2v2'
 YYiso(9)='i/yy00l.x65z04a0o2v2'
 YYiso(10)='i/yy00l.x59z06a0o2v2'
 YYiso(11)='i/yy00l.x53z08a0o2v2'
 $END
 $ISET2
 NYYiso=11
 YYiso(1)='i/yy00l.x76997z00001a2o2v2'
 YYiso(2)='i/yy00l.x7697z0001a2o2v2'
 YYiso(3)='i/yy00l.x7688z0004a2o2v2'
 YYiso(4)='i/yy00l.x767z001a2o2v2'
 YYiso(5)='i/yy00l.x758z004a2o2v2'
 YYiso(6)='i/yy00l.x749z007a2o2v2'
 YYiso(7)='i/yy00l.x74z01a2o2v2'
 YYiso(8)='i/yy00l.x71z02a2o2v2'
 YYiso(9)='i/yy00l.x65z04a2o2v2'
 YYiso(10)='i/yy00l.x59z06a2o2v2'
 YYiso(11)='i/yy00l.x53z08a2o2v2'
 $END
 $ISET4
 NYYiso=11
 YYiso(1)='i/yy00l.x76997z00001a4o2v2'
 YYiso(2)='i/yy00l.x7697z0001a4o2v2'
 YYiso(3)='i/yy00l.x7688z0004a4o2v2'
 YYiso(4)='i/yy00l.x767z001a4o2v2'
 YYiso(5)='i/yy00l.x758z004a4o2v2'
 YYiso(6)='i/yy00l.x749z007a4o2v2'
 YYiso(7)='i/yy00l.x74z01a4o2v2'
 YYiso(8)='i/yy00l.x71z02a4o2v2'
 YYiso(9)='i/yy00l.x65z04a4o2v2'
 YYiso(10)='i/yy00l.x59z06a4o2v2'
 YYiso(11)='i/yy00l.x53z08a4o2v2'
 $END
""" % (alpha_over_iron, targetZ, feh, age)

    params_file = open(tmp_dir + "/YY.nml", "w")
    params_file.write(content)
    params_file.close()

    previous_cwd = os.getcwd()
    os.chdir(tmp_dir)
    call([YY_dir + "/YYmix2" ])
    os.chdir(previous_cwd)

    #names = ["EEP", "M/Mo", "LogTeff", "LogG", "LogL/Lo", "U", "B", "V", "R", "I", "J", "H", "Ks", "Kp", "D51"]
    #isochrone = ascii.read(tmp_dir + "/tmp.iso", names=names)._data
    names = ["M/Msun", "logT", "logL/Ls", "logg", "Mv", "U-B", "B-V", "V-R", "V-I", "V-J", "V-H", "V-K", "V-L", "V-M", "#(x=-1)", "#(x=1.35)", "#(x=3)"]
    isochrone = ascii.read(tmp_dir + "/output.txt", data_start=3, names=names).as_array()
    shutil.rmtree(tmp_dir)

    return isochrone



if __name__ == '__main__':
    ################################################################################
    #--- iSpec directory -------------------------------------------------------------
    if os.path.exists('/home/sblancoc/shared/iSpec/'):
        # avakas
        ispec_dir = '/home/sblancoc/shared/iSpec/'
    elif os.path.exists('/home/blanco/shared/iSpec/'):
        # vanoise
        ispec_dir = '/home/blanco/shared/iSpec/'
    else:
        ispec_dir = '/home/marble/shared/iSpec/'
    ################################################################################
    #age = 9.5
    #feh = -0.5
    #isochrone = interpolate_isochrone(ispec_dir, age, feh)

    import matplotlib.pyplot as plt
    ages = [ 1.  ,   1.25,   1.5 ,   1.75,   2.  ,   2.25,   2.5 ,   2.75,
             3.  ,   3.25,   3.5 ,   3.75,   4.  ,   4.25,   4.5 ,   4.75,
             5.  ,   5.5 ,   6.  ,   6.5 ,   7.  ,   7.5 ,   8.  ,   8.5 ,
             9.  ,   9.5 ,  10.  ,  10.5 ,  11.  ,  11.5 ,  12.  ,  12.5 ,
            13.  ,  13.5 ,  14.  ,  14.5 ,  15.  ]
    feh = 0.0
    #plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim((4000, 8000))
    ax.set_ylim((1, 5))
    for age in ages:
        isochrone = interpolate_isochrone(ispec_dir, age, feh)
        ax.plot(np.power(10, isochrone['logT']), isochrone['logg'], marker='', ls='-', label=str(age) + " Gyrs")
    ax.set_xlabel("$T_{eff}$ (K)")
    ax.set_ylabel("$log(g)$ (dex)")
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(loc='best',prop={'size':12}, ncol=3)
    ax.grid(True)
    plt.show()

    import pudb
    pudb.set_trace()
