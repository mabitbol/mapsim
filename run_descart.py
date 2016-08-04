import numpy as np
import glob
import subprocess as sb
import os

class RunDescart():

    def run(self):
        fdir = '/cache/mabitbol/ebex/'
        datadir = fdir + 'bolo250/bolos/'
        bolodirs = glob.glob(datadir+'*')
        copyscript = fdir + 'copy_destriper'
        descart = './pcart'
        params = 'default_params.ini'
        for bd in bolodirs:
            sb.check_call([copyscript, bd])
        for bd in bolodirs:
            os.chdir(bd)
            sb.call(['mpirun', '-np', '4', descart, params])
        return 


if __name__ == "__main__":
    mapmaker = RunDescart()
    mapmaker.run()
