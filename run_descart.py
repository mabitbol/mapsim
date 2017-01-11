import numpy as np
import glob
import subprocess as sb
import os

class RunDescart():

    def run(self):
        fdir = '/cache/mabitbol/ebex/'
        datadir = fdir + 'ebexbolos250/bolo_segments/'
        bolodirs = glob.glob(datadir+'*')
        copyscript = fdir + 'copy_destriper'
        descart = './pcart'
        params = 'default_params.ini'
        for bd in bolodirs:
            sb.check_call([copyscript, bd])
        for bd in bolodirs:
            os.chdir(bd)
            nfiles = len(glob.glob('data/*.fits'))
            nproc = self.calc_nproc(nfiles)
            if nproc:
                sb.call(['mpirun', '-np', nproc, descart, params])
            else:
                sb.call([descart, params])
        return 

    def calc_nproc(self, nfiles):
        if nfiles < 2:
            return '0'
        if nfiles < 25:
            return str(int(nfiles/2) * 2)
        else:
            return '24'


if __name__ == "__main__":
    mapmaker = RunDescart()
    mapmaker.run()
