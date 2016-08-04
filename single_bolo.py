import numpy as np
import timeit
import fits_writer
import glob
import deepdish as dd
from tools_glp import *

class SingleBolo():

    def run(self):
        pointing_dir = 'ebex_galaxy/'
        #pointing_dir = 'goodboloseg_galaxy/'
        idn = "bolo"
        freq = 250
        freq = str(freq)
        idn += freq

        nside = 512
        num_det = 1
        net = 1000.e-6
        fknee = 0.2
        alpha = 2.0

        savedir = make_dir(idn)
        datadir = make_dir_path(savedir, 'bolos')
        tic = timeit.default_timer()
        pointing_files = glob.glob(pointing_dir+'201*'+freq+'.h5')
        bolo_list = []
        for f in pointing_files:
            bolo_list.append(f.split('_')[2])
        bolo_list = list(set(bolo_list))
        for bolo in bolo_list:
            bfiles = [f for f in pointing_files if bolo in f]
            if len(bfiles) == 0:
                break
            bolodir = make_dir_path(datadir, bolo)
            bolodatadir = make_dir_path(bolodir, 'data')
            offsetdir = make_dir_path(bolodir, 'offsets')
            delete_txt(bolodir)
            for bf in bfiles:
                bolosegname = bolodatadir+bf.split('/')[-1].split('.')[0]
                data, pointing, N = self.load_pointing(bf)
                offsets = [0,0]
                mjds = get_timing(pointing['time'])
                hdulist = fits_writer.make_fits(N, mjds, data, pointing, offsets, num_det, net, fknee, alpha)
                fits_writer.write_copy_fits(hdulist, bolosegname, bolodir)
        print "time ", timeit.default_timer()-tic
        return

    def load_pointing(self, pfile):
        # TOD in Kelvin
        data = dd.io.load(pfile)
        valid = data[0]['valid']
        times = data[0]['times'][valid]
        tod = data[0]['data'][valid].astype(np.double)
        lats = data[0]['lats'][valid]
        lons = data[0]['lons'][valid]
        hwps = data[0]['hwps'][valid]
        valid = valid[valid]
        N = len(times)
        pointing = {}
        pointing['time'] = times
        pointing['theta'] = pi/2. - lats
        pointing['phi'] = lons
        pointing['hwp'] = hwps % (2.*pi)
        pointing['roll'] = np.zeros(N)
        pointing['valid'] = valid
        return tod, pointing, N


if __name__ == "__main__":
    write = SingleBolo()
    write.run()


