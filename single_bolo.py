import numpy as np
import timeit
import fits_writer
import glob
import deepdish as dd
from tools_glp import *

class SingleBolo():

    def run(self):
        pointing_dir = 'ebex_all_prep/'
        idn = "ebexbolos"
        freq = 250
        freq = str(freq)
        idn += freq

        nside = 512
        num_det = 1
        net = 1000.e-6
        fknee = 0.2
        alpha = 2.2

        savedir = make_dir(idn)
        datadir = make_dir_path(savedir, 'bolo_segments')
        tic = timeit.default_timer()
        pointing_files = glob.glob(pointing_dir+'201*'+freq+'.h5')

        bad = np.load('bad_bolos.npy')
        for bolo in bad:
            for pf in pointing_files:
                if bolo in pf:
                    pointing_files.remove(pf)

        bolo_list = []
        segment_list = []
        for f in pointing_files:
            segbolo_str = f.split('/')[-1].split('_')
            segment_list.append(segbolo_str[0])
            bolo_list.append(segbolo_str[1])
        segment_list = list(set(segment_list))
        bolo_list = list(set(bolo_list))
        for k, bolo in enumerate(bolo_list):
            if k % 10 == 0:
                print k
            bfiles = [f for f in pointing_files if bolo in f]
            bolo_dir = make_dir_path(datadir, bolo)
            bolodatadir = make_dir_path(bolo_dir, 'data')
            offsetdir = make_dir_path(bolo_dir, 'offsets')
            delete_txt(bolo_dir)
            for segment in segment_list:
                boloseg_files = [f for f in bfiles if segment in f]
                if len(boloseg_files) > 0:
                    for boloseg in boloseg_files:
                        boloseg_filename = bolodatadir+boloseg.split('/')[-1].split('.')[0]
                        data, pointing, N = self.load_pointing(boloseg)
                        offsets = [0,0]
                        mjds = get_timing(pointing['time'])
                        hdulist = fits_writer.make_fits(N, mjds, data, pointing, offsets, num_det, net, fknee, alpha)
                        fits_writer.write_copy_fits(hdulist, boloseg_filename, bolo_dir)
        print "time ", timeit.default_timer()-tic
        return

    def load_pointing(self, pfile):
        # TOD in Kelvin
        data = dd.io.load(pfile)
        valid = data[0]['valid']
        times = data[0]['times'][valid]
        tod = data[0]['data'][valid].astype(np.double)
        tod -= np.mean(tod)

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


