import numpy as np
import pylab as pl
import healpy as hp
import timeit
import pink_noise
import fits_writer
from pylab import cos, sin, pi, e, arctan2, sqrt
import glob
import deepdish as dd

from tools_glp import *

from leap.lib.units.angles import *
from leap.lib.geometry import coordinates
from leap.lib.mapping import naive_binning
from leap.lib.numerical import fourier_analysis
from leap.resources.constants import sample_rates
from simple_fitter import psd_fit_errors
rate = sample_rates.bolo

class glp_sim():

    def run(self):
        pointing_dir = 'ebex250sgal/'
        idn = "ebexasnoise"
        freq = 250

        nside = 512
        num_det = 1
        freq = str(freq)
        idn += freq
        savedir = make_dir(idn)
        datadir = make_dir(idn, 'data')+'/'
        make_dir(idn, 'offsets')
        save_name = datadir+idn+"_nside"+str(nside)

        tic = timeit.default_timer()
        delete_txt(savedir)
        pointing_files = glob.glob(pointing_dir+'201*'+freq+'.h5')
        skymap = create_cmb(nside, True, freq)

        bad = np.load('bad_bolos.npy')
        for bolo in bad:
            for pf in pointing_files:
                if bolo in pf:
                    pointing_files.remove(pf)

        print len(pointing_files)
        for i, pfile in enumerate(pointing_files):
            if i % 100 == 0:
                print i
            good, data, pointing, N = self.load_pointing(pfile)
            if good:
                offsets = [0,0]
                mjds = get_timing(pointing['time'])
                net = 1000.e-6
                fknee = 0.2
                alpha = 2.5
                skydata = self.get_tod(nside, skymap, pointing)
                data += skydata
                hdulist = fits_writer.make_fits(N, mjds, data, pointing, offsets, num_det, net, fknee, alpha)
                fits_writer.write_copy_fits(hdulist, save_name+str(i), savedir)
        print "time ", timeit.default_timer()-tic
        return

    def load_pointing(self, pfile):
        # TOD in Kelvin
        good = True
        data = dd.io.load(pfile)
        valid = data[0]['valid']
        times = data[0]['times'][valid]
        tod = data[0]['data'][valid].astype(np.double)
        calib = data[0]['calib'][valid]
        tod *= calib
        #ip = data[0]['ip'][valid]
        #tod -= ip
        tod -= np.mean(tod)

        mask = np.abs(tod) > 20.
        if np.sum(mask) > 1:
            good = False

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
        return good, tod, pointing, N

    def get_tod(self, nside, skymap, pointing):
        hwp = pointing['hwp']
        pix = hp.ang2pix(nside, pointing['theta'], pointing['phi'])
        I = skymap[0][pix]
        Q = skymap[1][pix]
        U = skymap[2][pix]
        tod = I + Q*cos(hwp) + U*sin(hwp)
        return tod

if __name__ == "__main__":
    simulation = glp_sim()
    simulation.run()


