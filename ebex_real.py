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
        pointing_dir = 'ebex_galaxy/'
        #pointing_dir = 'goodboloseg_full/'
        idn = "ebex_corr"
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

        print len(pointing_files)
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
            #good, data, pointing, N, net, fknee, alpha = self.load_pointing_noise(pfile)
            if good:
                offsets = [0,0]
                mjds = get_timing(pointing['time'])
                net = 1000.e-6
                fknee = 0.2
                alpha = 2.5
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


    def load_pointing_noise(self, pfile):
        # TOD in Kelvin
        data = dd.io.load(pfile)
        valid = data[0]['valid']
        times = data[0]['times']
        tod = data[0]['data'].astype(np.double)
        lats = data[0]['lats']
        lons = data[0]['lons']
        hwps = data[0]['hwps']

        N = len(times)
        mask = np.abs(tod) > 20.
        if np.sum(mask) > 1:
            return False, 0, 0, 0, 0, 0, 0
        else:
            nets = 0
            fknees = 0
            alphas = 0
            m = 2**14
            L = int(np.floor(N/m))
            if L > 0:
                for i in range(L):
                    x = i*m
                    y = (i+1)*m
                    bad, net, fknee, alpha = self.get_fit(times[x:y], tod[x:y])
                    if bad:
                        valid[x:y] = 0 
            else:
                bad, nets, fknees, alphas = self.get_fit(times, tod)
                if bad:
                    return False, 0, 0, 0, 0, 0, 0
            N = len(times[valid])
            if N < 2**13:
                return False, 0, 0, 0, 0, 0, 0
            
            pointing = {}
            pointing['time'] = times[valid]
            pointing['theta'] = pi/2. - lats[valid]
            pointing['phi'] = lons[valid]
            pointing['hwp'] = hwps[valid] % (2.*pi)
            pointing['roll'] = np.zeros(N)
            pointing['valid'] = valid[valid]
            return True, tod[valid], pointing, N, nets, fknees, alphas

    def get_fit(self, times, tod):
        times_p = times - times[0]
        z = np.polyfit(times_p, tod, 3)
        poly = np.poly1d(z)
        fitpoly = poly(times_p)
        reduced_ts = tod - fitpoly
        ts_win = fourier_analysis.window(reduced_ts)
        psd_f = fourier_analysis.psd_func_real
        psd_clean = psd_f(ts_win, rate, truncate=False)
        params = self.fit(psd_clean, ts_win.std())
        return self.is_bad(*params)
        
    def fit(self, psd, std_in, band=10.):
        freq = psd[0][2:-1]
        psd = 0.5*psd[1][2:-1]
        ind = np.where(freq<band)[0]
        f = freq[ind]
        p = psd[ind]
        return psd_fit_errors(f, p, std_in)

    def is_bad(self, fknee, white, alpha, fkstd, whitestd, alphas):
        bad = False
        white = np.sqrt(white)
        whitestd = 0.5*whitestd/white
        if fknee > 5.0:
            bad = True
        if white > 10000.e-6 or white < 400.e-6:
            bad = True
        if alpha > 5.0:
            bad = True
        return bad, white, fknee, alpha

    def round_valid(self, valid):
        olen = 400
        if len(valid) == olen:
            if np.sum(valid) != olen:
                valid[:] = 0
                return
        if len(valid) > olen:
            if np.sum(valid[:olen]) != olen:
                valid[:olen] = 0
            start = np.where(valid[olen:] == 1)[0]
            if len(start) > 0:
                start = start[0] + olen
                self.round_valid(valid[start:])
            return 
        if len(valid) < olen:
            valid[:] = 0
        return


if __name__ == "__main__":
    simulation = glp_sim()
    simulation.run()


