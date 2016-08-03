import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import cm
import pylab as pl
import healpy as hp
import timeit
import pink_noise
import fits_writer
import seaborn as sns
from pylab import cos, sin, pi, e, arctan2, sqrt

import glob
import deepdish as dd

from tools_glp import *

from leap.lib.units.angles import *
from leap.lib.geometry import coordinates
from leap.lib.mapping import naive_binning
from leap.lib.numerical import fourier_analysis
from leap.resources.constants import sample_rates
rate = sample_rates.bolo

class glp_sim():

    def run(self):
        pointing_dir = 'ebex_galaxy/'
        idn = "ebex"
        freq = 150

        nside = 512
        num_det = 1
        freq = str(freq)
        plmaps = False
        idn += freq
        savedir = make_dir(idn)
        datadir = make_dir(idn, 'data')+'/'
        make_dir(idn, 'offsets')
        save_name = datadir+idn+"_nside"+str(nside)

        tic = timeit.default_timer()
        delete_txt(savedir)
        T, Q, U, hits = np.zeros((4,hp.nside2npix(nside)))
        pointing_files = glob.glob(pointing_dir+'201*'+freq+'.h5')
        net, fknee, alpha = self.get_noise(freq)
        for i, pfile in enumerate(pointing_files):
            data, pointing, N, good = self.load_pointing(pfile)
            if good:
                offsets = [0, 0]
                mjds = get_timing(pointing['time'])
                hdulist = fits_writer.make_fits(N, mjds, data, pointing, offsets, num_det, net, fknee, alpha)
                fits_writer.write_copy_fits(hdulist, save_name+str(i), savedir)

                if plmaps:
                    T1, Q1, U1, hits1 = self.make_maps(nside, pointing, data)
                    T += T1
                    Q += Q1
                    U += U1
                    hits += hits1
        if plmaps:
            T /= hits
            Q /= hits
            U /= hits
            self.plot_maps(T, Q, U, hits, savedir)
        print "time ", timeit.default_timer()-tic
        return

    def get_noise(self, freq):
        if freq == '150':
            net = 500.e-6
        if freq == '250':
            net = 1000.e-6
        fknee = 0.2
        alpha = 2.0
        return net, fknee, alpha

    def load_pointing(self, pfile):
        data = dd.io.load(pfile)
        good = True
        valid = data[0]['valid']
        times = data[0]['times'][valid]
        # TOD in Kelvin
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
        return tod, pointing, N, good

    def load_pointing_noise(self, pfile):
        data = dd.io.load(pfile)
        good = True
        valid = data[0]['valid']
        times = data[0]['times']
        N = len(times)
        # TOD in Kelvin
        tod = data[0]['data'].astype(np.double)
        mask = np.abs(tod) > 20.
        if np.sum(mask) > 1:
            good = False
        else:
            m = 2**16 
            L = int(np.floor(N/m))
            r = N-L*m
            for i in range(L):
                x = i*m
                y = (i+1)*m
                bad = get_fit(times[x:y], tod[x:y])
                if bad:
                    good = False
                    break
            lats = data[0]['lats']
            lons = data[0]['lons']
            hwps = data[0]['hwps']
            
            pointing = {}
            pointing['time'] = times
            pointing['theta'] = pi/2. - lats
            pointing['phi'] = lons
            pointing['hwp'] = hwps % (2.*pi)
            pointing['roll'] = np.zeros(N)
            pointing['valid'] = valid
        return tod, pointing, N, good

    def get_fit(times, tod):
        times_p = times - times[0]
        z = np.polyfit(times_p, tod, 3)
        poly = np.poly1d(z)
        fitpoly = poly(times_p)
        reduced_ts = tod - fitpoly
        ts_win = fourier_analysis.window(reduced_ts)
        psd_f = fourier_analysis.psd_func_real
        psd_clean = psd_f(ts_win, rate, truncate=False)
        params = self.fit(psd_clean, ts_win.std())
        return is_bad(params)
        
    def fit(self, psd, std_in, band=10.):
        freq = psd[0][2:-1]
        psd = 0.5*psd[1][2:-1]
        ind = np.where(freq<band)[0]
        f = freq[ind]
        p = psd[ind]
        params = psd_fit_errors(f, p, std_in)
        fknee, white, alpha, fkstd, whitestd, alphastd = params
        white = np.sqrt(white)
        whitestd = 0.5*whitestd/white
        return fknee, white, alpha, fkstd, whitestd, alphastd

    def is_bad(params):
        fknee, white, alpha, fkstd, whitestd, alphastd = params
        bad = False
        if fknee > 1.0 or fknee == 0.:
            bad = True
        if white > 2000.e-6:
            bad = True
        if alpha > 3.0 or alpha == 0.:
            bad = True
        return bad

    def make_maps(self, nside, pointing, data):
        print "adding maps"
        binning_func = naive_binning.add_to_hit_and_unnormalized_signal_map
        Nl = hp.nside2npix(nside)
        valid = pointing['valid']
        hwp = pointing['hwp'][valid]

        I = np.zeros(Nl, np.double) 
        hits = np.zeros(Nl, np.uint64)
        lat = pi/2. - pointing['theta'][valid]
        lon = pointing['phi'][valid]

        Q = np.zeros(Nl, np.double)
        U = np.zeros(Nl, np.double)
        hits2 = np.zeros(Nl, np.uint64)
        qdata = 2.0 * data[valid] * cos(hwp)
        udata = 2.0 * data[valid] * sin(hwp)
        binning_func(hits, I, data[valid], lon, lat)
        binning_func(hits2, Q, qdata, lon, lat)
        binning_func(hits2, U, udata, lon, lat)
        return I, Q, U, hits

    def plot_maps(self, I, Q, U, hits, idn):
        mask = hits==0
        hits[mask] = hp.UNSEEN
        I[mask] = hp.UNSEEN
        Q[mask] = hp.UNSEEN
        U[mask] = hp.UNSEEN
        np.savez(idn+"in_maps", hits=hits, I=I, Q=Q, U=U)
        return
        


if __name__ == "__main__":
    simulation = glp_sim()
    simulation.run()


