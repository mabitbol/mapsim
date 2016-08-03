import numpy as np
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
from leap.resources.constants import sample_rates
rate = sample_rates.bolo

class glp_sim():

    def run(self):
        nside = 512
        num_det = 1
        dt = 1/rate
        freq = 150
        freq = str(freq)

        net = 500.0               #uK sqrt(s)
        fknee = 0.2              #Hz
        alpha = 2.
        beam = 8.0              #arcmins
        beam = from_arcmin(beam)

        self.cmb_on = True
        self.plmaps = True
        idn = "ebexsimhwpflagsnew"
        idn += freq
        self.savedir = make_dir(idn)
        save_name = self.savedir+idn+"_nside"+str(nside)+"_net"+str(int(net))+\
                    "_fk"+str(int(fknee))
        ##########################################################

        tic = timeit.default_timer()
        skymap = create_cmb(nside, True, freq)
        T, Q, U, hits = np.zeros((4,hp.nside2npix(nside)))
        delete_txt(self.savedir)
        
        pointing_dir = 'ebex_galaxy_flags/'
        pointing_files = glob.glob(pointing_dir+'201*'+freq+'.h5')
        for i, pfile in enumerate(pointing_files):
            pointing, N, good = self.load_pointing(pfile)
            if good:
                data, offsets = self.get_tod(nside, skymap, pointing, dt, net, fknee, alpha)
                mjds = get_timing(pointing['time'])
                hdulist = fits_writer.make_fits(N, mjds, data, pointing, offsets, num_det, net, fknee, alpha)
                fits_writer.write_copy_fits(hdulist, save_name+str(i), self.savedir)

                #T1, Q1, U1, hits1 = self.make_maps(nside, pointing, data, num_det)
                #T += T1
                #Q += Q1
                #U += U1
                #hits += hits1
        #T /= hits
        #Q /= hits
        #U /= hits
        #self.plot_maps(T, Q, U, hits, self.savedir)
        print "time ", timeit.default_timer()-tic
        return

    def load_pointing(self, pfile):
        pointing = {}
        data = dd.io.load(pfile)
        valid = data[0]['valid']
        times = data[0]['times'][valid]
        lats = data[0]['lats'][valid]
        lons = data[0]['lons'][valid]
        hwps = data[0]['hwps'][valid]
        N = len(times)
        
        if np.sum(valid) > 0.5*N:
            good = True
        else:
            good = False
        valid = valid[valid] #LOLLLL
        if good:
            roll = np.zeros(N)
            pointing['time'] = times
            pointing['theta'] = pi/2. - lats
            pointing['phi'] = lons
            pointing['valid'] = valid
            pointing['roll'] = roll
            pointing['hwp'] = hwps % (2.*pi)
        return pointing, N, good

    def get_tod(self, nside, skymap, pointing, dt, net, fknee, a):
        print "calculating tod and noise"
        theta = pointing['theta']
        phi = pointing['phi']
        roll = pointing['roll']
        hwp = pointing['hwp']
        pix = hp.ang2pix(nside, theta, phi)
        I = skymap[0][pix]
        Q = skymap[1][pix]
        U = skymap[2][pix]
        tod = I + Q*cos(hwp) + U*sin(hwp)

        N = len(tod)
        fft_len = int(2**(2+np.floor(np.log2(N))))
        fnoise = pink_noise.pink_noise(length=fft_len, delta_t=dt, net=net, fknee=fknee, alpha=a, show_plot=False)
        x = np.random.randint(1,N/2)
        tod1 = tod + fnoise[x:N+x]
        offsets = [0,0]
        return tod1, offsets
    
    def make_maps(self, nside, pointing, data, num_det):
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
        print "saving maps"
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


