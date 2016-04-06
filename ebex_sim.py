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
from leap.resources.constants import sample_rates
rate = sample_rates.bolo

class glp_sim():

    def run(self):
        nside = 512
        num_det = 2
        dt = 1/rate

        self.polarization = False
        hwp_hz = 1.26*pi/e
        hwp_speed = 2.*pi*hwp_hz

        net = 25.               #uK sqrt(s)
        fknee = 0.05             #Hz
        alpha = 2.
        beam = 8.0              #arcmins
        beam = from_arcmin(beam)

        self.cmb_on = True
        self.plmaps = True
        idn = "ebextest"
        self.savedir = make_dir(idn)
        save_name = self.savedir+idn+"_nside"+str(nside)+"_net"+str(net)+\
                    "_fk"+str(int(fknee))
        ##########################################################

        tic = timeit.default_timer()
        skymap = create_cmb(nside, self.polarization)
        T, Q, U, hits = np.zeros((4,hp.nside2npix(nside)))
        delete_txt(self.savedir)
        
        pointing_dir = 'galaxy_cut/'
        pointing_files = glob.glob(pointing_dir+'2012-*')
        for i, pfile in enumerate(pointing_files):
            if self.polarization:
                print "no pol"
                #T += T1
                #Q += Q1
                #U += U1
                #hits += hits1
                return
            else:
                pointing, N = self.load_pointing(pfile)
                data, offsets = self.get_tod(nside, skymap, pointing, dt,\
                            net, fknee, alpha)
                mjds = get_timing(pointing['time'])
                flags = np.ones(N)
                hdulist = fits_writer.make_fits(N, mjds, data, pointing,\
                            offsets, num_det, flags, net, fknee, alpha)
                fits_writer.write_copy_fits(hdulist, save_name+str(i), self.savedir)
                T1, hits1 = self.make_maps(nside, pointing, data, num_det)
                T += T1
                hits += hits1
        T /= hits
        if self.polarization:
            Q /= hits
            U /= hits
        self.plot_maps(T, Q, U, hits, self.savedir)
        print "time ", timeit.default_timer()-tic
        return

    def load_pointing(self, pfile):
        data = dd.io.load(pfile)
        times = data[0]['times']
        lats = data[0]['lats']
        lons = data[0]['lons']
        
        pointing = {}
        pointing['time'] = times
        pointing['theta'] = pi/2. - lats
        pointing['phi'] = lons

        N = len(times)
        roll = np.zeros(N)
        pointing['roll'] = roll
        if self.polarization:
            # fix for polarization
            #hwp = (hwp_speed*time)%(2*pi)   
            pointing['hwp'] = roll
        else:
            pointing['hwp'] = roll
        return pointing, N

    def get_tod(self, nside, skymap, pointing, dt, net, fknee, a):
        print "calculating tod and noise"
        theta = pointing['theta']
        phi = pointing['phi']
        roll = pointing['roll']
        hwp = pointing['hwp']
        pix = hp.ang2pix(nside, theta, phi)
        if self.polarization:
            #I = hp.get_interp_val(skymap[0], theta, phi)
            #Q = hp.get_interp_val(skymap[1], theta, phi)
            #U = hp.get_interp_val(skymap[2], theta, phi)
            I = skymap[0][pix]
            Q = skymap[1][pix]
            U = skymap[2][pix]
            alpha = 0.5*arctan2(U, Q)
            Pin = sqrt(Q**2 + U**2)/I
            tod = I*(1 + Pin*cos(4*hwp - 2*alpha + 2*roll))
            #same as tod = I + Qcos(4hwp+2roll) + Usin(4hwp+2roll)
        else:
            #tod = hp.get_interp_val(skymap, theta, phi)
            tod = skymap[pix]
        N = len(tod)
        fft_len = int(2**(2+np.floor(np.log2(N))))
        fnoise = pink_noise.pink_noise(length=fft_len, delta_t=dt, net=net, fknee=fknee, alpha=a, show_plot=False)
        x = np.random.randint(1,N/2)
        tod1 = tod + fnoise[x:N+x]
        x = np.random.randint(N/2,fft_len-N-1)
        tod2 = tod + fnoise[x:N+x]
        data = [tod1, tod2]
        offsets = [[0,0],[0,0]]
        return data, offsets
    
    def make_maps(self, nside, pointing, data, num_det):
        print "adding maps"
        binning_func = naive_binning.add_to_hit_and_unnormalized_signal_map
        Nl = hp.nside2npix(nside)
        dmoney = (data[0] + data[1])/2.

        I = np.zeros(Nl, np.double) 
        hits = np.zeros(Nl, np.uint64)
        lat = pi/2. - pointing['theta']
        lon = pointing['phi']
        if self.polarization:
            Q = np.zeros(Nl, np.double)
            U = np.zeros(Nl, np.double)
            hits2 = np.zeros(Nl, np.uint64)
            qdata = 0.5*dmoney*cos(4*pointing['hwp']+2*pointing['roll'])
            udata = 0.5*dmoney*sin(4*pointing['hwp']+2*pointing['roll'])
            binning_func(hits, I, dmoney, lon, lat)
            binning_func(hits2, Q, qdata, lon, lat)
            binning_func(hits2, U, udata, lon, lat)
            return I, Q, U, hits
        else:
            binning_func(hits, I, dmoney, lon, lat)
            return I, hits

    def plot_maps(self, I, Q, U, hits, idn):
        print "plotting maps"
        #cmap = cm.YlGnBu_r
        cmap = cm.jet       #lol
        cmap.set_under("0.5")

        mask = hits==0
        hits[mask] = hp.UNSEEN
        I[mask] = hp.UNSEEN
        Q[mask] = hp.UNSEEN
        U[mask] = hp.UNSEEN
        np.savez(idn+"in_maps", hits=hits, I=I, Q=Q, U=U)

        pl.figure()
        hp.mollview(hits, cmap=cmap, unit="counts")
        pl.title("Hit map")
        pl.savefig(idn+"hitmap")
        pl.close()

        pl.figure()
        hp.mollview(I, cmap=cmap, unit="$\mu K$")
        pl.title("T map")
        pl.savefig(idn+"Tmap")
        pl.close()

        if self.polarization:
            pl.figure()
            hp.mollview(Q, cmap=cmap, unit="$\mu K$")
            pl.title("Q map")
            pl.savefig(idn+"Qmap")
            pl.close()

            pl.figure()
            hp.mollview(U, cmap=cmap, unit="$\mu K$")
            pl.title("U map")
            pl.savefig(idn+"Umap")
            pl.close()
        return
        


if __name__ == "__main__":
    simulation = glp_sim()
    simulation.run()


