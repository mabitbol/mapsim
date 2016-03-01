import numpy as np
import matplotlib as mpl
#mpl.use('gtkagg')
from matplotlib import cm
import pylab as pl
import healpy as hp
import timeit
import pink_noise
import fits_writer
import seaborn as sns
from pylab import cos, sin, pi, e, arctan2, sqrt

from tools_glp import *

from leap.lib.units.angles import *
from leap.lib.geometry import coordinates

class glp_sim():

    def run(self):
        nside = 512
        num_det = 2
        fsamp = 11. #Hz
        dt = 1/fsamp 
        days = 10
        lat = from_degrees(-75)
        omega = from_degrees(0.5)          #rad per second

        self.polarization = True
        hwp_hz = 1.26*pi/e
        hwp_speed = 2.*pi*hwp_hz

        net = 25.               #uK sqrt(s)
        fknee = 0.5             #Hz
        alpha = 1.5
        beam = 8.0              #arcmins
        beam = from_degrees(beam/60.)

        self.cmb_on = True
        self.plmaps = True
        idn = "pred"
        self.savedir = make_dir(idn)
        save_name = self.savedir+idn+"_nside"+str(nside)+"_days"+str(days)+\
                    "_net"+str(net)+"_fk"+str(int(fknee))
        ##########################################################

        tic = timeit.default_timer()
        N_sec = 24*60*60
        N = int(N_sec / dt) 
        skymap = create_cmb(nside, self.polarization)
        T, Q, U, hits = np.zeros((4,hp.nside2npix(nside)))
        delete_txt(self.savedir)
        for i in range(days):
            if self.polarization:
                T1, Q1, U1, hits1 = self.run_on_day(N, omega, hwp_speed,\
                                dt, lat, beam, nside, skymap, net,\
                                fknee, alpha, num_det, save_name, i)
                T += T1
                Q += Q1
                U += U1
                hits += hits1
            else:
                T1, hits1 = self.run_on_day(N, omega, hwp_speed, dt,\
                            lat, beam, nside, skymap, net, fknee,\
                            alpha, num_det, save_name, i)
                T += T1
                hits += hits1
        T /= hits
        if self.polarization:
            Q /= hits
            U /= hits
        self.plot_maps(T, Q, U, hits, self.savedir)
        print "time ", timeit.default_timer()-tic
        return

    def run_on_day(self, N, omega, hwp_speed, dt, lat, beam, nside,\
                   skymap, net, fknee, alpha, num_det, save_name, i):
        start = i*N    
        pointing = self.get_pointing_day(start, omega, hwp_speed,\
                        dt, lat, beam)
        data, offsets = self.get_tod(nside, skymap, pointing, dt,\
                            net, fknee, alpha)
        mjds = get_timing(pointing['time'])
        flags = np.ones(N)
        hdulist = fits_writer.make_fits(N, mjds, data, pointing,\
                    offsets, num_det, flags, net, fknee, alpha)
        fits_writer.write_copy_fits(hdulist, save_name+str(i), self.savedir)
        return self.make_maps(nside, pointing, data, num_det)

    def get_pointing_day(self, start, omega, hwp_speed, dt, lat, beam):
        #start on equator. elevation goes continuously 0 to pi/2 in a day
        #az does circles about elevation.  earth rotates around
        print "creating pointing"
        lst_rad_per_sec = 7.292117834034606e-05  #equals omega earth
        N_sec = 24*60*60
        N = int(N_sec / dt)

        pointing = {}
        time = start+np.arange(N)*dt      
        lsts = time*lst_rad_per_sec
        az = omega*time
        az %= (2*pi)
        el = np.ones(N)*from_degrees(65)
        ra, dec = my_hor_to_eq(az,el,lat,lsts)
        ra, dec = coordinates.eq_to_gal(ra,dec)
        hwp = (hwp_speed*time)%(2*pi)   
        roll = np.zeros(N)
        pointing['time'] = time
        pointing['theta'] = pi/2. - dec
        pointing['phi'] = ra
        pointing['roll'] = roll
        if self.polarization:
            pointing['hwp'] = hwp
        else:
            pointing['hwp'] = roll
        return pointing

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
        #fnoise2 = pink_noise.pink_noise(length=fft_len, delta_t=dt, net=net, fknee=fknee, alpha=a, show_plot=False)
        x = np.random.randint(1,N/2)
        tod1 = tod + fnoise[x:N+x]
        x = np.random.randint(N/2,fft_len-N-1)
        tod2 = tod + fnoise[x:N+x]
        data = [tod1, tod2]
        offsets = [[0,0],[0,0]]
        return data, offsets
    
    def make_maps(self, nside, pointing, data, num_det):
        print "adding maps"
        Nl = hp.nside2npix(nside)
        pix = hp.ang2pix(nside, pointing['theta'], pointing['phi'])
        length = len(pix)
        dmoney = data[0] + data[1]
        if self.polarization:
            qdata = 0.5*dmoney*cos(4*pointing['hwp']+2*pointing['roll'])
            udata = 0.5*dmoney*sin(4*pointing['hwp']+2*pointing['roll'])
            I = np.bincount(pix, weights=dmoney, minlength=Nl)
            Q = np.bincount(pix, weights=qdata, minlength=Nl)
            U = np.bincount(pix, weights=udata, minlength=Nl)
            hits = np.bincount(pix, weights=np.ones(length)*num_det, minlength=Nl)
            return I, Q, U, hits
        else:
            I = np.bincount(pix, weights=dmoney, minlength=Nl)
            hits = np.bincount(pix, weights=np.ones(length)*num_det, minlength=Nl)
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



#el = np.ones(N)*pi/2.
#el_speed = (pi/2.)/(24.*60.*60.)
#el = el_speed*time
#el %= pi 
#x = el>(pi/2.)
#el[x] = pi - el[x]
#el %= (pi - 2*res)    #minus res here so that we dont point straight up.
#x = el>(pi/2. - res)
#el[x] = (pi - 2*res) - el[x]
