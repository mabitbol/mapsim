import numpy as np
import matplotlib as mpl
mpl.use('gtkagg')
from matplotlib import cm
import pylab as pl
import healpy as hp
import ephem
import timeit
import pink_noise
import fits_writer
import seaborn as sns
from leap.lib.tools import leap_multiprocessing
from leap.lib.timing import progress_indicator

from pyprimes import factors
from pylab import cos, sin, pi, e, arccos, arcsin, arctan2, sqrt, array, floor, ceil

from leap.lib.units.angles import *
from leap.lib.geometry import coordinates

class glp_sim():

    def run(self):
        nside = 256
        num_det = 2
        fsamp = 100. #Hz
        dt = 1/fsamp 
        days = 6
        lat = from_degrees(-75.0)
        omega = from_degrees(1.0)          #rad per second
        self.fullsky = False

        self.polarization = True
        hwp_hz = 1.26*pi/e
        hwp_speed = 2.*pi*hwp_hz

        pol_frac = 0.1
        signal = 50.0              
        num_gs = 10

        net = 50               #uK sqrt(s)
        fknee = 0.5             #Hz
        alpha = 2.0
        beam = 8.0              #arcmins
        beam = from_degrees(beam/60.)

        self.cmb_on = True
        self.plmaps = True
        idn = "multired_"
        save_name = idn+str(self.polarization)+"_nside"+\
                    str(nside)+"_cmb_"+str(self.cmb_on)+"_days"+\
                    str(days)+"_fullsky_"+str(self.fullsky)
        ##########################################################

        tic = timeit.default_timer()
        N_sec = 24*60*60
        N = int(N_sec / dt) 
        nl = hp.nside2npix(nside)
        T, Q, U, hits = np.zeros((4,hp.nside2npix(nside)))

        num_proc = days
        run_par = True

        pool = leap_multiprocessing.Pool(num_proc, run_parallel=run_par,\
                    maxtasksperchild=1, use_queue=True, shuffle=True)

        for i in range(days):
            pool.add_job(self.run_on_day, args=[N, omega, hwp_speed, dt,\
                         lat, beam, nside, net, fknee,\
                         alpha, num_det, save_name, i])
        if run_par:
            pool.start()
        for job in progress_indicator.wrapper(pool.jobs, verb='doing day'):
            T1, Q1, U1, hits1 = job.get(np.inf)
            print "doing jobs1"
            T += T1
            Q += Q1
            U += U1
            hits += hits1
        pool.end()
        self.plot_maps(self.T, self.Q, self.U, self.hits, idn)
        print "time ", timeit.default_timer()-tic
        return

    def run_on_day(self, N, omega, hwp_speed, dt, lat, beam, nside,\
                   net, fknee, alpha, num_det, save_name, i):
        start = i*N    
        pointing = self.get_pointing_day(start, omega, hwp_speed,\
                        dt, lat, beam)
        skymap = self.create_cmb(nside)
        data, offsets = self.get_tod(nside, pointing, skymap, dt,\
                            net, fknee, alpha)
        mjds = self.get_timing(pointing['time'])
        flags = np.ones(N)
        hdulist = fits_writer.make_fits(N, mjds, data, pointing,\
                    offsets, num_det, flags)
        fits_writer.write_copy_fits(hdulist, save_name+str(i))
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
        ra, dec = self.my_hor_to_eq(az,el,lat,lsts)
        ra, dec = coordinates.eq_to_gal(ra,dec)
        hwp = (hwp_speed*time)%(2*pi)   
        roll = np.zeros(N)
        pointing['time'] = time
        pointing['theta'] = pi/2. - dec
        pointing['phi'] = ra
        pointing['roll'] = roll
        pointing['hwp'] = hwp
        return pointing

    def my_hor_to_eq(self, az, el, lat, lsts):
        dec = arcsin(sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az))
        argument = (sin(el) - sin(lat) * sin(dec)) / (cos(lat) * cos(dec))
        argument = pl.clip(argument, -1.0, 1.0)
        H = arccos(argument)
        flag = sin(az) > 0
        H[flag] = 2.0*pi - H[flag]
        ra = lsts - H
        ra %= 2*pi
        return ra,dec

    def get_timing(self, time):
        today_djd = float(ephem.Date('2015/01/01 00:00:00'))
        djds = today_djd+time/(60.0*60.0*24.0)
        return djds+15019.5         #converted from time in seconds to MJD

    def create_cmb(self, nside):
        print "creating map"
        cmbfname = "planck_data/cmb_map512.fits"
        cmb = hp.read_map(cmbfname, field=(0,1,2))
        cmb = hp.ud_grade(cmb, nside)
        #cmb = hp.sphtfunc.smoothing(cmb, fwhm=beam, pol=True)
        cmbnew = np.array(cmb)
        norm = 10**6
        cmbnew[0] *= norm
        cmbnew[1] *= norm
        cmbnew[2] *= norm
        if self.polarization:
            return cmbnew
        else:
            return cmbnew[0]

    def get_tod(self, nside, pointing, skymap, dt, net, fknee, a):
        print "calculating tod and noise"
        theta = pointing['theta']
        phi = pointing['phi']
        roll = pointing['roll']
        hwp = pointing['hwp']
        pix = hp.ang2pix(nside, theta, phi)
        if self.polarization:
            I = skymap[0][pix]
            Q = skymap[1][pix]
            U = skymap[2][pix]
            alpha = 0.5*arctan2(U, Q)
            Pin = sqrt(Q**2 + U**2)/I
            tod = I*(1 + Pin*cos(4*hwp - 2*alpha + 2*roll))
            #same as tod = I + Qcos(4hwp+2roll) + Usin(4hwp+2roll)
        else:
            tod = skymap[pix]
        N = len(tod)
        fft_len = int(2**(2+np.floor(np.log2(N))))
        fnoise = pink_noise.pink_noise(length=fft_len, delta_t=dt, net=net, fknee=fknee, alpha=a, show_plot=False)
        #fnoise2 = pink_noise.pink_noise(length=fft_len, delta_t=dt, net=net, fknee=fknee, alpha=a, show_plot=False)
        x = np.random.randint(1,2*N)
        tod1 = tod + fnoise[x:N+x]
        x = np.random.randint(2*N,3*N)
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
            mask = hits==0
            I = I/hits
            Q = Q/hits
            U = U/hits
            I[mask] = 0.
            Q[mask] = 0.
            U[mask] = 0.
            return I, Q, U, hits
        else:
            I = np.bincount(pix, weights=dmoney, minlength=Nl)
            hits = np.bincount(pix, weights=np.ones(length)*num_det, minlength=Nl)
            mask = hits==0
            I = I/hits
            I[mask] = 0.
            return I, hits
        return

    def plot_maps(self, I, Q, U, hits, idn):
        print "plotting maps"
        cmap = cm.YlGnBu_r
        #cmap = cm.jet       #lol
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

        #pl.show() 
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
