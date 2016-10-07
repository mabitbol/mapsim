import numpy as np
import timeit
import fits_writer
import pink_noise
import glob
import deepdish as dd
from tools_glp import *
from pylab import cos, sin, pi, e, arctan2, sqrt

class SingleBolo():

    def run(self):
        pointing_dir = 'ebex_galaxy/'
        #pointing_dir = 'goodboloseg_full/'
        idn = "ebexgalbolos"
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

        skymap = create_cmb(nside, True, freq)
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

                data, offsets = self.get_tod(nside, skymap, pointing, dt, net, fknee, alpha)

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

if __name__ == "__main__":
    write = SingleBolo()
    write.run()


