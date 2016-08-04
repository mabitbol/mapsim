import numpy as np
import healpy as hp
from astropy.io import fits
import glob
import os

class BuildMaps():

    def run(self):
        fdir = '/cache/mabitbol/ebex/'
        datadir = fdir + 'bolo250/bolos/'
        bolodirs = glob.glob(datadir+'*')
        nside = 512
        signals = np.zeros(hp.nside2npix(nside))
        hitmap = np.zeros(hp.nside2npix(nside))
        for bd in bolodirs:
            os.chdir(bd)
            destriped = hp.read_map('map.fits')
            hits = self.gethits()
            mask = hits > 0
            signals[mask] += destriped[mask] * hits[mask]
            hitmap[mask] += hits[mask]
        hp.write_map(fdir+'signalmap.fits', signals)
        hp.write_map(fdir+'hitmap.fits', hitmap)
        mask = hitmap > 0 
        totalmap = np.zeros(hp.nside2npix(nside))
        totalmap[mask] = signals[mask] / hitmap[mask]
        hp.write_map(fdir+'totalmap.fits', totalmap)
        return 

    def gethits(self):
        x = np.zeros(hp.nside2npix(512))
        with open('map.hits') as f:
            k = 0
            for line in f:
                x[k] = int(line.split(' ')[-1].split('\n')[0])
                k += 1
        return x
            

if __name__ == "__main__":
    build = BuildMaps()
    build.run()
