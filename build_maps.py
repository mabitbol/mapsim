import numpy as np
import healpy as hp
from astropy.io import fits
import glob
import os

class BuildMaps():

    def run(self):
        fdir = '/cache/mabitbol/ebex/'
        datadir = fdir + 'ebexbolos250/bolo_segments/'
        bolodirs = glob.glob(datadir+'*')
        nside = 512
        bad = []
        testmap = np.zeros((3, hp.nside2npix(nside)))
        testhitmap = np.zeros((3, hp.nside2npix(nside)))
        for bd in bolodirs:
            if (bd+'/map.fits') in glob.glob(bd+'/*'):
                destriped = hp.read_map(bd+'/map.fits', field=(0,1,2), verbose=False)
                cov = hp.read_map(bd+'/cov.fits', field=(0,1,2),verbose=False)
                mask = destriped[0] != hp.UNSEEN
                k=0
                testmap[k][mask] += destriped[k][mask] / cov[k][mask]
                testhitmap[k][mask] += 1./cov[k][mask]
                k=1
                testmap[k][mask] += destriped[k][mask] / cov[k][mask]
                testhitmap[k][mask] += 1./cov[k][mask]
                k=2
                testmap[k][mask] += destriped[k][mask] / cov[1][mask]
                testhitmap[k][mask] += 1./cov[1][mask]
            else:
                bad.append(bd)

        print len(bad)
        print len(bolodirs)
        mask = testhitmap[0] > 0
        coadded = np.zeros((3, hp.nside2npix(nside)))
        for k in range(3):
            coadded[k][mask] = testmap[k][mask] / testhitmap[k][mask]
            coadded[k][~mask] = hp.UNSEEN

        hp.write_map(datadir+'../totalmap.fits', coadded)
        hp.write_map(datadir+'../covmap.fits', testhitmap)
        hp.write_map(datadir+'../rawmap.fits', testmap)
        return 


    def gethits(self, bd):
        x = np.zeros(hp.nside2npix(512))
        with open(bd+'/map.hits') as f:
            k = 0
            for line in f:
                x[k] = int(line.split(' ')[-1].split('\n')[0])
                k += 1
        return x
            

if __name__ == "__main__":
    build = BuildMaps()
    build.run()
