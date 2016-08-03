import numpy as np
import pylab as pl
import healpy as hp
import seaborn as sns
from scipy.io import FortranFile
from astropy.io import fits
import glob
from leap.lib.numerical import fourier_analysis as fa
from leap.lib.time_domain_processing import filter_wrapper, filtering
from leap.lib.mapping import naive_binning


#make weights and used add to hit and weights
#mapper lives in lip mapping wrapper, multiprocess_build_map
#moving weights on 20 minute scale

def filter_bolo(f):
    binning_func = naive_binning.add_to_hit_and_weight_and_unnormalized_signal_map
    Nl = hp.nside2npix(512)
    filtmap = np.zeros(Nl, np.double)
    hitmap = np.zeros(Nl, np.uint64)
    weightmap = np.zeros(Nl, np.double)
    fd = fits.open(f)
    data = fd[1].data['CH1'].astype(np.double)
    lons = (fd[1].data['RA'])*np.pi/180.
    lats = (fd[1].data['DEC'])*np.pi/180.
    valid = fd[1].data['FLAG1'].astype('bool')
    fd.close()
    signal = filter_wrapper.do_T_filtering(data, 'data', deconvolve_bolo=False)
    weights = filtering.moving_weights(signal, window_size=1200.)
    binning_func(hitmap, weightmap, filtmap, signal[valid], weights[valid], lons[valid], lats[valid])
    return filtmap, hitmap, weightmap

def run():
    dir1 = "/cache/mabitbol/ebex/goodboloseg250/data/"
    ffiles = glob.glob(dir1+'goodboloseg*.fits')
    Nl = hp.nside2npix(512)
    total_signal = np.zeros(Nl, np.double)
    total_weight = np.zeros(Nl, np.double)
    total_hits = np.zeros(Nl, np.uint64)
    for f in ffiles:
        filteredmap, hitsmap, weightmap = filter_bolo(f)
        total_signal += filteredmap
        total_weight += weightmap
        total_hits += hitsmap
    hp.write_map('filteredsignal.fits', total_signal)
    hp.write_map('filteredhits.fits', total_hits)
    hp.write_map('filteredweights.fits', total_weight)

run()

