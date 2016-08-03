import numpy as np
import pylab as pl
import healpy as hp
import os
import shutil
from astropy.io import fits
from pylab import pi, array
from leap.lib.units.angles import *


def make_fits(length, mjds, data, pointing, offsets, num_det, w, fk, al):
    print "making fits"
    columns = make_columns(mjds, pointing, data)
    primary_hdu = fits.PrimaryHDU()
    table_binary_hdu = fits.BinTableHDU.from_columns(columns)
    second = second_hdu(length, num_det, w, fk, al)
    third = third_hdu(offsets, num_det)
    return fits.HDUList([primary_hdu, table_binary_hdu, second, third])


def make_columns(time, pointing, cmb_tod):
    print "making columns"
    N = len(time)
    empt = np.zeros(N)
    times = fits.Column(name="MJD", format="D", array=time)
    azimuth = fits.Column(name="AZIMUTH", format="D", array=empt)
    elev = fits.Column(name="ELEVATIO", format="D", array=empt)
    ra_deg = fits.Column(name="RA", format="D", array=to_degrees(pointing['phi']))
    dec_deg = fits.Column(name="DEC", format="D", array=to_degrees(pi/2-pointing['theta']))
    roll_deg = fits.Column(name="PARANGLE", format="D", array=to_degrees(pointing['roll']))
    flag = fits.Column(name="FLAG1", format="I", array=pointing['valid'].astype(int))
    bolo = fits.Column(name="CH1", format="E", array=cmb_tod)
    hwp_angle = fits.Column(name="HWP", format="D", array=to_degrees(pointing['hwp']))
    columns = fits.ColDefs([times, azimuth, elev, ra_deg, dec_deg, roll_deg, flag, hwp_angle, bolo])
    return columns


def second_hdu(len_time, num_det, w, fk, al):
    shell = pl.ones((1,num_det))
    snum = str(num_det)
    start = fits.Column(name="START", format='J', array=pl.array([0]))
    end = fits.Column(name="END", format='J', array=pl.array([len_time-1]))
    sigma = fits.Column(name="SIGMA", format=snum+'D', array=shell*w)
    alpha = fits.Column(name="ALPHA", format=snum+'D', array=shell*(-al))
    fknee = fits.Column(name="FKNEE", format=snum+'D', array=shell*fk)
    col2 = fits.ColDefs([start, end, sigma, alpha, fknee])
    return fits.BinTableHDU.from_columns(col2)


def third_hdu(offsets, num_det):
    az = [offsets[0]]
    el = [offsets[1]]
    bolo_on = pl.ones(num_det)
    onoff = fits.Column(name="ONOFF", format="I", array=bolo_on.astype(int))
    azoff = fits.Column(name="AZOFF", format="E", array=az)
    eloff = fits.Column(name="ELOFF", format="E", array=el)
    gain = fits.Column(name="GAIN", format="E", array=bolo_on)
    delta_p = fits.Column(name="DELTA_PARANGLE", format="E", array=pl.zeros(num_det))
    col3 = fits.ColDefs([onoff, azoff, eloff, gain, delta_p])
    return fits.BinTableHDU.from_columns(col3)


def write_copy_fits(hdulist, name, savedir):
    print "saving files"
    descart_dir = "/home/mabitbol/mapmaking/simulator/"
    descart_fname = name+".fits"
    file_fname = savedir+"files_test.txt"
    fname = os.path.join(descart_dir, descart_fname)
    hdulist.writeto(fname,clobber=True)
    descart = open(os.path.join(descart_dir, file_fname), 'a+')
    descart.write(os.path.join(descart_dir,descart_fname)+'\n')
    descart.close()
    return

