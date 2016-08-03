import numpy as np
import healpy as hp
import ephem
import os
from pylab import cos, sin, pi, arccos, arcsin, arctan2, array, clip


def my_hor_to_eq(az, el, lat, lsts):
    dec = arcsin(sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az))
    argument = (sin(el) - sin(lat) * sin(dec)) / (cos(lat) * cos(dec))
    argument = clip(argument, -1.0, 1.0)
    H = arccos(argument)
    flag = sin(az) > 0
    H[flag] = 2.0*pi - H[flag]
    ra = lsts - H
    ra %= 2*pi
    return ra,dec

def get_timing(time):
    today_djd = float(ephem.Date('2015/01/01 00:00:00'))
    djds = today_djd+time/(60.0*60.0*24.0)
    return djds+15019.5         #converted from time in seconds to MJD

def create_cmb(nside, polon, freq='0'):
    print "creating map"
    if freq == '0':
        cmbfname = "planck_data/cmb_map512.fits"
        cmb = hp.read_map(cmbfname, field=(0,1,2))
        cmb = hp.ud_grade(cmb, nside)
        #cmb = hp.sphtfunc.smoothing(cmb, fwhm=beam, pol=True)
        cmbnew = np.array(cmb)
        norm = 10**6
        cmbnew[0] *= norm
        cmbnew[1] *= norm
        cmbnew[2] *= norm
        if polon:
            return cmbnew
        else:
            return cmbnew[0]
    else:
        cmbfname = 'planck_data/ebex_'+freq+'.fits'
        cmbT = hp.read_map(cmbfname)
        if polon:
            if freq != '150':
                print "no pol maps at this frequencye"
            Q = hp.read_map('planck_data/ebex_150Q.fits')
            U = hp.read_map('planck_data/ebex_150U.fits')
            return [cmbT, Q, U]
        else:
            return cmbT

def make_dir(dir_name, name=""):
    #harddir = "simdata/"
    harddir = "/cache/mabitbol/ebex/"
    path = harddir+dir_name+"/"+name
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    return path

def delete_txt(savedir):
    #descart_dir = "/home/mabitbol/mapmaking/simulator/"
    descart_dir = "/cache/mabitbol/ebex/"
    fname = descart_dir+savedir+"files_test.txt"
    try:
        os.remove(fname)
    except OSError:
        if os.path.isfile(fname):
            raise
    return

