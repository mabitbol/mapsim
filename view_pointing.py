import numpy as np
import healpy as hp
import matplotlib as mpl
mpl.use('Agg')
import pylab as pl
import glob
import deepdish as dd
from leap.lib.units.angles import *
import os

def plot_hitmap():
    pointing_dir = "ebex_pointing/"
    segments = glob.glob(pointing_dir+'segment*')
    for seg in segments:
        boards = glob.glob(seg+'/board*')
        for board in boards:
            bolos = glob.glob(board+'/2012-*')
            for bolo in bolos:
                data = read_bolo_pointing(bolo)
                plot_bolo_data(data)




# deepdish pointing files save as list with dictionaries in them.
# should change this to use h5py asap
# [segname, bolo_name, times, valids, lats, lons, int: frequency]
def read_bolo_pointing(fname):
    data = dd.io.load(fname)
    segment = str(data[0]['segname'])
    bolo_name = str(data[1]['bolo_name'])
    times = data[2]['times']
    # 1 is good for valid
    valid = data[3]['valids'].astype('bool')
    lats = data[4]['lats']
    lons = data[5]['lons']
    freq = str(data[6])
    return times, valid, lats, lons, freq, segment, bolo_name 


def plot_bolo_data(data):
    times, valid, lats, lons, freq, segment, bolo_name = data
    pl.figure()
    pl.plot(times[valid], to_degrees(lats[valid])+90., 'b,')
    pl.plot(times[~valid], to_degrees(lats[~valid])+90., 'r,')
    pl.plot(times[valid], to_degrees(lons[valid]), 'g,')
    pl.plot(times[~valid], to_degrees(lons[~valid]), 'r,')
    pl.xlim([0,10000])
    pl.ylabel('pointing [deg]')
    pl.xlabel('time [s]')
    pl.show()



def get_galcut(times, lats, lons, valid):
    gc = lats < from_degrees(5.)
    gc &= lats > from_degrees(-5.)
    gc &= valid
    return times[gc], lats[gc], lons[gc]

def save_data(fname, times, lats, lons):
    dd.io.save(fname+'.h5', [{"times": times, "lats": lats,\
               "lons": lons}])
    return
    

def galaxy_cut():
    pointing_dir = "ebex_pointing_hwp/"
    galcut_dir = "galaxy_cut_hwp"
    make_dir(galcut_dir)
    segments = glob.glob(pointing_dir+'segment*')
    for seg in segments:
        boards = glob.glob(seg+'/board*')
        for board in boards:
            bolos = glob.glob(board+'/201*')
            for bolo in bolos:
                data  = read_bolo_pointing(bolo)
                times, valid, lats, lons, freq, segname, bolo_name = data
                times, lats, lons = get_galcut(times, lats, lons, valid)
                if len(times):
                    fname = galcut_dir+'/'+segname+'_'+bolo_name+'_'+freq
                    save_data(fname, times, lats, lons)
                else:
                    print "no valid galaxy pointing"
    return 



def make_dir(dir_name, name=""):
    path = dir_name+"/"+name
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    return path


galaxy_cut()
