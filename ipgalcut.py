import numpy as np
import healpy as hp
import glob
import deepdish as dd
import os
from leap.lib.geometry.coordinates import angular_distance

# deepdish pointing files save as list with dictionaries in them.
# should change this to use h5py asap
# [segname, bolo_name, times, valids, lats, lons, int: frequency]
def read_bolo_pointing(fname):
    data = dd.io.load(fname)
    segment = str(data[0]['segname'])
    bolo_name = str(data[1]['bolo_name'])
    times = data[2]['times']
    tod = data[3]['data']
    # 1 is good for valid
    valid = data[4]['valids']
    lats = data[5]['lats']
    lons = data[6]['lons']
    hwps = data[7]['hwps']
    calib = data[8]['calib']
    ip = data[9]['ip']
    freq = str(data[10])
    return times, tod, valid, lats, lons, hwps, calib, ip, freq, segment, bolo_name, len(times)

def get_galcut(times, tod, lats, lons, valid, hwps, calib, ip):
    sun_lats, sun_lons = get_sun_bl(times)
    dist = angular_distance(sun_lats, sun_lons, lats, lons)
    gc = dist > (np.pi/4)
    return times[gc], tod[gc], lats[gc], lons[gc], hwps[gc], valid[gc], calib[gc], ip[gc]

def save_data(fname, times, tod, lats, lons, hwps, valid, calib, ip):
    dd.io.save(fname+'.h5', [{"times": times, "data": tod, "lats": lats,\
               "lons": lons, "hwps": hwps, "valid":valid, "calib":calib, "ip":ip}])
    return

def galaxy_cut():
    pointing_dir = "ebex250ipraw/"
    galcut_dir = "ebex250suncut"
    make_dir(galcut_dir)
    segments = glob.glob(pointing_dir+'segment*')
    for seg in segments:
        boards = glob.glob(seg+'/board*')
        for board in boards:
            bolos = glob.glob(board+'/201*')
            for bolo in bolos:
                data  = read_bolo_pointing(bolo)
                times, tod, valid, lats, lons, hwps, calib, ip, freq, segname, bolo_name, N = data
                for k in range(N):
                    times0, tod0, lats0, lons0, hwps0, valid0, calib0, ip0 = \
                        get_galcut(times[k], tod[k], lats[k], lons[k], valid[k], hwps[k], calib[k], ip[k])
                    if len(times0[valid0]) > 400:
                        fname = galcut_dir+'/'+segname+'_'+bolo_name+'_'+str(k)+'_'+freq
                        save_data(fname, times0, tod0, lats0, lons0, hwps0, valid0, calib0, ip0)
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

def get_sun_bl(times):
    a = 5.2467176609e-15
    b = -1.41356413129e-05
    c = 9520.70604379
    lons = a*times*times + b*times + c

    a1 = 1.50793236866e-15
    b1 = -4.27070840547e-06
    c1 = 3018.41141079
    lats = a1*times*times + b1*times + c1
    return lats, lons


galaxy_cut()
