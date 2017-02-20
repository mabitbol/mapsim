import numpy as np
import deepdish as dd
import matplotlib as mpl
#mpl.rcsetup.all_backends
#mpl.use('gtkagg')
from matplotlib import pyplot as pl
import glob as glob

def plot_tod():
    pointing_files = glob.glob('ebex_all_prep/201*250.h5')
    bad = np.load('bad_bolos.npy')
    for bolo in bad:
        for pf in pointing_files:
            if bolo in pf:
                pointing_files.remove(pf)
    bolo_list = []
    segment_list = []
    for f in pointing_files:
        segbolo_str = f.split('/')[-1].split('_')
        segment_list.append(segbolo_str[0])
        bolo_list.append(segbolo_str[1])
    segment_list = list(set(segment_list))
    bolo_list = list(set(bolo_list))

    for bolo in bolo_list:
        bfiles = [f for f in pointing_files if bolo in f]
        for segment in segment_list:
            boloseg_files = [f for f in bfiles if segment in f]
            if len(boloseg_files):
                pl.figure()
                for bseg in boloseg_files:
                    times, valid, tod, lats = load_pointing(bseg)
                    times -= 1356912000
                    mask = np.abs(lats) < 3.*np.pi/180.
                    pl.plot(times[valid & mask]/60., tod[valid & mask], 'g.')
                    pl.plot(times[valid & ~mask]/60., tod[valid & ~mask], 'b.')
                    pl.plot(times[~valid & mask]/60., tod[~valid & mask], 'k,')
                    pl.plot(times[~valid & ~mask]/60., tod[~valid & ~mask], 'r,')
                    #green is valid and on galaxy
                    #blue is valid and off galaxy
                    #black is not valid and on galaxy
                    # red is not valid and off galaxy
                pl.grid()
                pl.xlabel('Time since 12/31/12 00:00:00 [mins]')
                pl.ylabel('TOD [K]')
                pl.title(bolo+' '+segment)
                pl.savefig('todplots250/'+bolo+segment) 
                pl.close()



def load_pointing(pfile):
    data = dd.io.load(pfile)
    valid = data[0]['valid']
    times = data[0]['times']
    lats = data[0]['lats']
    tod = data[0]['data'].astype(np.double)
    return times, valid, tod, lats


plot_tod()


