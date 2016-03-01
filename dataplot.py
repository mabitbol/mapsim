import numpy as np
import healpy as hp
import pylab as pl
import seaborn as sns
from matplotlib import cm

def foo1():
    nside = 128
    cmb = hp.read_map("planck_data/cmb_map512.fits", field=(0,1,2))
    cmb = hp.ud_grade(cmb, 256)
    norm = 1e6
    T = cmb[0]*norm
    Q = cmb[1]*norm
    U = cmb[2]*norm

    fdir = "simdata/nonoise/"
    #fw = np.load("analysis/white_in_maps.npz")
    fr = np.load(fdir+"in_maps.npz")
    hits = fr['hits']
    mask = hits<1

    Thit = T
    Thit[mask] = hp.UNSEEN
    #Thit = hp.ud_grade(Thit, nside)

    cmap = cm.jet
    cmap.set_under('0.5')

    pl.figure()
    hp.mollview(Thit, unit='$\mu K$', cmap=cmap)
    pl.title('CMB only')
    pl.savefig('cmbonly')
    pl.close()

    #Tw = fw['I']
    Tr = fr['I']
    #Tw = hp.ud_grade(fw['I'], nside)  
    #Tr = hp.ud_grade(fr['I'], nside)  

    #pl.figure()
    #hp.mollview(Tw, unit='$\mu K$')
    #pl.title('T - white noise')
    #pl.savefig('white_t')
    #pl.close()
    
    pl.figure()
    hp.mollview(Tr, unit='$\mu K$')
    pl.title('T: red noise')
    pl.savefig('red_t')
    pl.close()

    #cmbw_r = Tw - Thit
    #cmbw_r[mask] = hp.UNSEEN
    cmbr_r = Tr - Thit
    cmbr_r[mask] = hp.UNSEEN
    
    #pl.figure()
    #hp.mollview(cmbw_r, unit='$\mu K$')
    #pl.title('white - cmb residuals')
    #pl.savefig('cmbw_r')
    #pl.close()

    pl.figure()
    hp.mollview(cmbr_r, unit='$\mu K$')
    pl.title('red - cmb residuals')
    pl.savefig('cmbr_r')
    pl.close()

    #pl.figure()
    #cmbw_frac = cmbw_r / Thit
    #cmbw_frac[mask] = hp.UNSEEN
    #hp.mollview(cmbw_frac, unit='fraction', min=-1, max=1)
    #pl.title('(white - cmb) / cmb fraction')
    #pl.savefig('cmbw_frac')
    #pl.close()

    pl.figure()
    cmbr_frac = cmbr_r / Thit
    cmbr_frac[mask] = hp.UNSEEN
    hp.mollview(cmbr_frac, unit='fraction', min=-10, max=10)
    pl.title('(red - cmb) / cmb fraction')
    pl.savefig('cmbr_frac')
    pl.close()


    Thit = hp.ud_grade(Thit, nside)
    hits = hp.ud_grade(hits, nside)
    mask = hits<1
    desT = hp.read_map(fdir+"map.fits", field=0)
    nvT = hp.read_map(fdir+"naive_map.fits", field=0)

    pl.figure()
    hp.mollview(desT, unit='$\mu K$')
    pl.title('destriped')
    pl.savefig('des_T')
    pl.close()

    descmb_r = desT - Thit
    descmb_r[mask] = hp.UNSEEN
    nvcmb_r = nvT - Thit
    nvcmb_r[mask] = hp.UNSEEN
     
    descmb_frac = descmb_r / Thit
    descmb_frac[mask] = hp.UNSEEN
    
    nvcmb_frac = nvcmb_r / Thit
    nvcmb_r[mask] = hp.UNSEEN

    pl.figure()
    hp.mollview(descmb_r, unit="$\mu K$")
    pl.title('destriped - cmb residuals')
    pl.savefig('descmb_r')
    pl.close()

    pl.figure()
    hp.mollview(cmbr_frac, unit='fraction', min=-10, max=10)
    pl.title('(destriped - cmb) / cmb fraction')
    pl.savefig('descmb_frac')
    pl.close()

    return


foo1()
