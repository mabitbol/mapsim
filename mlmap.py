import numpy as np
from pylab import *
import healpy as hp
import pink_noise as pn
import seaborn as sns
from matplotlib import cm

def mlmaptest():
    Ntot = 3
    nl = 10
    nlsq = nl*nl
    ind = np.arange(nl)
    x = (-(2./nl)**2)*ind*ind + (4./nl)*ind
    y = (-(2./nl)**2)*ind*ind + (4./nl)*ind
    sig = np.matrix(x).T*np.matrix(y)
    #matshow(sig)
    #show()

    sigfl = np.array(sig)
    sigfl = sigfl.reshape(nlsq)
    pointing = np.arange(nlsq*Ntot)%(nlsq)
    data = sigfl[pointing]
    #plot(data)
    #show()

    recon = np.zeros([nl,nl])
    for i in range(Ntot):
        recon += data[i*nlsq:(i+1)*nlsq].reshape(nl,nl)
    recon /= Ntot
    #matshow(recon)
    #show()

    std = 1.
    fk = 0.5
    #noise = np.random.normal(0,std,Ntot*nlsq)
    fnoise = pn.pink_noise(Ntot*nlsq, delta_t=1., net=std, fknee=fk, alpha=2, show_plot=False)
    tod = data+fnoise

    recon = np.zeros([nl,nl])
    for i in range(Ntot):
        recon += tod[i*nlsq:(i+1)*nlsq].reshape(nl,nl)
    recon /= np.float(Ntot)
    #matshow(recon)
    #show()

    print "calculating N"
    it = 1000
    #noise = np.random.normal(0,std,Ntot*nlsq)
    #N = np.matrix(noise).T*np.matrix(noise)
    N = np.matrix(fnoise).T*np.matrix(fnoise)
    for i in range(it-1):
        fnoise = pn.pink_noise(Ntot*nlsq, delta_t=1., net=std, fknee=fk, alpha=2, show_plot=False)
        #noise = np.random.normal(0,std,Ntot*nlsq)
        N += np.matrix(fnoise).T*np.matrix(fnoise)
    N /= np.float(it)
    #matshow(N)
    #show()

    print np.linalg.det(N)
    print "calculating Ninv"
    Ninv = np.linalg.inv(N)
    #matshow(Ninv)
    #show()

    idn = np.matmul(N,Ninv)
    matshow(idn)
    show()

    Pmat = np.diag(np.ones(nlsq))
    Pmat = np.concatenate((Pmat,Pmat,Pmat),axis=1)
    Pmat = np.matrix(Pmat)
    
    LH = np.matmul(Pmat, np.matmul(Ninv,Pmat.T))
    LHinv = np.linalg.inv(LH)
    idn = np.matmul(LH, LHinv)
    matshow(idn)
    show()

    todm = np.matrix(tod).T
    RH = np.matmul(Pmat, np.matmul(Ninv, todm))

    mmap = np.matmul(LHinv, RH)
    mmap = mmap.reshape(nl,nl)

    matshow(mmap)
    show()
    
    #matshow(mmap-recon)
    #show()

mlmaptest()
