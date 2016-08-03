#! /usr/bin/python

import numpy as np
import lmfit
from numpy import std, log, exp, arange, pi
from numpy.fft import fft
euler_mascheroni = 0.57721566490153286060
A1 = np.sqrt(pi*pi/6.)

"""This is the core of the fitting function. It contains
the procedure for fitting a 1/f model to a PSD in log linear
space. Corrected for the chi^2(2 dof) distribution."""

def _residuals(p, f, logP):
    fk = p["f_knee"].value
    log_w = p["log_white"].value
    alpha = p["alpha"].value
    model=log_w+log(1.+(fk/f)**alpha)
    return (logP-model+euler_mascheroni)/A1

def _residuals_white(p, f, logP):
    log_w = p["log_white"].value
    model = log_w*np.ones(len(f))
    return (logP-model+euler_mascheroni)/A1

def psd_fit(f, psd, white=1., ftol=1e-8):
    log_Pj = log(psd)
    params = lmfit.Parameters()
    params.add("f_knee", value=0.1, min=0.0, max=25.)
    params.add("log_white", value=log(white), vary=True)
    params.add("alpha", value=2., min=0.0, max=20.)
    out = lmfit.minimize(_residuals, params, args=(f, log_Pj), ftol=ftol, scale_covar=True)
    fk = params["f_knee"].value
    log_w = params["log_white"].value
    alpha = params["alpha"].value
    if (fk>10.) or (alpha>10.):
        params = lmfit.Parameters()
        params.add("log_white", value=log(white), vary=True)
        out = lmfit.minimize(_residuals_white, params, args=(f, log_Pj))
        log_w = params["log_white"].value
        fk = 0.
        alpha = 0.
    w = exp(log_w)
    return fk, w, alpha, out.redchi


def psd_fit_errors(f, psd, white=1., ftol=1e-8):
    log_Pj = log(psd)
    params = lmfit.Parameters()
    params.add("f_knee", value=0.1, min=0.0, max=25.)
    params.add("log_white", value=log(white), vary=True)
    params.add("alpha", value=2., min=0.0, max=20.)
    out = lmfit.minimize(_residuals, params, args=(f, log_Pj), ftol=ftol, scale_covar=True)
    fk = params["f_knee"].value
    log_w = params["log_white"].value
    alpha = params["alpha"].value
    fkstd = params["f_knee"].stderr
    log_wstd = params["log_white"].stderr
    alphastd = params["alpha"].stderr
    if (fk>10.) or (alpha>10.):
        params = lmfit.Parameters()
        params.add("log_white", value=log(white), vary=True)
        out = lmfit.minimize(_residuals_white, params, args=(f, log_Pj))
        log_w = params["log_white"].value
        fk = 0.
        alpha = 0.
        fkstd = 0.
        alphastd = 0.
        log_wstd = params["log_white"].stderr
    w = exp(log_w)
    wstd = abs(log_wstd*w)
    return fk, w, alpha, fkstd, wstd, alphastd


def fullresid(p, f, psd, counter):
    fk = p["f_knee"].value
    w = p["white"].value
    a = p["alpha"].value
    model = w*(1. + (fk/f)**a)
    return (psd - model) / (psd/np.sqrt(counter))

def fullresid_w(p, f, psd, counter):
    w = p["white"].value
    model = w*np.ones(len(f))
    return (psd - model) / (psd/np.sqrt(counter))
    
def psd_fit_errors_all(f, psd, counter, ftol=1e-8):
    winit = get_init(f, psd)
    params = lmfit.Parameters()
    params.add("f_knee", value=0.1, min=0.0, max=25.)
    params.add("white", value=winit, vary=True)
    params.add("alpha", value=2., min=0.0, max=20.)
    out = lmfit.minimize(fullresid, params, args=(f, psd, counter), ftol=ftol, scale_covar=True)
    fk = params["f_knee"].value
    w = params["white"].value
    alpha = params["alpha"].value
    fkstd = params["f_knee"].stderr
    wstd = params["white"].stderr
    alphastd = params["alpha"].stderr
    if (fk>10.) or (alpha>10.):
        params = lmfit.Parameters()
        params.add("white", value=winit, vary=True)
        out = lmfit.minimize(fullresid_w, params, args=(f, psd, counter), ftol=ftol, scale_covar=True)
        w = params["white"].value
        fk = 0.
        alpha = 0.
        fkstd = 0.
        alphastd = 0.
        wstd = params["white"].stderr
    return fk, w, alpha, fkstd, wstd, alphastd

def get_init(f, psd):
    ind = (f>1) & (f<10)
    return np.mean(psd[ind])


