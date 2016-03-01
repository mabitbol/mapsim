import math
import numpy as np
import pylab as pl
import timeit

def pink_noise(length, 
              delta_t = 1e-3,  #[sec]
              net = 250.0,     #[uK sqrt(sec)]
              fknee = 10.0,     #1/f knee [Hz] 
              alpha = 1,     #1/f^alpha exponent
              offset = 0,      #DC offset [uK]
              show_plot = True):

    tic = timeit.default_timer()

    N = length
    nyquist = 1/(2*delta_t)     #Hz
    delta_nu = 1/(N*delta_t)    #Hz
    norm = N/delta_t
    ##Creating frequency array for the fft
    frequency = delta_nu*np.arange(N)
    frequency[N/2+1:] = -1*frequency[N/2-1:0:-1]

    ##Noise model
    whitenoise = net*math.sqrt(2)  #uK /sqrt(Hz) 
    model = np.zeros(N/2+1)
    modelps = np.zeros([2,N])
    modelps[0,:] = frequency
    if fknee == 0:
        model[1:] = whitenoise
        modelps[1,0] = offset**2
        modelps[1,1:] = (whitenoise**2)/2.
        noise = np.sqrt(1./delta_t)*np.random.normal(loc=offset, scale=net, size=N)
    else:
        model[1:] = whitenoise*np.sqrt(1.0+(fknee/frequency[1:N/2+1])**alpha)
        modelps[1,0] = offset**2
        modelps[1,1:] = ((whitenoise*np.sqrt(1.0+(fknee/np.absolute(frequency[1:]))**alpha))**2)/2.
        
        fnoise = np.zeros(N) + 1j*np.zeros(N)
        real = math.sqrt(norm)*model[1:N/2]*np.random.standard_normal(N/2-1)/2.
        imag = math.sqrt(norm)*model[1:N/2]*np.random.standard_normal(N/2-1)/2.
        fnoise[1:N/2] = real + 1j*imag
        fnoise[N/2+1:] = real[::-1] - 1j*imag[::-1]
        ##Setting zero and nyquist frequency elements of the fft by hand
        nyqcomp = math.sqrt(norm)*model[N/2]*np.random.standard_normal(1)/math.sqrt(2)
        fnoise[N/2] = nyqcomp
        fnoise[0] = offset
        ##Noise time stream
        noise = np.real(np.fft.ifft(fnoise))

    ##Power spectrum results
    ps = (1./norm)*abs(np.fft.fft(noise))**2
    ps[1:N/2] = ps[1:N/2] + ps[N:N/2:-1]
    ps = np.sqrt(ps[0:N/2+1])
    freq = frequency[0:N/2+1]
    time = delta_t*np.arange(N)

    toc = timeit.default_timer()

    print "elapsed:", toc-tic

    if show_plot:
        ##Plotting noise time stream
        pl.figure(1)
        pl.xlabel("Time (s)")
        pl.ylabel("$\mu$K")
        pl.title("Noise Realization")
        pl.plot(time[:(100/delta_t)],noise[:(100/delta_t)])
        pl.grid()
        pl.show()
        ##Plotting noise power spectrum
        pl.figure(2)
        pl.xlabel("Frequency (Hz)")
        pl.ylabel("$\mu K^2/Hz$")
        pl.title("Power Spectrum of Noise Realization")
        pl.loglog(freq[1:],ps[1:]**2)
        pl.loglog(freq[1:],model[1:]**2)
        pl.grid()
        pl.show()

    return noise


def pink_noise2(length, 
              delta_t = 1e-3,  #[sec]
              net = 250.0,     #[uK sqrt(sec)]
              fknee = 10.0,     #1/f knee [Hz] 
              alpha = 1,     #1/f^alpha exponent
              offset = 0,      #DC offset [uK]
              show_plot = True,
              power = False):

    tic = timeit.default_timer()

    N = length
    nyquist = 1./(2.*delta_t)     #Hz
    delta_nu = 1./(N*delta_t)    #Hz
    ##Creating frequency array for the fft
    frequency = delta_nu*np.arange(N)
    frequency[N/2+1:] = -1*frequency[N/2-1:0:-1]

    ##Noise model
    whitenoise = net*math.sqrt(2)
    model = np.zeros(N/2+1)
    modelps = np.zeros([2,N])
    modelps[0,:] = frequency
    if fknee == 0:
        model[1:] = whitenoise
        modelps[1,0] = offset**2
        modelps[1,1:] = (whitenoise**2)/2.
    else:
        model[1:] = whitenoise*np.sqrt(1.0+(fknee/frequency[1:N/2+1])**alpha)
        modelps[1,0] = offset**2
        modelps[1,1:] = ((whitenoise*np.sqrt(1.0+(fknee/np.absolute(frequency[1:]))**alpha))**2)/2.

    ##Frequency domain noise
    ##Reverse fft to create time domain noise
    fnoise = np.zeros(N) + 1j*np.zeros(N)
    real = math.sqrt(delta_nu)*model[1:N/2]*np.random.standard_normal(N/2-1)/2.
    imag = math.sqrt(delta_nu)*model[1:N/2]*np.random.standard_normal(N/2-1)/2.
    fnoise[1:N/2] = real + 1j*imag
    fnoise[N/2+1:] = real[::-1] - 1j*imag[::-1]
    ##Setting zero and nyquist frequency elements of the fft by hand
    nyqcomp = math.sqrt(delta_nu)*model[N/2]*np.random.standard_normal(1)/math.sqrt(2)
    fnoise[N/2] = nyqcomp
    fnoise[0] = offset
    ##Noise time stream
    noise = np.real(np.fft.ifft(fnoise))

    ##Power spectrum results
    #Note PS is not the density.  FFT is not normalized by 1/N
    ps = (1./(delta_nu))*abs(np.fft.fft(noise))**2
    ps[1:N/2] = ps[1:N/2] + ps[N:N/2:-1]
    ps = np.sqrt(ps[0:N/2+1])
    freq = frequency[0:N/2+1]
    time = delta_t*np.arange(N)

    toc = timeit.default_timer()

#    print "elapsed:", toc-tic

    if show_plot == True:
        ##Plotting noise time stream
        pl.figure(1)
        pl.xlabel("Time (s)")
        pl.ylabel("$\mu$K")
        pl.title("Noise Realization")
        #pl.plot(time[:(100/delta_t)],noise[:(100/delta_t)])
        pl.plot(time,noise)
        pl.grid()
        pl.show()
        ##Plotting noise power spectrum
        pl.figure(2)
        pl.xlabel("Frequency (Hz)")
        pl.ylabel("$\mu K/\sqrt{Hz}$")
        pl.title("Amplitude Spectrum (sqrt(PS)) of Noise Realization")
        pl.loglog(freq[1:],ps[1:])
        pl.loglog(freq[1:],np.sqrt(modelps[1,1:N/2+1]))
        pl.grid()
        pl.show()

    if power:
        return noise, freq, ps
    else:
        return noise
