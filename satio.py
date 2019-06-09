# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 14:02:05 2019

@author: smrak@bu.edu
"""
import pysatCDF
import numpy as np
from glob import glob
import os
import xarray
from datetime import datetime, timedelta
from scipy import signal
import matplotlib.pyplot as plt
from apexpy import Apex

def readDMSPtxt(filename = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\Eclipse2017\\data\\Hairston\\data\\',
                save=False, fout=None):

    data = np.loadtxt(filename, skiprows=3)
    yy = int(os.path.split(filename)[1][6:8])
    secinday = data[:,1]
    dt = [datetime.strptime(str(yy) + str(data[i][0])[-5:-2] + str(timedelta(seconds=t)), '%y%j%H:%M:%S') for i,t in enumerate(secinday)]
    t64 = np.array([str(np.datetime64(t)) for t in dt])
    dimension = np.array(['R', 'I', 'galt', 'glat', 'glon', 'mlat', 'mlt', 
                          'vx', 'vy', 'vz', 'RMSx', 'SigmaY', 'SigmaZ', 
                          'Ni', 'O', 'He', 'H', 'Ti', 'Te', 'pts'])
    D = xarray.Dataset({'rpa': (('time', 'data'), data[:,2:])}, coords = {'time': t64, 'data': dimension})
    
    if save:
        if fout is None:
            fout = filename + '.nc'
        D.to_netcdf(fout, mode='w', group='dmsp')
        return
    else:
        return D

def readSWARMascii(folder=None, sv=None):
    fnlist = sorted(glob(folder+sv+'*'))
    data = {}
    dimension = []
    for fn in fnlist:
        fnn = os.path.split(fn)[1]
        name = fnn.rsplit('.')[0].replace(sv+'_', "")#split('_')[1:]
        dimension.append(name)
        if 't1' not in locals():
            t1 = np.genfromtxt(fn, dtype='U30', usecols=0)
            dt = np.array([t.replace('/', 'T') for t in t1], dtype='datetime64[ms]')
            gloc = np.nan * np.ones((dt.size, 2))
        if 'glat' in name:
            gloc[:,0] = np.genfromtxt(fn, dtype=np.float32, usecols=1)
        elif 'glon' in name:
            gloc[:,1] = np.genfromtxt(fn, dtype=np.float32, usecols=1)
        elif '16hz_ne' in name:
            if 't16' not in locals():
                t16 = np.genfromtxt(fn, dtype='U30', usecols=0)
                dt16 = np.array([t.replace('/', 'T') for t in t16], dtype='datetime64[ms]')
            if '16hz_ne_d' in name:
                ned = np.genfromtxt(fn, dtype=np.float32, usecols=1)
            else:
                ne = np.genfromtxt(fn, dtype=np.float32, usecols=1)
    
    D = xarray.Dataset()
    D['t'] = dt
    D['glat'] = gloc[:,0]
    D['glon'] = gloc[:,1]
    D['t16'] = dt16
    D['ne'] = ne
    D['ned'] = ned
    
    return D

def readSWARMcdf(fn):
    typ = os.path.split(fn)[1].split('_')[4]
    with pysatCDF.CDF(fn) as cdf:
        d = cdf.data
    if typ == 'HM':#in ('a', 'b', 'A', 'B'):
        D = {'time': d['Timestamp'],
              'glat': d['Latitude'],
              'glon': d['Longitude'],
              'galt': d['Height'],
              'mlat': d['AACGMLat'],
              'mlon': d['AACGMLon'],
              'Te': d['T_elec'],
              'V': d['U_SC'],
              'Ni': d['n'],
              'flags': d['Flagbits']}
    elif typ == 'FP':#in ('c', 'C'):
        D = {'time': d['Timestamp'],
             'glat': d['Latitude'],
             'glon': d['Longitude'],
             'galt': d['Height'],
             'Ni': d['Density'],
             'j': d['Current']}
    return D

def magConjugated(glat, glon, galtkm, conj_galtkm=None):
    if conj_galtkm is None:
        conj_galtkm = galtkm
    A = Apex()
    mlat, mlon = A.convert(lon=np.unwrap(glon, 180), lat = glat, source='geo', dest='apex', height = galtkm)
    Cglat, Cglon = A.convert(lat=-mlat, lon=mlon, source='apex', dest='geo', height = conj_galtkm)
    
    return Cglat, Cglon

def runningMedian(x, N):
    shift = int(N/2)
    iterate = np.arange(shift, x.size- shift)
    y = np.nan * np.zeros(x.size)
    for i in iterate:
        y[i] = np.nanmedian(x[i-shift : i+shift])
    return y

def butter_hpf(highcut, fs, order):
    """
    Sebastijan Mrak
    Design the Butterwoth response highpass filter with N-th order and 
    3db cutoff frequency 'highcut' in Hz.
    Output are the poles 'b' and zeroes 'a' of the filter
    """
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high, btype='highpass')
    w, h = signal.freqz(b, a, worN=1000)
    return b, a

def butter_lpf(fc, fs, order):
    """
    Sebastijan Mrak
    Design the Butterwoth response highpass filter with N-th order and 
    3db cutoff frequency 'highcut' in Hz.
    Output are the poles 'b' and zeroes 'a' of the filter
    """
    nyq = 0.5 * fs
    high = fc / nyq
    b, a = signal.butter(order, high, btype='lowpass', analog=False)
    w, h = signal.freqz(b, a, worN=1000)
    return b, a

def bpf(y, lowcut, highcut, fs=1, order=5, plot=False):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
#    mid = ((highcut - lowcut)/3)*2 / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    idx = abs(w-high).argmin()
    y_filt = signal.lfilter(b, a, y)
    if plot:
        plt.figure()
        plt.semilogx(w, 20*np.log10(np.abs(h)))
        plt.plot(w[idx], 20*np.log10(np.abs(h[idx])), 'xr')
        plt.ylim([-60,5])
        plt.title('Magnitude-normalized Butterworth filter frequency response')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        ################################################
    
        plt.figure()
        plt.semilogx(w[1:], gd)
        plt.plot(w[idx], gd[idx], 'xr')
        plt.title('LPF group delay')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    return y_filt , gd[idx]

def hpf(y, fc=0.1, order=5, fs=1,plot=False, group_delay=False, verbatim=False):
    """
    Sebastijan Mrak
    Filter the input data 'y' with desired HP filter.  
    """
    b, a = butter_hpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    if verbatim: print ('Group delay of the filter is '+ str(gd[-1])+' samples.')
    if plot:
        plt.figure()
        plt.semilogx(w, 20*np.log10(np.abs(h)))
        plt.ylim([-60,5])
        plt.title('Magnitude-normalized Butterworth filter frequency response')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        ################################################
    
        plt.figure()
        plt.semilogx(w[1:], gd)
        plt.title('LPF group delay')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    if group_delay:
        return y_filt, gd[-1]
    else:
        return y_filt

def lpf(y, fc=0.1, order=5, fs=1, plot=False, group_delay=False, verbatim=False):
    b, a = butter_lpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    if verbatim: print ('Group delay of the filter is '+ str(gd[1])+' samples.')
    if plot:
        plt.figure()
        plt.semilogx(w, np.log10(np.abs(h)))
        plt.ylim([-60,5])

        plt.title('Magnitude-normalized butterworth filter frequency response')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        ################################################

        plt.figure()
        plt.semilogx(w[1:], gd)
        plt.title('LPF group delay')
        plt.xlabel('Frequency [radians / second]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    if group_delay:
        return y_filt, gd[1]
    else:
        return y_filt