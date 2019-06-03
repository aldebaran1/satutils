# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 14:02:05 2019

@author: smrak@bu.edu
"""

import numpy as np
from glob import glob
import os
import xarray
from datetime import datetime, timedelta
from apexpy import Apex

def readDMSPtxt(filename = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\Eclipse2017\\data\\Hairston\\data\\',
                save=False, fout=None):

    data = np.loadtxt(filename, skiprows=3)
    yydd = int(data[0][0])
    secinday = data[:,1]
    dt = [datetime.strptime(str(yydd)+str(timedelta(seconds=t)), '%Y%j%H:%M:%S') for t in secinday]
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

def magConjugated(glat, glon, galtkm, conj_galtkm=None):
    if conj_galtkm is None:
        conj_galtkm = galtkm
    A = Apex()
    mlat, mlon = A.convert(lon=np.unwrap(glon, 180), lat = glat, source='geo', dest='apex', height = galtkm)
    Cglat, Cglon = A.convert(lat=-mlat, lon=mlon, source='apex', dest='geo', height = conj_galtkm)
    
    return Cglat, Cglon