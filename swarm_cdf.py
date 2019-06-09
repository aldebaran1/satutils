#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:54:32 2019

@author: smrak
"""

import satio
import plot1d
import plot2d
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

fn = '/home/smrak/Documents/LWS2019/28may17/swarm/SW_EXTD_EFIA_LP_HM_20170528T000000_20170528T235959_0101.cdf'
tlim = [datetime(2017,5,28,0,50), datetime(2017,5,28,2)]

D = satio.readSWARMcdf(fn)
time = D['time']
idt = (time >= tlim[0]) & (time <= tlim[1])
dt = time[idt].astype("datetime64[ms]")
glat = D['glat'][idt]
glon = D['glon'][idt]
galt = D['galt'][idt]
Ni = D['Ni'][idt]
Ni_d = satio.hpf(D['Ni'], fs=2, order=6, fc=0.1)
Ni_d = Ni_d[idt]

ts = np.int64(dt)

plot1d.plotSWARM(D, tlim=tlim,nticks=6, xlabel=['mlat','mlon'])
plot2d.plotDMSPtrack(glat=glat, glon=glon, galt=galt, nticks=8, dt=dt, tick_labels=True,
                     northb=1, original=0,
                     mlat_labels=0, mlon_labels=0)