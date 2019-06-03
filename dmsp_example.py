# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 14:04:36 2019

@author: smrak@bu.edu
"""

import satio
import plot1d
import plot2d
from apexpy import Apex
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from cartomap import geogmap as gm
import cartopy.crs as ccrs


def _runningMedian(x, N):
    shift = int(N/2)
    iterate = np.arange(shift, x.size- shift)
    y = np.nan * np.zeros(x.size)
    for i in iterate:
        y[i] = np.nanmedian(x[i-shift : i+shift])
    return y
fn = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\LWS2019\\08sep17\\Sebastijan_8Sept2017\\data\\f18_rl172510149.txt'
tlim = [datetime(2017,9,8,2,30), datetime(2017,9,8,3,30)]
d = satio.readDMSPtxt(fn)
dt = d.time.values.astype("datetime64[s]")
ts = np.int64(dt)
ddt = np.array([datetime.utcfromtimestamp(t) for t in ts])

idt = (ddt >= tlim[0]) & (ddt <= tlim[1])
ts = ts[idt]
ddt = ddt[idt]
mlt = d.rpa.sel(data='mlt').values[idt]
mlat = d.rpa.sel(data='mlat').values[idt]
Ni = d.rpa.sel(data='Ni').values[idt]
O = d.rpa.sel(data='O').values[idt]
I = d.rpa.sel(data='I').values[idt]
R = d.rpa.sel(data='R').values[idt]
vx = d.rpa.sel(data='vx').values[idt]
vy = d.rpa.sel(data='vy').values[idt]
vz = d.rpa.sel(data='vz').values[idt]

glat = d.rpa.sel(data='glat').values[idt]
glon = d.rpa.sel(data='glon').values[idt]
galt = d.rpa.sel(data='galt').values[idt]
#plot1d.plotIDM(XD=d, tlim=tlim, title='F18 -- IDM', nticks=5)
#plot1d.plotRPA(XD=d, tlim=tlim, title='F18 -- RPA', nticks=5, oxygen=0)

plot2d.plotDMSPtrack(glon, glat, galt, conjugated=True, nticks=20, original=False, dt=ddt,tick_labels=True)
#
#fig, ax = gm.plotCartoMap(latlim=[-70, 70], lonlim=[-150, -50], projection='stereo',
#                          background_color='w', grid_linewidth=1,states=False,
#                          meridians=np.arange(0,361,30), 
#                          parallels=np.arange(-90,91,30),
#                          apex=True, mlon_cs='mlon', date=tlim[1], 
#                          mlon_levels=np.arange(0,361,40),
#                          mlon_colors='b', mgrid_style='-', mgrid_width=1,
#                          mlat_levels=np.arange(-90,91,30),
#                          mlat_colors='b', 
#                          mlon_labels=True,
#                          mlat_labels=True)
#
#plt.plot(np.unwrap(glon,180), glat, 'm', lw=3, transform=ccrs.PlateCarree())
##1: 
#A = Apex()
#mlat,mlon = A.convert(lon=np.unwrap(glon,180), lat=glat, source='geo', dest='apex', height = galt)
##2:
#Cglat, Cglon = A.convert(lat=-mlat,lon=mlon,source='apex',dest='geo',height = galt)
#plt.plot(Cglon, Cglat, '--m', lw=3, transform=ccrs.PlateCarree())
#plt.plot(-80,50, 'xr', ms=10, lw=5, transform=ccrs.PlateCarree())
#1: 
#A = Apex()
#mlat,mlon = A.convert(lon=-80, lat=50, source='geo', dest='apex', height = 0)
#Cglat, Cglon = A.convert(lat=-mlat,lon=mlon,source='apex',dest='geo',height=0)
#plt.plot(Cglon,Cglat, 'xr', ms=10, lw=5, transform=ccrs.PlateCarree())

#plt.plot(np.unwrap(glon,180), glat, 'r', lw=3, transform=ccrs.PlateCarree())