#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:45:51 2019

@author: smrak
"""
import os
import satio
from glob import glob
from gpstec import gpstec
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartomap.geogmap as cm

SAVE = 0
odir = '/home/smrak/Documents/gem2019/201709/tec0-50/'

fnlist = [#'/media/smrak/figures/gpstec/2017/gps/04sep17/conv_0904T0000-0905T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/05sep17/conv_0905T0000-0906T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/06sep17/conv_0906T0000-0907T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/07sep17/conv_0907T0000-0908T0000.h5',
          '/media/smrak/figures/gpstec/2017/gps/08sep17/conv_0908T0000-0909T0000.h5',
          '/media/smrak/figures/gpstec/2017/gps/09sep17/conv_0909T0000-0910T0000.h5'
#          '/media/smrak/figures/gpstec/2017/gps/10sep17/conv_0910T0000-0911T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/11sep17/conv_0911T0000-0912T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/12sep17/conv_0912T0000-0913T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/13sep17/conv_0913T0000-0914T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/14sep17/conv_0914T0000-0915T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/15sep17/conv_0915T0000-0916T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/16sep17/conv_0916T0000-0917T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/17sep17/conv_0917T0000-0918T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/18sep17/conv_0918T0000-0919T0000.h5']
]

if 'D' not in locals():
    D = gpstec.merge_time(fnlist)

obstimes = [datetime(2017,9,8,4,45)]
tlim = [datetime(2017,9,15,0), datetime(2017,9,16,0)]
#dt = 60
#obstimes = []
#t = tlim[0]
#while t <= tlim[1]:
#    obstimes.append(t)
#    t += timedelta(minutes=dt)

# DMSP

dmspd = '/home/smrak/Documents/LWS2019/8sep17/dmsp/data/'
svlist = sorted(glob(dmspd + '*.txt'))

dmspt = [datetime(2017,9,8,4,45), datetime(2017,9,8,5,15)]
conjugated = 1

dmspfn = svlist[12]
DMSP = satio.readDMSPtxt(dmspfn)
dt = DMSP.time.values.astype("datetime64[s]")
ts = np.int64(dt)
ddt = np.array([datetime.utcfromtimestamp(t) for t in ts])
idt = (ddt >= dmspt[0]) & (ddt <= dmspt[1])
ts = ts[idt]
ddt = ddt[idt]

glat = DMSP.rpa.sel(data='glat').values[idt]
glon = np.unwrap(DMSP.rpa.sel(data='glon').values[idt], 180)
galt = DMSP.rpa.sel(data='galt').values[idt]
if conjugated:
    glat, glon = satio.magConjugated(glat=glat, glon=glon, galtkm=galt, conj_galtkm=galt)

Ti = DMSP.rpa.sel(data='Ti').values[idt]
Ti[Ti < 0] = np.nan
Ti = satio.runningMedian(Ti, 20)
Te = satio.runningMedian(DMSP.rpa.sel(data='Te').values[idt], 60)
Ni = satio.runningMedian(DMSP.rpa.sel(data='Ni').values[idt], 30)
vx = DMSP.rpa.sel(data='vx').values[idt]
vy = satio.runningMedian(DMSP.rpa.sel(data='vy').values[idt], 60)
vz = DMSP.rpa.sel(data='vz').values[idt]

latlim = [-70, 70]
lonlim = [-170, -30]
tecclim = [0, 30]

for it in obstimes:
#D = gpstec.readFromHDF(TECFN)
    tectime = D['time']
    xgrid = D['xgrid']
    ygrid = D['ygrid']
    idt_tec = abs(tectime - it).argmin()
    tecim = D['tecim'][idt_tec-2 : idt_tec+1]
    tecim = np.nanmean(tecim, 0)
    
    fig, ax = cm.plotCartoMap(latlim=latlim, lonlim=lonlim, projection='stereo',
                              meridians=None, parallels=None,
                              grid_linewidth=1, states = False,
                              title='{} -- {}'.format(os.path.split(dmspfn)[1][:3].upper(), it),
                              background_color='grey',
    #                          nightshade=True, ns_dt=it,
                              apex=True, mlat_levels=[-80,-60,-40,-20,0,20,40,60,80,90],
                              mlat_colors='m', mgrid_width=1, mgrid_style='--',
                              mlon_levels=np.arange(0,361,40), mlat_labels=False,
                              mlon_colors='m', mlon_labels=False,
#                              igrf=True, date=it,
#                              decl_levels=[0], decl_colors='w')
                              )
    
    im = plt.pcolormesh(xgrid, ygrid, tecim.T, cmap='Greys_r',#cmap='nipy_spectral', 
                        transform=ccrs.PlateCarree())
    posn = ax.get_position()
    cax = fig.add_axes([posn.x0+posn.width+0.01, posn.y0, 0.02, posn.height])
    fig.colorbar(im, cax=cax, label='TEC [TECu]')
    im.set_clim(tecclim)
    # Cross-track
    maskleft = (vy < 0)
    ax.scatter(glon[maskleft], glat[maskleft], s=-vy[maskleft], c=vz[maskleft], 
               marker=1, cmap='bwr', 
               vmin=-5e2, vmax=5e2, 
               transform=ccrs.PlateCarree())
    maskright = (vy > 0)
    dmspim = ax.scatter(glon[maskright], glat[maskright], c=vz[maskright], s=vy[maskright], 
               marker=0, cmap='bwr', 
               vmin=-5e2, vmax=5e2, 
               transform=ccrs.PlateCarree())
    posn = ax.get_position()
    cax = fig.add_axes([posn.x0, posn.y0-0.04, posn.width, 0.02])#, orientation='horizontal')
    fig.colorbar(dmspim, cax=cax, label='vz [m/s]', orientation='horizontal')
    im.set_clim(tecclim)
    
    # cross-track legend
    ax.scatter(lonlim[0], latlim[1]-10, s=1000, 
               marker='_', c = 'r', lw=2,
#               vmin=-5e2, vmax=5e2, 
               transform=ccrs.PlateCarree())
    
    # DMSP track
    ax.scatter(glon, glat, c=Ti, vmin=1500, vmax=3000, marker='.', s=2,
               cmap = 'viridis', transform=ccrs.PlateCarree())
    
    if SAVE:
        if not os.path.exists(odir):
            import subprocess
            subprocess.call('mkdir -p {}'.format(odir), shell=True, timeout=2)
        plt.savefig(odir+'{}.png'.format(it.strftime('%m%d_%H%M')), dpi=200)
        plt.close(fig)
