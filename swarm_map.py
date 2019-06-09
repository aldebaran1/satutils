#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:55:27 2019

@author: smrak
"""
import os
from apexpy import Apex
from gpstec import gpstec
import satio
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartomap.geogmap as cm

def _runningMax(x, N):
    n2 = int(N/2)
    iterate = np.arange(n2, x.size-n2)
    y = np.nan * np.copy(x)
    for i in iterate:
        chunk = x[i-n2:i+n2]
        y[i] = np.nanmax(abs(chunk))
    return y

SAVE = 1
sv = 'a'
fn = '/home/smrak/Documents/LWS2019/28may17/swarm/SW_EXTD_EFIA_LP_HM_20170528T000000_20170528T235959_0101.cdf'
fn = '/home/smrak/Documents/LWS2019/28may17/swarm/SW_EXTD_EFIC_LP_FP_20170528T000058_20170528T235959_0201.cdf'
fn = '/home/smrak/Documents/LWS2019/28may17/swarm/SW_EXTD_EFIC_LP_HM_20170527T000000_20170527T235959_0101.cdf'
fn = '/home/smrak/Documents/LWS2019/8sep17/swarm/SW_EXTD_EFIA_LP_HM_20170907T000000_20170907T235959_0101.cdf'
#fn = '/home/smrak/Documents/LWS2019/8sep17/swarm/SW_EXTD_EFIA_LP_HM_20170908T000000_20170908T235959_0101.cdf'
#fn = '/home/smrak/Documents/LWS2019/8sep17/swarm/SW_EXTD_EFIC_LP_HM_20170907T000000_20170907T235959_0101.cdf'
#fn = '/home/smrak/Documents/LWS2019/8sep17/swarm/SW_EXTD_EFIC_LP_HM_20170908T000000_20170908T235959_0101.cdf'
odir = '/home/smrak/Documents/gem2019/20170907/swarm{}/'.format(sv)
#tlim = [datetime(2017,5,27,23,20),
#        datetime(2017,5,28,0,0)]
#it = datetime(2017,5,27,23,45)
tlim = [datetime(2017,9,7,0,0),
        datetime(2017,9,7,6,0)]
it = datetime(2017,9,7,2,0)

typ = os.path.split(fn)[1].split('_')[4]
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
if typ == 'HM':# in('a', 'b'):
    Te = D['Te'][idt]
else:
    j = D['j'][idt] /1e3

envelope = _runningMax(abs(Ni_d), 10)
#plt.figure()
#plt.semilogy(dt, envelope, 'b')
#plt.plot(dt, Ni_d, 'b')
#plt.plot(dt, envelope, 'r')
#plt.plot(dt, -envelope, 'r')

latlim = [-70, 70]
lonlim = [-170, -30]
lonlim = [-150, -10]
tecclim = [0, 20]

latlim1 = [-10, 70]
lonlim1 = [-150, -50]

tecfn = ['/media/smrak/figures/gpstec/2017/gps/27may17/conv_0527T0000-0528T0000.h5',
         '/media/smrak/figures/gpstec/2017/gps/28may17/conv_0528T0000-0529T0000.h5',]
tecfn = ['/media/smrak/figures/gpstec/2017/gps/07sep17/conv_0907T0000-0908T0000.h5',
         '/media/smrak/figures/gpstec/2017/gps/08sep17/conv_0908T0000-0909T0000.h5',]

tecrange = 2
Dt = gpstec.merge_time(tecfn)
tectime = Dt['time']
xgrid = Dt['xgrid']
ygrid = Dt['ygrid']
idt_tec = abs(tectime - it).argmin()
tecim = Dt['tecim'][idt_tec - tecrange : idt_tec+tecrange-1]
tecim = np.nanmean(tecim, 0)

fig = plt.figure(figsize=[12,7])
ax0 = plt.subplot(121, projection=ccrs.Stereographic(central_longitude=(sum(lonlim)/2)))
ax1 = plt.subplot(122, projection=ccrs.Stereographic(central_longitude=(sum(lonlim)/2)))
ax0 = cm.plotCartoMap(latlim=latlim, lonlim=lonlim, projection='stereo',
                          meridians=None, parallels=None, ax=ax0,
                          grid_linewidth=1, states = False,
                          title='TEC: {} -- \n {}'.format(tectime[idt_tec - tecrange], tectime[idt_tec + tecrange - 1]),
                          background_color='grey',
                          apex=True, mlat_levels=[-80,-60,-40,-20,0,20,40,60,80,90],
                          mlat_colors='w', mgrid_width=0.5, mgrid_style='--',
                          mlon_levels=np.arange(0,361,40), mlat_labels=False,
                          mlon_colors='w', mlon_labels=False,
                          )
ax1 = cm.plotCartoMap(latlim=latlim1, lonlim=lonlim1, projection='stereo',
                          meridians=None, parallels=None, ax=ax1,
                          grid_linewidth=1, states = False,
                          title='Swarm{}: {} -- \n {}'.format(sv.upper(),tlim[0], tlim[1]),
                          background_color='grey',
                          apex=True, mlat_levels=[-80,-60,-40,-20,0,20,40,60,80,90],
                          mlat_colors='w', mgrid_width=0.5, mgrid_style='--',
                          mlon_levels=np.arange(0,361,40), mlat_labels=False,
                          mlon_colors='w', mlon_labels=False,
                          )

im = ax0.pcolormesh(xgrid, ygrid, tecim.T, cmap='Greys_r',#cmap='nipy_spectral', 
                    transform=ccrs.PlateCarree())
im.set_clim(tecclim)
imne = ax0.scatter(glon-0.5, glat, #s=-vy[maskleft], 
           c=np.log10(Ni), s = envelope / 1e2,
           marker=0, cmap='autumn', 
           vmin=4, vmax=6, 
           transform=ccrs.PlateCarree())
posn = ax0.get_position()
cax2 = fig.add_axes([posn.x0, posn.y0-0.04, posn.width, 0.02])#, orientation='horizontal')
fig.colorbar(imne, cax=cax2, label='Ni [/cm$^3$]', orientation='horizontal')
if typ == 'HM':
    temask = (Te>500) & (Te<10000)
    imte = ax0.scatter(glon[temask]+0.5, glat[temask], #s=-vy[maskleft], 
               c=Te[temask], s = 50,
               marker=1, cmap='cool', 
               vmin=2000, vmax=5e3, 
               transform=ccrs.PlateCarree())
else:
    imte = ax0.scatter(glon + 0.5, glat, #s=-vy[maskleft], 
           c=j, s = 50,
           marker=1, cmap='cool', 
           vmin=-50, vmax=0, 
           transform=ccrs.PlateCarree())
##############################################################################
im11 = ax1.pcolormesh(xgrid, ygrid, tecim.T, cmap='Greys_r',
                    transform=ccrs.PlateCarree())
A = Apex()
mlat, mlon = A.convert(lon=glon, lat=glat, source='geo', dest='apex',)
#2: Conjugated --> back to geog
Cglat, Cglon = A.convert(lat=-mlat, lon=mlon, source='apex', dest='geo',)
m2 = mlat <= 0
imn = ax1.scatter(glon[~m2], glat[~m2], #s=-vy[maskleft], 
           c=np.log10(Ni[~m2]), s = envelope[~m2] / 1e2,
           marker=0, cmap='autumn', 
           vmin=4, vmax=6, 
           transform=ccrs.PlateCarree())
ims = ax1.scatter(Cglon[m2], Cglat[m2],
           c=np.log10(Ni[m2]), s = envelope[m2] / 1e2,
           marker=1, cmap='autumn', 
           vmin=4, vmax=6, 
           transform=ccrs.PlateCarree())
posn = ax1.get_position()
if typ == 'HM':
    temaskn = (Te[~m2]>500) & (Te[~m2]<10000)
    temasks = (Te[m2]>500) & (Te[m2]<10000)
    imnt = ax1.scatter(glon[~m2][temaskn], glat[~m2][temaskn], #s=-vy[maskleft], 
           c=Te[~m2][temaskn]/1e3, s = 2,
           marker=0, cmap='cool', 
           vmin=2, vmax=5, 
           transform=ccrs.PlateCarree())
    imst = ax1.scatter(Cglon[m2][temasks], Cglat[m2][temasks],
           c=Te[m2][temasks]/1e3, s = 2,
           marker=0, cmap='cool', 
           vmin=2, vmax=5, 
           transform=ccrs.PlateCarree())
    cax3 = fig.add_axes([posn.x0, posn.y0-0.04, posn.width, 0.02])
    fig.colorbar(imte, cax=cax3, label='Te [K]', orientation='horizontal')
else:
    imnt = ax1.scatter(glon[~m2], glat[~m2], #s=-vy[maskleft], 
           c=j[~m2], s = 2,
           marker=0, cmap='cool', 
           vmin=-50, vmax=0, 
           transform=ccrs.PlateCarree())
    imst = ax1.scatter(Cglon[m2], Cglat[m2],
           c=j[m2], s = 2,
           marker=0, cmap='cool', 
           vmin=-50, vmax=0, 
           transform=ccrs.PlateCarree())
    cax3 = fig.add_axes([posn.x0, posn.y0-0.04, posn.width, 0.02])
    fig.colorbar(imte, cax=cax3, label='I [uA]', orientation='horizontal')
cax = fig.add_axes([posn.x0+posn.width+0.01, posn.y0, 0.02, posn.height])
fig.colorbar(im, cax=cax, label='TEC [TECu]')
im11.set_clim(tecclim)



if SAVE:
    if not os.path.exists(odir):
        import subprocess
        subprocess.call('mkdir -p {}'.format(odir), shell=True, timeout=2)
    plt.savefig(odir+'{}.png'.format(it.strftime('%m%d_%H%M')), dpi=200)
    plt.close(fig)