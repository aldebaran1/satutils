#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:37:19 2019

@author: smrak
"""
import os
import numpy as np
from gpstec import gpstec
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartomap.geogmap as cm

fnlist = [#'/media/smrak/figures/gpstec/2017/gps/06sep17/conv_0906T0000-0907T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/07sep17/conv_0907T0000-0908T0000.h5',
          '/media/smrak/figures/gpstec/2017/gps/08sep17/conv_0908T0000-0909T0000.h5',]
#          '/media/smrak/figures/gpstec/2017/gps/09sep17/conv_0909T0000-0910T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/10sep17/conv_0910T0000-0911T0000.h5',]

#fnlist = ['/media/smrak/figures/gpstec/2017/gps/26may17/conv_0526T0000-0527T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/27may17/conv_0527T0000-0528T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/28may17/conv_0528T0000-0529T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/29may17/conv_0529T0000-0530T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/30may17/conv_0530T0000-0531T0000.h5',
#          '/media/smrak/figures/gpstec/2017/gps/31may17/conv_0531T0000-0601T0000.h5',]

SAVE = 1
odir = '/home/smrak/Documents/gem2019/201709/jro/'
D = gpstec.merge_time(fnlist)

tlim = [datetime(2017,5,26,1), datetime(2017,6,1,0)]
tlim = [datetime(2017,9,8,0), datetime(2017,9,8,10)]
dt = 5
obstimes = []
t = tlim[0]
while t <= tlim[1]:
    obstimes.append(t)
    t += timedelta(minutes=dt)
    
latlim = [-40, 40]
lonlim = [-130, -30]
tecclim = [0, 50]

for it in obstimes:
#D = gpstec.readFromHDF(TECFN)
    tectime = D['time']
    xgrid = D['xgrid']
    ygrid = D['ygrid']
    idt_tec = abs(tectime - it).argmin()
    tecim = D['tecim'][idt_tec-1 : idt_tec]
    tecim = np.nanmean(tecim, 0)
    
    fig, ax = cm.plotCartoMap(latlim=latlim, lonlim=lonlim, projection='stereo',
                              meridians=None, parallels=None,figsize=(10,8),
                              grid_linewidth=1, states = False,
                              title='{}'.format(it), 
                              background_color='grey',
    #                          nightshade=True, ns_dt=it,
                              apex=True, mlat_levels=[-90,-60,-40,-20,0,20,40,60,80,90],
                              mlat_colors='w', mgrid_width=1, mgrid_style='--',
                              mlon_levels=np.arange(0,361,40), mlat_labels=False,
                              mlon_colors='w', mlon_labels=False,
#                              igrf=True, date=it,
#                              decl_levels=[0], decl_colors='w')
                              )
    
    im = plt.pcolormesh(xgrid, ygrid, tecim.T, cmap='nipy_spectral',#cmap='nipy_spectral', 
                        transform=ccrs.PlateCarree())
    posn = ax.get_position()
    cax = fig.add_axes([posn.x0+posn.width+0.01, posn.y0, 0.02, posn.height])
    fig.colorbar(im, cax=cax, label='TEC [TECu]')
    im.set_clim(tecclim)
    
#    plt.tight_layout()
    
    if SAVE:
        if not os.path.exists(odir):
            import subprocess
            subprocess.call('mkdir -p {}'.format(odir), shell=True, timeout=2)
        plt.savefig(odir+'{}.png'.format(it.strftime('%m%d_%H%M')), dpi=100)
        plt.close(fig)
    
#    break