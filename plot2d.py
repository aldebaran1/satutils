# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 22:05:41 2019

@author: smrak@bu.edu
"""
import numpy as np
import datetime
from apexpy import Apex
import matplotlib.pyplot as plt
from cartomap import geogmap as gm
import cartopy.crs as ccrs

def plotDMSPtrack(glon, glat, galt, dt=None,
                  latlim=[-70, 70], lonlim=[-150, -50],
                  meridians = None, 
                  parallels = None,
                  mlon_levels = np.arange(0,361,40), 
                  mlat_levels = np.arange(-90,91,30),
                  original = True, northb = False,
                  conjugated = False, conjugated_altkm = None,
                  conjugated_linestyle = '--',
                  lw=3, color='m', nticks=None, tick_labels=False,
                  mlon_labels=True,
                  mlat_labels=True):
    
    fig, ax = gm.plotCartoMap(latlim=[-70, 70], lonlim=[-150, -50], projection='stereo',
                          background_color='w', grid_linewidth=1,states=False,
                          meridians=meridians, 
                          parallels=parallels,
                          apex=True, mlon_cs='mlon',
                          mlon_levels=mlon_levels,
                          mlon_colors='b', mgrid_style='-', mgrid_width=1,
                          mlat_levels=mlat_levels,
                          mlat_colors='b', 
                          mlon_labels=mlon_labels,
                          mlat_labels=mlat_labels)
    if conjugated_altkm is None:
        conjugated_altkm = galt
    if original:
        if np.nanmax(abs(glon)) > 180:
            glon = np.unwrap(glon, 180)
        plt.plot(glon, glat, c=color, lw=lw, transform=ccrs.PlateCarree())
            
    if conjugated:
        
        #1: To mag
        A = Apex()
        mlat,mlon = A.convert(lon=np.unwrap(glon,180), lat=glat, source='geo', dest='apex', height = conjugated_altkm)
        #2: Conjugated --> back to geog
        Cglat, Cglon = A.convert(lat=-mlat, lon=mlon, source='apex', dest='geo', height = conjugated_altkm)
        #3 Plot conjugaed track
        plt.plot(Cglon, Cglat, conjugated_linestyle, c=color, lw=lw, transform=ccrs.PlateCarree())
    
    if northb:
        glon = np.unwrap(glon, 180)
        #1: To mag
        A = Apex()
        mlat, mlon = A.convert(lon=glon, lat=glat, source='geo', dest='apex', height = conjugated_altkm)
        #2: Conjugated --> back to geog
        Cglat, Cglon = A.convert(lat=-mlat, lon=mlon, source='apex', dest='geo', height = conjugated_altkm)
#        m1 = mlat >= 0
        m2 = mlat <= 0
        glon[m2] = Cglon[m2]
        glat[m2] = Cglat[m2]
        plt.plot(glon, glat, c=color, lw=lw, transform=ccrs.PlateCarree())
#        plt.plot(Cglon[m2], Cglat[m2], c=color, lw=lw, transform=ccrs.PlateCarree())
    
    if nticks is not None:
        ixticks = np.linspace(0, glon.size - 1, nticks).astype(np.int16)
        if original or northb:
            plt.scatter(np.unwrap(glon, 180)[ixticks], glat[ixticks], marker=1, c=color, s=50,
                        transform=ccrs.PlateCarree())
        if conjugated:
            plt.scatter(np.unwrap(Cglon, 180)[ixticks], Cglat[ixticks], marker=1, c=color, s=50,
                        transform=ccrs.PlateCarree())
        if tick_labels:
            if dt is None: 
                print ('Ther is trajectory time information. Set-up the "dt" argument if you want tick labels.')
            else:
                if not isinstance(dt[0], datetime.datetime):
                    dt = dt.astype(datetime.datetime)
                assert isinstance(dt[0], datetime.datetime)
                fmt = np.array([t.strftime("%H:%M") for t in dt[ixticks]])
                if original or northb:
                    fx = np.unwrap(glon, 180)[ixticks]+2.5
                    fy = glat[ixticks]-1
                    for i,txt in enumerate(fmt):
                        if (fx[i] > lonlim[0]) and (fx[i] < lonlim[1]) and (fy[i] > latlim[0]) and (fy[i] < latlim[1]):
                            plt.text(fx[i], fy[i], txt, color=color, transform=ccrs.PlateCarree())
                if conjugated:
                    fx = np.unwrap(Cglon, 180)[ixticks]+2.5
                    fy = Cglat[ixticks]-1
                    for i,txt in enumerate(fmt):
                        if (fx[i] > lonlim[0]) and (fx[i] < lonlim[1]) and (fy[i] > latlim[0]) and (fy[i] < latlim[1]):
                            plt.text(fx[i], fy[i], txt, color=color, transform=ccrs.PlateCarree())
                            
    return fig, ax