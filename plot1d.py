# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 16:44:35 2019

@author: smrak@bu.edu
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def _runningMedian(x, N):
    shift = int(N/2)
    iterate = np.arange(shift, x.size- shift)
    y = np.nan * np.zeros(x.size)
    for i in iterate:
        y[i] = np.nanmedian(x[i-shift : i+shift])
    return y

def plotIDM(XD, tlim=None, xlabel='mlat,mlt', remove_outliers=False, 
            eps=30, figsize=[10,8], ms=3, title='', nticks=None):
    
    ts = np.int64(XD.time.values.astype("datetime64[s]"))
    dt = np.array([datetime.utcfromtimestamp(t) for t in ts])
    if tlim is None:
        idt = np.ones(dt.size, dtype=bool)
    else:
        idt = (dt >= tlim[0]) & (dt <= tlim[1])
    
    ts = ts[idt]
    dt = dt[idt]
    vx = XD.rpa.sel(data='vx').values[idt]
    vy = XD.rpa.sel(data='vy').values[idt]
    vz = XD.rpa.sel(data='vz').values[idt]
    
    vxm = vx == -9999.
    vym = vy == -9999.
    vzm = vz == -9999.
    vx[vxm] = np.nan
    vy[vym] = np.nan
    vz[vzm] = np.nan
    if remove_outliers:
        vx = _runningMedian(vx, eps)
        vy = _runningMedian(vy, eps)
        vz = _runningMedian(vz, eps)
        
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(311)
    ax1.set_title(title)
    ax1.plot(ts, vx, 'b', ms=ms)
    ax1.set_ylabel('vx [m/s]')
    
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax2.plot(ts, vy, '.b', ms=ms)
    ax2.set_ylabel('vy [m/s]')
    
    ax3 = fig.add_subplot(313, sharex=ax1)
    ax3.plot(ts, vz, '.b', ms=ms)
    ax3.set_ylabel('vz [m/s]')
    
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.setp( ax2.get_xticklabels(), visible=False)
    
    if nticks is not None:
        ticks = np.linspace(0, dt.size-1, nticks).astype(np.int16)
        ax1.set_xticks(ts[ticks])
    
    xti = np.where(np.isin(ts, ax1.get_xticks()))[0]
    xlabels = np.empty(xti.size, dtype='U30')
    for i,ix in enumerate(xti):
        xlabels[i] = dt[ix].strftime('%H:%M') # UT
    for label in xlabel.split(','):
        for i,ix in enumerate(xti):
            arg = XD.rpa.sel(data=label).values[idt]
            xlabels[i] += '\n{}'.format(arg[ix])
    ax3.set_xticklabels(xlabels)
    
    return fig
    
def plotRPA(XD, tlim=None, xlabel='mlat,mlt', remove_outliers=False, 
            eps=30, figsize=[10,8], ms=3, title='', nticks=None,
            hydrogen=False, oxygen=False, helium=False):
    
    ts = np.int64(XD.time.values.astype("datetime64[s]"))
    dt = np.array([datetime.utcfromtimestamp(t) for t in ts])
    if tlim is None:
        idt = np.ones(dt.size, dtype=bool)
    else:
        idt = (dt >= tlim[0]) & (dt <= tlim[1])
    
    ts = ts[idt]
    dt = dt[idt]
    Ni = XD.rpa.sel(data='Ni').values[idt]
    O = XD.rpa.sel(data='O').values[idt]
    H = XD.rpa.sel(data='H').values[idt]
    He = XD.rpa.sel(data='He').values[idt]
    Te = XD.rpa.sel(data='Te').values[idt]
    Ti = XD.rpa.sel(data='Ti').values[idt]
    
    Nim = Ni <= 0
    Om = O == -9999.
    Hm = H == -9999.
    Hem = He == -9999.
    Tim = (Ti == -9999.) | (Ti <= 0)
    Tem = (Te == -9999.) | (Te <= 0)
    
    Ni[Nim] = np.nan
    O[Om] = np.nan
    H[Hm] = np.nan
    He[Hem] = np.nan
    Ti[Tim] = np.nan
    Te[Tem] = np.nan
    
    if remove_outliers:
        Ni = _runningMedian(Ni, eps)
        Ti = _runningMedian(Ti, eps)
        Te = _runningMedian(Te, eps)
    
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(311)
    ax1.set_title(title)
    ax1.plot(ts, Ni, '.b', ms=ms)
    if hydrogen:
        ax1.plot(ts, Ni*H, '.k', ms=1, label='H$^+$')
    if oxygen:
        ax1.plot(ts, Ni*O, '.g', ms=1, label='O$^+$')
    if helium:
        ax1.plot(ts, Ni*He, '.', c='orange', ms=1, label='He$^+$')
    ax1.set_ylabel('Ni [/cm$^3$]')
    try:
        ax1.get_legend_handles_labels()[1][0]
        ax1.legend()
    except:
        pass
    
    
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax2.plot(ts, Ti, '.b', ms=ms)
    ax2.set_ylabel('Ti [K]')
    
    ax3 = fig.add_subplot(313, sharex=ax1)
    ax3.plot(ts, Te, '.b', ms=ms)
    ax3.set_ylabel('Te [K]')
    
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.setp( ax2.get_xticklabels(), visible=False)
    
    if nticks is not None:
        ticks = np.linspace(0, dt.size-1, nticks).astype(np.int16)
        ax1.set_xticks(ts[ticks])
    
    xti = np.where(np.isin(ts, ax1.get_xticks()))[0]
    xlabels = np.empty(xti.size, dtype='U30')
    for i,ix in enumerate(xti):
        xlabels[i] = dt[ix].strftime('%H:%M') # UT
    for label in xlabel.split(','):
        for i,ix in enumerate(xti):
            arg = XD.rpa.sel(data=label).values[idt]
            xlabels[i] += '\n{}'.format(arg[ix])
    ax3.set_xticklabels(xlabels)
    
    return fig