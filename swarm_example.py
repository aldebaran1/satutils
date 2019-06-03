# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:13:07 2019

@author: smrak@bu.edu
"""

import satio
import plot2d
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
folder = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\LWS2019\\28may17\\swarm_ascii\\'
sv = 'swarmA'
tlim = np.array([datetime(2017,5,28,2), datetime(2017,5,28,3)], dtype='datetime64[ns]')
D = satio.readSWARMascii(folder=folder, sv=sv)
idt = (D.t.values >= tlim[0]) & (D.t.values <= tlim[1])
glon = D.glon.values[idt]
glat = D.glat.values[idt]
plot2d.plotDMSPtrack(glat=glat, glon=glon, galt=350, conjugated=True, original=False)
