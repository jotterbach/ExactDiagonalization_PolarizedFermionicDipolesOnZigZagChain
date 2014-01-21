# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 09:09:36 2013

@author: jotterbach
"""

from numpy import *
import matplotlib.pyplot as plt


angle_full = load('FigureData/Spectrum_N6_J0-1_vdd1_angle.npy')
eigval = load('FigureData/Spectrum_N6_J0-1_vdd1_EigVals.npy')
number_particles = 6
hopping = .1
#angle = angle_full[:,0]
angle = linspace(-1,1,41)
print eigval.shape


f101=plt.figure(101)
noEV_plot=21

plt.subplot(2,1,1)
colorarray = array(['blue', 'red', 'green', 'magenta', 'cyan', 'black'])
for m in range(0,noEV_plot):
    plt.plot(angle, eigval[:,m],'ro')#,color=colorarray[remainder(m,6)], linewidth=2.5, linestyle='-')
plt.xlabel(r'angle $ \theta/\pi $')
plt.ylabel(r'energy $ \omega $')
title1 = 'many-body spectrum, N='+str(number_particles)+', J='+str(hopping)
plt.title(r'many-body spectrum, N='+str(number_particles)+', J='+str(hopping))
plt.autoscale(enable=True, axis='both', tight=True)


plt.subplot(2,1,2)
colorarray = array(['blue', 'red', 'green', 'magenta', 'cyan', 'black'])
for m in range(0,noEV_plot-1):
    plt.plot(angle, abs(eigval[:,m]-eigval[:,m+1]),'ro')#color=colorarray[remainder(m,6)], linewidth=2.5, linestyle='-')
#plt.yscale('log')
plt.xlabel(r'angle $ \theta/\pi $')
plt.ylabel(r'gap $ \Delta $')
plt.title(r'many-body gap, N='+str(number_particles)+', J='+str(hopping))
plt.autoscale(enable=True, axis='both', tight=True)

plt.tight_layout()
plt.show()

fname = 'Figures/Spectrum_N'+str(number_particles)+'_J'+str(hopping).replace('.','-')+'.eps'
print fname

#f101.savefig(fname)#, dpi=None, facecolor='w', edgecolor='w',
             #orientation='landscape', papertype=None, format=None,
             #transparent=False, bbox_inches=None, pad_inches=0.1)