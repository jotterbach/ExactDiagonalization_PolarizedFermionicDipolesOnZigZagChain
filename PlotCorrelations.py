# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 08:59:58 2013

@author: jotterbach
"""

from numpy import *
import matplotlib.pyplot as plt
from ED_HalfFilling import EigSys_HalfFilling
from DotProduct import scalar_prod

srcfile = 'FigureData/Correlations_N6_J1_vdd-50_Theta0'
local_density = load(srcfile+'_local_density.npy')
g2 = load(srcfile+'_g2.npy')

idx = 0
idx_list = list()
for m in srcfile:
    if m == '_':
        idx_list = idx_list + [idx]
    idx += 1        

number_particles = int(srcfile[idx_list[0]+2:idx_list[1]])
if srcfile[idx_list[1]+2] == '-':
    hopping = float(srcfile[idx_list[1]+3:idx_list[2]].replace('-','.'))
    hopping = -hopping
else:
    hopping = float(srcfile[idx_list[1]+2:idx_list[2]].replace('-','.'))

if srcfile[idx_list[2]+4] == '-':
    interaction_strength = float(srcfile[idx_list[2]+5:idx_list[3]].replace('-','.'))
    interaction_strength = -interaction_strength
else:
    interaction_strength = float(srcfile[idx_list[2]+4:idx_list[3]].replace('-','.'))

if srcfile[idx_list[3]+6] == '-':
    angle = float(srcfile[idx_list[3]+7:srcfile.__len__()].replace('-','.'))
    angle = -angle
else:
    angle = float(srcfile[idx_list[3]+6:srcfile.__len__()].replace('-','.'))    
print 'number particles:', number_particles, '\nhopping:', hopping, '\nangle:', angle


vNN = interaction_strength
angle_fraction = angle
gamma = 2*pi/3
v1 = vNN*(1-3*cos(pi - gamma/2 - pi*angle_fraction)**2) / (-.25*(2+6*cos(gamma)))
v2 = vNN*(1-3*cos(gamma/2 - pi*angle_fraction)**2)/ (-.25*(2+6*cos(gamma)))
print 'symmetric nn-int:', (v1+v2)/2, '\nantisymmetric nn-int:', (v1-v2)/2

f101=plt.figure(101)
 #plt.plot(range(0,eigval.shape[0]), eigval.real)
plt.subplot(2,1,1)
plt.plot(range(0,2*number_particles), local_density,color="blue", linewidth=2.5, linestyle='-', marker='s', ms=10, mfc='r')
plt.xlabel(r'site index $i$')
plt.ylabel(r'$\langle n_i \rangle$')
plt.title(r'local density, N='+str(number_particles)+', J='+str(hopping)+' $\Theta/\pi=$'+str(angle))
plt.ylim(-.1, 1.1)
plt.xlim(0,2*number_particles-1)
    
plt.subplot(2,1,2)
plt.plot(range(0,2*number_particles), g2,color="blue", linewidth=2.5, linestyle='-', marker='s', ms=10, mfc='r')
plt.xlabel(r'site index $i$')
plt.ylabel(r'$\langle n_0 n_i \rangle$')
plt.title(r'density-density correlation, N='+str(number_particles)+', J='+str(hopping)+' $\Theta/\pi=$'+str(angle))
plt.ylim(-.1, .6)

#plt.subplot(3,1,3)
#plt.plot(range(0,2*number_particles), sf,color="blue", linewidth=2.5, linestyle='-', marker='s', ms=10, mfc='r')
#plt.xlabel(r'site index $i$')
#plt.ylabel(r'$\langle S^+_0 S^-_i\rangle $')
#plt.title(r'superfluid correlation, N='+str(number_particles)+', J='+str(hopping)+' $\Theta/\pi=$'+str(angle))

plt.autoscale(enable=True, axis='both', tight=True)
plt.tight_layout()



fname = 'Figures/Correlations_N'+str(number_particles)+'_J'+str(hopping).replace('.','-')+'_vdd'+str(interaction_strength).replace('.','-')+'_Theta'+str(angle).replace('.','-')+'.eps'
print fname

f101.savefig(fname)#, dpi=None, facecolor='w', edgecolor='w',
             #orientation='landscape', papertype=None, format=None,
             #transparent=False, bbox_inches=None, pad_inches=0.1)
             
plt.show()