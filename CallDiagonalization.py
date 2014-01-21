# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 19:54:52 2014

@author: jotterbach
"""

import DiagonalizationMethods as dm
import numpy as np
import sys


method_id = str(sys.argv[1])
if method_id != 'corr' and len(sys.argv) == 9:
    number_particles = int(sys.argv[2])
    noEV = int(sys.argv[3])
    gamma = float(sys.argv[4])
    arg4 = float(sys.argv[5])
    range_min = float(sys.argv[6])
    range_max = float(sys.argv[7])
    range_points = int(sys.argv[8])
    arg5 = np.linspace(range_min, range_max, range_points)
elif method_id == 'corr' and len(sys.argv) == 8:
    number_particles = int(sys.argv[2])
    noEV = int(sys.argv[3])
    gamma = float(sys.argv[4])
    angle = float(sys.argv[5])
    hopping = float(sys.argv[6])
    degeneracy = int(sys.argv[7])
else:
    print 'method call or specified parameters not right.'
    print 'following methods available: angle, int, hop, corr'

print str(sys.argv)
if method_id == 'angle':
    "AngleSpectrum(number_particles, noEV, gamma, hopping, angle)"
    dm.AngleSpectrum(number_particles, noEV, gamma*np.pi, arg4, arg5)
elif method_id == 'int':
    "InteractionSpectrum(number_particles, noEV, gamma, angle, interaction_strength)"
    dm.InteractionSpectrum(number_particles, noEV, gamma*np.pi, arg4, arg5)
elif method_id == 'hop':
    "HoppingSpectrum(number_particles, noEV, gamma, angle, hopping)"
    dm.HoppingSpectrum(number_particles, noEV, gamma*np.pi, arg4, arg5)
elif method_id == 'corr':
    dm.DensityCorrelations(number_particles, noEV, gamma, angle, hopping, degeneracy)
else:
    print 'oh snap! something went wrong'
