# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:47:42 2013

@author: jotterbach
"""

from os.path import *
from numpy import *
from scipy.special import *
from ED_HalfFilling import EigSys_HalfFilling
from DotProduct import scalar_prod


sx = array( [ (0, 1), (1, 0) ] )
sy = array( [ (0, -1j), (1j, 0) ] )
sz = array( [ (1, 0), (0, -1) ] )
id2 = array( [ (1, 0), (0 ,1) ] )



def loc_den(basisstates, number_particles, number_states, eigvec, site_idx):

    operator = zeros((2*number_particles,2,2))
    for m in range(0,2*number_particles):
        operator[m,:,:]=id2
    
    operator2 = operator.copy()
    operator2[site_idx,:,:] = real(id2+sz)/2
    local_density = 0
    for m in range(0,number_states):
        local_density = local_density+(dot(eigvec[m].conj(),eigvec[m])*scalar_prod(basisstates[m,:], basisstates[m,:], operator2))
    return site_idx, local_density

def pair_corr(basisstates, number_particles, number_sites, number_states, eigvec, site_idx):
    operator = zeros((number_sites,2,2))
    for m in range(0,number_sites):
        operator[m,:,:]=id2
        
    operator2 = operator.copy()
    operator2[0,:,:] = real(id2+sz)/2
    operator2[site_idx,:,:] = real(id2+sz)/2
    g2 = 0
    for m in range(0,number_states):
        g2 = g2+(dot(eigvec[m].conj(),eigvec[m])* scalar_prod(basisstates[m,:], basisstates[m,:], operator2))
    return site_idx, g2

