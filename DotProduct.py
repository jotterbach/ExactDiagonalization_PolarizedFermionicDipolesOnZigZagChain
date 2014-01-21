# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 13:21:38 2013

@author: jotterbach
"""

from numpy import *
from scipy.special import *
import scipy.sparse as sp
import scipy.sparse.linalg as la
import itertools

def scalar_prod(state1, state2, operator):
    '''calculate the scalar product of many-body states'''
    
    number_sites = state1.shape[0]
    
    local_ev1 = empty((number_sites,2))
    local_ev2 = empty((number_sites,2))
    for m in range(0,number_sites):
        if state1[m]==0:
            local_ev1[m,:] = array([0,1])
        elif state1[m]==1:
            local_ev1[m,:] = array([1,0])
        if state2[m]==0:
            local_ev2[m,:] = array([0,1])
        elif state2[m]==1:
            local_ev2[m,:] = array([1,0])
            
    expv=1;
    for m in range(0,number_sites):
        expv = expv * local_ev1[m,:].dot(dot(operator[m,:,:],local_ev2[m,:].transpose()))
        
    return expv
    