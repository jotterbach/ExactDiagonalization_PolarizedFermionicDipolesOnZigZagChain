# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:21:19 2013

@author: jotterbach
"""

from numpy import *
from scipy.special import *
import scipy.sparse as sp
import scipy.sparse.linalg as la
import itertools
import time
import re

def Find(pat, text):
    match = re.search(pat, text)
    if match:
#        print match.group()
        return match.group()
    else:
        print 'not found'


def nn_hopping(basis, number_particles, number_sites):
    '''create hopping Hamiltonian for nearest neighbor hopping'''
    
    b = basis
    number_states = b.shape[0]
    number_sites = b.shape[1]


    '''create hashtags '''
    idx = linspace(1,number_sites,number_sites)
    T = zeros((number_states,1))
    print 'calculating the hashes'
    for m in range(number_states):
        T[m] = int((Find(r'\d[.\s\d]+\d',str(b[m,:]))).replace(r' ','').replace('.',''),2)
        print T[m]
    
    ham = sp.lil_matrix((number_states, number_states), dtype=float)
    tic1 = time.time()
    for m in range(0,number_states):
        print m
        for n in range(0,number_sites):
            right = b[m,:].copy()
            left = b[m,:].copy()
            if b[m,mod(n,number_sites)]==1 and b[m,mod(n+1,number_sites)]==0:
                right[mod(n,number_sites)] = 0
                right[mod(n+1,number_sites)] = 1
#                print b[m,:], right
                T_right = int((Find(r'\d[.\s\d]+\d',str(right))).replace(r' ','').replace('.',''),2)
                p_right = argwhere(T == T_right).flatten()
#                print p_right
                ham[m,p_right[0]] = -1
            if b[m,mod(n,number_sites)]==1 and b[m,mod(n-1,number_sites)]==0:
                left[mod(n,number_sites)] = 0
                left[mod(n-1,number_sites)] = 1
#                print b[m,:], left
                T_left = int((Find(r'\d[.\s\d]+\d',str(left))).replace(r' ','').replace('.',''),2)
                p_left = argwhere(T == T_left).flatten()
                ham[m,p_left[0]] = -1
    toc1 = time.time()
    print "generation time with hashtags:", toc1-tic1

    return ham
    
    
def nnn_hopping(basis, number_particles, number_sites):
    '''create next-nearest-neighbor hopping Hamiltonian for nearest neighbor hopping'''
    
    b = basis
    number_states = b.shape[0]
    
    '''create hashtags '''
    idx = linspace(1,number_sites,number_sites)
    T = zeros((number_states,1))
    print 'calculating the hashes'
    for m in range(number_states):
        T[m] = int((Find(r'\d[.\s\d]+\d',str(b[m,:]))).replace(r' ','').replace('.',''),2)
        print T[m]
            
    ham = sp.lil_matrix((number_states, number_states), dtype=float)
    tic1 = time.time()
    for m in range(0,number_states):
        print m
        for n in range(0,number_sites):
            right = b[m,:].copy()
            left = b[m,:].copy()
            if b[m,mod(n,number_sites)]==1 and b[m,mod(n+2,number_sites)]==0:
                right[mod(n,number_sites)] = 0
                right[mod(n+2,number_sites)] = 1
                T_right = int((Find(r'\d[.\s\d]+\d',str(right))).replace(r' ','').replace('.',''),2)
                p_right = argwhere(T == T_right).flatten()
                ham[m,p_right[0]] = -1
            if b[m,mod(n,number_sites)]==1 and b[m,mod(n-2,number_sites)]==0:
                left[mod(n,number_sites)] = 0
                left[mod(n-2,number_sites)] = 1
                T_left = int((Find(r'\d[.\s\d]+\d',str(left))).replace(r' ','').replace('.',''),2)
                p_left = argwhere(T == T_left).flatten()
                ham[m,p_left[0]] = -1
    toc1 = time.time()
    print "generation time with hashtags:", toc1-tic1
                        
    return ham
    
def int_ham(basis, number_particles, number_sites, gamma, vNN, vNNN, angle_fraction, idx_state, independent_v1_v2):
    if independent_v1_v2 ==1:
        v1 = vNN*(cos(pi*angle_fraction)+sin(pi*angle_fraction))
        v2 = vNN*(cos(pi*angle_fraction)-sin(pi*angle_fraction))
        v3 = vNNN
    else:
        v1 = vNN*(1-3*cos(pi - gamma/2 - pi*angle_fraction)**2) / (-.25*(2+6*cos(gamma)))
        v2 = vNN*(1-3*cos(gamma/2 - pi*angle_fraction)**2)/ (-.25*(2+6*cos(gamma)))
        v3 = vNNN*(1-3*cos(pi/2 - pi*angle_fraction)**2)/sqrt(2*(1-cos(gamma)))/ (-.25*(2+6*cos(gamma))) 
    
    v = 0
    idx = basis[idx_state,:].nonzero()[0]
    for m in range(0,idx.shape[0]):
        if (((idx[remainder(m+1,idx.shape[0])]-idx[m]==1) or
            (remainder(idx[remainder(m+1,idx.shape[0])]-idx[m],number_sites)==1)) and
            (remainder(idx[remainder(m+1,idx.shape[0])],2)==1)):
                v = v+v1
        if (((idx[remainder(m+1,idx.shape[0])]-idx[m]==1) or
            (remainder(idx[remainder(m+1,idx.shape[0])]-idx[m],number_sites)==1)) and
            (remainder(idx[remainder(m+1,idx.shape[0])],2)==0)):
                v = v+v2
        if ((idx[remainder(m+1,idx.shape[0])]-idx[m]==2) or
            (remainder(idx[remainder(m+1,idx.shape[0])]-idx[m],number_sites)==2)):
                v = v+v3
    return idx_state, v
        
        
