# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 09:41:25 2013

@author: jotterbach
"""

from os.path import *
from numpy import *
from scipy.special import *
import scipy.sparse as sp
import scipy.linalg as la_ns
import scipy.sparse.linalg as la
import itertools
from hamiltonian_generator import nn_hopping
from hamiltonian_generator import nnn_hopping
from hamiltonian_generator import int_ham
import time
from multiprocessing import Process
from multiprocessing import *
from multiprocessing import Pool
from ParallelizationTools import info

def EigSys_HalfFilling(number_particles, number_sites, J, vNN, angle_fraction, noEV, spectrum, gamma, independent_v1_v2):
    ''' Hopping matrix element J and angle angle_fraction as fraction of pi
        gamma is angle of lattice, i.e. gamma = 2*pi/3 = 120° for single chain of 
        hexagonal lattice'''
    
    info('Diagonalization Module')    
    
    set_printoptions(precision=4, suppress=True, linewidth = 100)
    
    vNNN = 0*vNN
    if independent_v1_v2 ==1:
        v1 = vNN*(cos(pi*angle_fraction)+sin(pi*angle_fraction))
        v2 = vNN*(cos(pi*angle_fraction)-sin(pi*angle_fraction))
        v3 = vNNN
    else:
        v1 = vNN*(1-3*cos(pi - gamma/2 - pi*angle_fraction)**2) / (-.25*(2+6*cos(gamma)))
        v2 = vNN*(1-3*cos(gamma/2 - pi*angle_fraction)**2)/ (-.25*(2+6*cos(gamma)))
        v3 = vNNN*(1-3*cos(pi/2 - pi*angle_fraction)**2)/sqrt(2*(1-cos(gamma)))/ (-.25*(2+6*cos(gamma)))
    vp = (v1+v2)/2
    vm = (v1-v2)/2
    
    print ' symmetric nn-int: %s \n antisymmetric nn-int: %s \n nnn-int: %s \n hopping: %s' % ( vp, vm, v3, J)
    


    basis_name = 'HamiltonianData/BasisStates_N'+str(number_particles)+'_L'+str(number_sites)
    if isfile(basis_name+'.npy') == False:
        number_states = int(round(binom(number_sites, number_particles)))
        full_base = asarray(list(itertools.product(range(2), repeat=number_sites)), dtype=int)
        base = zeros((number_states,number_sites))
        midx = 0
        for m in range(0, full_base.shape[0]):
            if sum(full_base[m,:])==number_particles:
                print m, midx
                base[midx,:]=full_base[m,:]
                midx = midx+1
        save(basis_name, base)
        print 'saved new basis:', basis_name
    elif isfile(basis_name+'.npy') == True:
        print 'loading basis:', basis_name
        base = load(basis_name+'.npy')
    
    
    number_states = base.shape[0]


    nn_ham_name = 'HamiltonianData/NearestNeighbor_HoppingHamiltonian_N'+str(number_particles)+'_L'+str(number_sites)                    
    if isfile(nn_ham_name+'.npy') == False:
        nn_ham = nn_hopping(base, number_particles, number_sites)
        save(nn_ham_name, nn_ham)
        print 'saved new nearest-neighbor hopping Hamiltonian: ', nn_ham_name
    elif isfile(nn_ham_name+'.npy') == True:
        print 'loading nearest-neighbor hopping Hamiltonian: ', nn_ham_name
        nn_ham = load(nn_ham_name+'.npy')
    '''Commenting out the NNN hopping'''       
#    nnn_ham_name = 'HamiltonianData/NextNearestNeighbor_HoppingHamiltonian_N'+str(number_particles)+'_L'+str(number_sites)                    
#    if isfile(nnn_ham_name+'.npy') == False:
#        nnn_ham = nnn_hopping(base, number_particles, number_sites)
#        save(nnn_ham_name, nnn_ham)
#        print 'saved new next-nearest-neighbor hopping Hamiltonian: ', nnn_ham_name
#    elif isfile(nnn_ham_name+'.npy') == True:
#        print 'loading next-nearest-neighbor hopping Hamiltonian: ', nnn_ham_name
#        nnn_ham = load(nnn_ham_name+'.npy')
        
    ham = J*nn_ham#+0*(J/sqrt(2*(1-cos(gamma))))*nnn_ham

    if (ham.todense() == (ham.todense()).transpose()).all() != True:
        print (ham.todense() == (ham.todense()).transpose()).all()
        print 'Hopping is non-symmetric!'
        return
        
    print 'generating interaction Hamiltonian'
    if vNN == 0:
        print 'free fermions, no interactions'
    else:
        for idx_state in range(0,number_states):
            hidx, v = int_ham(base, number_particles, number_sites, gamma, vNN, vNNN, angle_fraction, idx_state, independent_v1_v2)
            ham[hidx, hidx] = v
        
            
    if spectrum == 0:
        print 'calculating eigenvectors and spectrum'
        eig_val, eig_vec = la.eigsh(ham, k=noEV, which='SA', maxiter=10000, return_eigenvectors=True)
    
        index_ev = eig_val.ravel().argsort()
        eigval = eig_val[index_ev]
        eigvec = eig_vec[:,index_ev]

        return eigval, eigvec, base
        
    elif spectrum == 1:

        print 'calculating spectrum'
        eig_val = la.eigsh(ham, k=noEV, which='SA', maxiter = 10000, return_eigenvectors=False)
        
        index_ev = eig_val.ravel().argsort()
        eigval = eig_val[index_ev]
        
        return angle_fraction, eigval, v1, v2, v3, J, vNN




