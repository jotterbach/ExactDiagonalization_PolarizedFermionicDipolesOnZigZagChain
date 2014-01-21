# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:09:05 2013

@author: jotterbach

"""


from numpy import *
from ED_HalfFilling import EigSys_HalfFilling
from DotProduct import scalar_prod
from multiprocessing import *
from multiprocessing import Pool
import matplotlib.pyplot as plt
from ParallelizationTools import info
from os.path import *
from scipy.special import *
from scipy.linalg import qr
from DotProduct import scalar_prod
from Correlation_Generator import *
from datetime import datetime


''' define the datestamp for the filenames '''
date = str(datetime.now())
now = date[0:10]+'_'+date[11:13]+'h'+date[14:16]+'m'

def AngleSpectrum(number_particles, noEV, gamma, hopping, angle):
    """
    AngleSpectrum(number_particles, noEV, gamma, hopping, angle):
    
    computes the energy eigenspectrum as a function of the angle of the dipoles
    with the chain axis given an unit interaction V and a hopping J

    parameters of the function:
        number_particles: number of particles in the problem
        noEV: number of eigenvalues being calculated
        gamma: opening angle of the zig-zag chain
        hopping: hopping parameter in units of interaction V
        angle: array containing the angles as a multiple of **PI**
    """

    ''' default values for other methods that are being called by the current
        function '''
    spectrum = 1 #ensures that the spectrum is calculated in EigSys_HalfFilling
    independet_v1_v2 = 1 #makes v1 and v2 independent of each other
    number_sites = 2*number_particles #condition for half-filling
    interaction_strength = 1 #unit of energy
#    number_particles = 6
#    noEV = 5*number_sites #degeneracy of GS requires noEV>number_sites
#    hopping = .1
#    gamma = 2*pi/3
#    angle = linspace(-.8,-.7,41)    

    ''' intialization of variables that will be stored for later use '''
    eigval = zeros((angle.shape[0],noEV), dtype = float)
    degeneracies = zeros((angle.shape[0],1))
    v1 = zeros((angle.shape[0],1))
    v2 = zeros((angle.shape[0],1))
    v3 = zeros((angle.shape[0],1))
    
    
    ''' actual method call '''    
    if __name__ == 'DiagonalizationMethods':
        info('main line')
        pool = Pool()
        ''' invocation of the eigenvalue procedure '''
        it = [pool.apply_async(EigSys_HalfFilling, (number_particles, number_sites, hopping, interaction_strength, angle[angle_idx], noEV, spectrum, gamma, independet_v1_v2)) for angle_idx in range(0,angle.shape[0])]
        for ridx in it:
            angle_idx = nonzero(angle == ridx.get()[0])        
            eigval[angle_idx,:]= ridx.get()[1]#floor(10*around(real(ridx.get()[1]),decimals = 2))/10
            degeneracies[angle_idx] = sum((eigval[angle_idx,:] == eigval[angle_idx,0]).astype(int))
            v1[angle_idx]=ridx.get()[2]
            v2[angle_idx]=ridx.get()[3]
            v3[angle_idx]=ridx.get()[4]
            print 'angle:', angle[angle_idx], '\nground-state degeneracy:', degeneracies[angle_idx]
                
        filename = 'FigureData/'+now+'_AngleSpectrum_N'+str(number_particles)+'_J'+str(hopping).replace('.','-')+'_vdd'+str(interaction_strength).replace('.','-')
        save(filename+'_EigVals', eigval)
        save(filename+'_angle', angle)
        print 'saved: '+filename

def InteractionSpectrum(number_particles, noEV, gamma, angle, interaction_strength):
    ''' computes the eigenvalue spectrum for a given angle 
        as a function of the interaction strength in units of J
        
        parameters of the function:
            number_particles: number of particles in the problem
            noEV: number of eigenvalues being calculated
            gamma: opening angle of the zig-zag chain
            angle: array containing the angles as a multiple of **PI**
            interaction_strength: interaction in units of hopping J
    '''
        
    ''' default values for other methods that are being called by the current
        function '''
    spectrum = 1 #ensures that the spectrum is calculated in EigSys_HalfFilling
    independent_v1_v2 = 1 #makes v1 and v2 independent of each other
    number_sites = 2*number_particles #condition for half-filling
    hopping = 1 #unit of energy


    ''' intialization of variables that will be stored for later use '''
    eigval = zeros((len(interaction_strength),noEV), dtype = float)
    v1 = zeros((interaction_strength.shape[0],1))
    v2 = zeros((interaction_strength.shape[0],1))
    v3 = zeros((interaction_strength.shape[0],1))

    ''' actual method call '''
    if __name__ == 'DiagonalizationMethods':
        info('main line')
        pool = Pool()
        ''' invocation of eigenvalue procedure '''
        it = [pool.apply_async(EigSys_HalfFilling, (number_particles, number_sites, hopping, interaction_strength[idx], angle, noEV, spectrum, gamma, independent_v1_v2)) for idx in range(len(interaction_strength))]
        for ridx in it:
            idx = nonzero(interaction_strength == ridx.get()[6])
            v1=ridx.get()[2]
            v2=ridx.get()[3]
            v3=ridx.get()[4]
            eigval[idx,:]= ridx.get()[1]#floor(10*around(real(ridx.get()[1]),decimals = 2))/10
            print 'interaction:', interaction_strength[idx], 'interaction constants: ', v1,v2,v3
       
    filename = 'FigureData/'+now+'_InteractionSpectrum_N'+str(number_particles)+'_J'+str(hopping).replace('.','-')+'_Theta'+str(angle).replace('.','-')
    save(filename+'_EigVals', eigval)
    save(filename+'_interaction',interaction_strength)
    print 'saved: '+filename


def HoppingSpectrum(number_particles, noEV, gamma, angle, hopping):
    ''' computes the eigenvalue spectrum for given interactions as a function 
        of the hopping in units of interaction V
    
        parameters of the function:
            number_particles: number of particles in the problem
            noEV: number of eigenvalues being calculated
            gamma: opening angle of the zig-zag chain
            angle: array containing the angles as a multiple of **PI**
            hopping: hopping in units of interaction V
    '''
    
    ''' default values for other methods that are being called by the current
        function '''
    spectrum = 1 #ensures that the spectrum is calculated in EigSys_HalfFilling
    independent_v1_v2 = 1 #makes v1 and v2 independent of each other
    number_sites = 2*number_particles #condition for half-filling
    interaction_strength = 1 #unit of energy
        
    ''' intialization of variables that will be stored for later use '''
    eigval = zeros((len(hopping),noEV), dtype = float)
    v1 = zeros((hopping.shape[0],1))
    v2 = zeros((hopping.shape[0],1))
    v3 = zeros((hopping.shape[0],1))    
    
    ''' actual method call '''
    if __name__ == 'DiagonalizationMethods':
        info('main line')
        pool = Pool()
        ''' invocation of eigenvalue procedure '''        
        it = [pool.apply_async(EigSys_HalfFilling, (number_particles, number_sites, hopping[idx], interaction_strength, angle, noEV, spectrum, gamma, independent_v1_v2)) for idx in range(len(hopping))]
        for ridx in it:
            idx = nonzero(hopping == ridx.get()[5])
            v1=ridx.get()[2]
            v2=ridx.get()[3]
            v3=ridx.get()[4]
            eigval[idx,:]= ridx.get()[1]
            print 'hopping:', hopping[idx], 'interactions: ', v1,v2,v3
       

    filename = 'FigureData/'+now+'_HoppingSpectrum-nnhopping_N'+str(number_particles)+'_vdd'+str(interaction_strength).replace('.','-')+'_Theta'+str(angle).replace('.','-')
    save(filename+'_EigVals', eigval)
    save(filename+'_hopping', hopping)
    print 'saved: '+filename

def DensityCorrelations(number_particles, noEV, gamma, angle, hopping, degeneracy):
    ''' computes the density correlation function for a given set of angle,
    interaction and hopping'''
    
    ''' default values for other methods that are being called by the current
        function '''
    spectrum = 0 #ensures that the spectrum AND the eigenvectors are calculated in EigSys_HalfFilling
    independent_v1_v2 = 1 #makes v1 and v2 independent of each other
    number_sites = 2*number_particles #condition for half-filling
    interaction_strength = 1 #unit of energy
    
    ''' function specific parameter initilaization '''
    eigval, eigvec, basisstates = EigSys_HalfFilling(number_particles, number_sites, hopping, interaction_strength, angle, noEV, spectrum, gamma, independent_v1_v2)
    eigval = around(real(eigval),decimals = 2)

    print '\nlow-energy spectrum: \n', eigval
    print 'GS degeneracy:', degeneracy
    eigvec = eigvec.astype(complex)

    if degeneracy > 1:
        print '\nOrthogonalizing GS manifold'
        eigvec_GS = zeros((eigvec.shape[0],degeneracy), dtype = complex)
        for m in range(degeneracy):
            eigvec_GS[:,m] = eigvec[:,m]
        Q, R = qr(eigvec_GS, mode = 'economic')    
        for m in range(degeneracy):
            eigvec[:,m] = Q[:,m] 
        del Q, R, eigvec_GS


    number_states = basisstates.shape[0]
        
    if __name__ == 'DiagonalizationMethods':
        ''' local density '''
        print '\nCalculating local density'
        local_density = zeros((2*number_particles,1), dtype = float)
        pool = Pool()
        for deg_idx in range(0,degeneracy):
            print 'state index: ', deg_idx
            it = [pool.apply_async(loc_den, (basisstates, number_particles, number_states, eigvec[:,deg_idx], site_idx)) for site_idx in range(0,2*number_particles)]
            for ridx in it:
                site_idx = ridx.get()[0]
                local_density[site_idx] += real(ridx.get()[1])/degeneracy
        
        ''' density-density correlation '''
        print '\nCalculating density-density correlations'
        g2 = zeros((number_sites,1), dtype = float)
        for deg_idx in range(0,degeneracy):
            print 'state index: ', deg_idx
            it = [pool.apply_async(pair_corr, (basisstates, number_particles, number_sites, number_states, eigvec[:,deg_idx], site_idx)) for site_idx in range(0,number_sites)]
            for ridx in it:
                site_idx = ridx.get()[0]
                g2[site_idx] += real(ridx.get()[1])/degeneracy
        
        filename='FigureData/'+now+'_Correlations_N'+str(number_particles)+'_J'+str(hopping).replace('.','-')+'_vdd'+str(interaction_strength).replace('.','-')+'_Theta'+str(angle).replace('.','-')
        save(filename+'_local_density', local_density)
        save(filename+'_g2', g2)
        
        print 'saved: '+filename
