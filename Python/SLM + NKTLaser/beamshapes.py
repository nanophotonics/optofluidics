#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 10:18:20 2016

@author: Andrei
"""
import scipy.special
import numpy as np
import scipy.misc
from scipy.ndimage.interpolation import zoom, rotate
from numpy import genfromtxt


def LG(l, p, w, slmpix):
    def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return(rho, phi)
    nx=slmpix[0]
    ny=slmpix[1]
    x=np.linspace(-(nx-1)/2,(nx-1)/2,nx)
    y=np.linspace(-(ny-1)/2,(ny-1)/2,ny)
    xm, ym = np.meshgrid(x, y)
    rho, phi = cart2pol(xm,ym)
    LG= \
        np.exp(-np.square(rho/w))* \
        np.power(np.sqrt(2)*rho/w,l)* \
        np.polyval(scipy.special.genlaguerre(p, l),2*np.power(rho/w,2))* \
        np.exp(1.j*l*phi)
    #edit by Marius
    #notice gives intensities greater 1 to SLM therefore normalise by maximum real part
    LG=LG/np.max(np.abs(np.real(LG)))
    return LG

def test(rep,slmpix):
    def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return(rho, phi)
    nx=slmpix[0]
    ny=slmpix[1]
    x=np.linspace(-(nx-1)/2,(nx-1)/2,nx)
    y=np.linspace(-(ny-1)/2,(ny-1)/2,ny)
    xm, ym = np.meshgrid(x, y)
    rho, phi = cart2pol(xm,ym)
    dist=np.exp(rho/rep*2*np.pi*1.j)
    return dist
    
def gausscomp(sigma, slmpix):
    x=np.linspace(-(slmpix[0]-1)/2,(slmpix[0]-1)/2,slmpix[0])
    xmat=np.tile(x,(slmpix[1],1))
    y=np.linspace(-(slmpix[1]-1)/2,(slmpix[1]-1)/2,slmpix[1])
    ymat=np.tile(y,(slmpix[0],1))
    ymat=np.transpose(ymat)
    gauss=np.exp((-np.square(xmat)-np.square(ymat))/(2*sigma*sigma))   
    return -(1-gauss+gauss[511,511])  #note the -ve is s.t. you have pi phase shift
    
def square(L,slmpix):
    intaim=np.zeros(slmpix)
    intaim[(slmpix[0]-L)/2:(slmpix[0]+L)/2,(slmpix[1]-L)/2:(slmpix[1]+L)/2]=intaim[(slmpix[0]-L)/2:(slmpix[0]+L)/2,(slmpix[1]-L)/2:(slmpix[1]+L)/2]+1
    return intaim
'''    
def pic(zoomf,slmpix):
    phases = genfromtxt('phases1.csv', delimiter=',')
    amp = genfromtxt('amplitude1.csv', delimiter=',')
    
    phases= zoom(phases, zoomf)
    phases[phases>1.5]=np.pi
    phases[phases<1.5]=0
    L1=np.shape(phases)[0]
    L2=np.shape(phases)[1]
    phaseaim=np.zeros(slmpix)
    phaseaim[(slmpix[0]-L1)/2:(slmpix[0]+L1)/2,(slmpix[1]-L2)/2:(slmpix[1]+L2)/2]=phases
    
    amp= zoom(amp, zoomf)
    amp[amp>0.5]=1
    amp[amp<0.5]=0
    L1=np.shape(amp)[0]
    L2=np.shape(amp)[1]
    ampaim=np.zeros(slmpix)
    ampaim[(slmpix[0]-L1)/2:(slmpix[0]+L1)/2,(slmpix[1]-L2)/2:(slmpix[1]+L2)/2]=amp
    return ampaim*np.exp(1.j*phaseaim)
'''

def pic(path, zoomf,slmpix):
    amp = genfromtxt(path, delimiter=',')    
    amp= zoom(amp, zoomf)
    amp=amp/np.amax(amp)
    L1=np.shape(amp)[0]
    L2=np.shape(amp)[1]
    ampaim=np.zeros(slmpix)
    ampaim[(slmpix[0]-L1)/2:(slmpix[0]+L1)/2,(slmpix[1]-L2)/2:(slmpix[1]+L2)/2]=amp
    return ampaim

    
def nothing(slmpix):
    intaim=np.ones(slmpix)
    return intaim
    
def angled(distribution, phy, phx, slmpix):
    x=np.linspace(-(slmpix[0]-1)/2,(slmpix[0]-1)/2,slmpix[0])
    xmat=np.tile(x,(slmpix[1],1))
    y=np.linspace(-(slmpix[1]-1)/2,(slmpix[1]-1)/2,slmpix[1])
    ymat=np.tile(y,(slmpix[0],1))
    ymat=np.transpose(ymat)
    return distribution*np.exp(1.j*phx*xmat/slmpix[0])*np.exp(1.j*phy*ymat/slmpix[1])

def displace(distribution, xmove, ymove, slmpix):
    xmove=int(xmove)
    ymove=int(ymove)
    distribution = np.roll(distribution,xmove,axis=0)
    distribution = np.roll(distribution,ymove,axis=1)
    if xmove<0:
        distribution[xmove:slmpix[0],:]=0
    else:
        distribution[0:xmove,:]=0
    if ymove<0:
        distribution[:,ymove:slmpix[1]]=0
    else:
        distribution[:,0:ymove]=0
    return distribution
def rotator(distribution, angle):
    return rotate(np.real(distribution), angle, reshape=False)+1.j*rotate(np.imag(distribution), angle, reshape=False)
    
def corrector(distribution, sigma, slmpix):
    return np.multiply(distribution,gausscomp(sigma,slmpix))