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

class beamshapes():
    def __init__(self,slmp):
        self.slmpix=slmp
        return


    def LG(self, l, p, w):
        def cart2pol(x, y):
            rho = np.sqrt(x**2 + y**2)
            phi = np.arctan2(y, x)
            return(rho, phi)
        nx=self.slmpix[0]
        ny=self.slmpix[1]
        x=np.linspace(-(nx-1)/2,(nx-1)/2,nx)
        y=np.linspace(-(ny-1)/2,(ny-1)/2,ny)
        xm, ym = np.meshgrid(x, y)
        rho, phi = cart2pol(xm,ym)
        LG= \
            np.exp(-np.square(rho/w))* \
            np.power(np.sqrt(2)*rho/w,l)* \
            np.polyval(scipy.special.genlaguerre(p, l),2*np.power(rho/w,2))* \
            np.exp(1.j*l*phi)
        return LG

    def test(self, rep):
        def cart2pol(x, y):
            rho = np.sqrt(x**2 + y**2)
            phi = np.arctan2(y, x)
            return(rho, phi)
        nx=self.slmpix[0]
        ny=self.slmpix[1]
        x=np.linspace(-(nx-1)/2,(nx-1)/2,nx)
        y=np.linspace(-(ny-1)/2,(ny-1)/2,ny)
        xm, ym = np.meshgrid(x, y)
        rho, phi = cart2pol(xm,ym)
        dist=np.exp(rho/rep*2*np.pi*1.j)
        return dist
    
    def gausscomp(self, sigma):
        x=np.linspace(-(self.slmpix[0]-1)/2,(self.slmpix[0]-1)/2,self.slmpix[0])
        xmat=np.tile(x,(self.slmpix[1],1))
        y=np.linspace(-(self.slmpix[1]-1)/2,(self.slmpix[1]-1)/2,self.slmpix[1])
        ymat=np.tile(y,(self.slmpix[0],1))
        ymat=np.transpose(ymat)
        gauss=np.exp((-np.square(xmat)-np.square(ymat))/(2*sigma*sigma))   
        return -(1-gauss+gauss[self.slmpix[0]-1,self.slmpix[1]-1])  #note the -ve is s.t. you have pi phase shift
    
    def square(self, L):
        intaim=np.zeros(self.slmpix)
        intaim[(self.slmpix[0]-L)/2:(self.slmpix[0]+L)/2,(self.slmpix[1]-L)/2:(self.slmpix[1]+L)/2]=intaim[(self.slmpix[0]-L)/2:(self.slmpix[0]+L)/2,(self.slmpix[1]-L)/2:(self.slmpix[1]+L)/2]+1
        return intaim

    def pic(self, path, zoomf):
        amp = genfromtxt(path, delimiter=',')    
        amp= zoom(amp, zoomf)
        amp=amp/np.amax(amp)
        L1=np.shape(amp)[0]
        L2=np.shape(amp)[1]
        ampaim=np.zeros(self.slmpix)
        ampaim[(self.slmpix[0]-L1)/2:(self.slmpix[0]+L1)/2,(self.slmpix[1]-L2)/2:(self.slmpix[1]+L2)/2]=amp
        return ampaim

    
    def nothing(self):
        intaim=np.ones(self.slmpix)
        return intaim
    
    def angled(self, distribution, phy, phx):
        x=np.linspace(-(self.slmpix[0]-1)/2,(self.slmpix[0]-1)/2,self.slmpix[0])
        xmat=np.tile(x,(self.slmpix[1],1))
        y=np.linspace(-(self.slmpix[1]-1)/2,(self.slmpix[1]-1)/2,self.slmpix[1])
        ymat=np.tile(y,(self.slmpix[0],1))
        ymat=np.transpose(ymat)
        return distribution*np.exp(1.j*phx*xmat/self.slmpix[0])*np.exp(1.j*phy*ymat/self.slmpix[1])
    
    def displace(self, distribution, xmove, ymove):
        xmove=int(xmove)
        ymove=int(ymove)
        distribution = np.roll(distribution,xmove,axis=0)
        distribution = np.roll(distribution,ymove,axis=1)
        if xmove<0:
            distribution[xmove:self.slmpix[0],:]=0
        else:
            distribution[0:xmove,:]=0
        if ymove<0:
            distribution[:,ymove:self.slmpix[1]]=0
        else:
            distribution[:,0:ymove]=0
        return distribution
    
    def rotator(self, distribution, angle):
        return rotate(np.real(distribution), angle, reshape=False)+1.j*rotate(np.imag(distribution), angle, reshape=False)
        
    def corrector(self, distribution, sigma):
        return np.multiply(distribution,self.gausscomp(sigma,self.slmpix))