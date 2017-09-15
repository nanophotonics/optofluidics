#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 10:18:20 2016

@author: Andrei
"""
import scipy.special
import numpy as np

def LG(l, p, w, nx, ny):
    def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return(rho, phi)
    
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
