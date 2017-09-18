#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 10:25:44 2016

@author: Andrei
"""
import numpy as np
import slmpy
import time
inverter=np.loadtxt('asinc.txt')



def separate(slmpix,spfrq):
    x=np.linspace(-(slmpix[0]-1)/2,(slmpix[0]-1)/2,slmpix[0])
    xmat=np.tile(x,(slmpix[1],1))
    y=np.linspace(-(slmpix[1]-1)/2,(slmpix[1]-1)/2,slmpix[1])
    ymat=np.tile(y,(slmpix[0],1))
    ymat=np.transpose(ymat)

    return (2*np.pi*spfrq/np.sqrt(2)*(xmat-ymat))
    
def genhol(distribution,separator,inverter):
    intaim=np.absolute(distribution)
    phaseaim=np.angle(distribution)
    intcorr=1-np.interp(intaim,inverter[:,0],inverter[:,1])
    tosend=np.multiply(intcorr,np.mod(separator+phaseaim,np.pi*2))
    return np.mod(tosend,2*np.pi)


def SLMprocess_old(toset,wavelength):
    #tosend is the desired phase!
    if wavelength <500:
        raise Exception('Wavelength outwith range')
    elif wavelength <= 550:
        #old versions are commented (I think incorrect!)
        #toset=(0.15+toset/(np.pi*22.3))*65535
        toset=(0.168+toset/(np.pi*22.3))*65535
    elif wavelength <= 600:
        #toset=(0.15+toset/(np.pi*19.9))*65535
        toset=(0.161+toset/(np.pi*19.9))*65535
    elif wavelength <= 650:
        #toset=(0.15+toset/(np.pi*17.5))*65535
        toset=(0.159+toset/(np.pi*17.5))*65535
    elif wavelength <= 700:
        #toset=(0.13+toset/(np.pi*15.5))*65535
        toset=(0.148+toset/(np.pi*15.5))*65535
    elif wavelength <= 750:
        #toset=(0.13+toset/(np.pi*14.6))*65535
        toset=(0.143+toset/(np.pi*14.6))*65535
    elif wavelength <= 800:
        #toset=(0.13+toset/(np.pi*13.3))*65535
        toset=(0.133+toset/(np.pi*13.3))*65535
    elif wavelength <= 810:
        toset=( (toset/np.pi)**5 * (-0.38455961) 
            + (toset/np.pi)**4 * (2.09392385) 
            + (toset/np.pi)**3 * (-1.24688497)
            + (toset/np.pi)**2 * (-10.03337945) 
            + (toset/np.pi) * (26.62612235))*65535*0.01
        #params: -0.38455961, 2.09392385, -1.24688497, -10.03337945, 26.62612235
    elif wavelength <= 820:
        toset=( (toset/np.pi)**5 * (-0.25994767) 
            + (toset/np.pi)**4 * (0.68186957) 
            + (toset/np.pi)**3 * (3.93830709)
            + (toset/np.pi)**2 * (-17.43210166) 
            + (toset/np.pi) * (30.11298442))*65535*0.01
        #params: -0.25994767, 0.68186957, 3.93830709, -17.43210166, 30.11298442
    elif wavelength <= 830:
        toset=( (toset/np.pi)**5 * (-0.53073512) 
            + (toset/np.pi)**4 * (1.96923007)
            + (toset/np.pi)**3 * (2.73958557)
            + (toset/np.pi)**2 * (-19.02779726) 
            + (toset/np.pi) * (32.28894872))*65535*0.01
        #params: -0.53073512, 1.96923007, 2.73958557, -19.02779726, 32.28894872
    elif wavelength <= 840:
        toset=( (toset/np.pi)**5 * (-0.61121743) 
            + (toset/np.pi)**4 * (2.76524192) 
            + (toset/np.pi)**3 * (0.35339914)
            + (toset/np.pi)**2 * (-16.59084706) 
            + (toset/np.pi) * (31.90949444))*65535*0.01
        #params: -0.61121743, 2.76524192, 0.35339914, -16.59084706, 31.90949444
    else:
        raise Exception('Wavelength outwith range')
    toset=toset.astype('uint16')
    toset1= np.bitwise_and(toset,0xff)
    toset2=np.bitwise_and(np.right_shift(toset,8),0xff)
    blank=np.zeros(np.shape(toset1))
    tosend=np.concatenate((np.expand_dims(toset1, axis=2),np.expand_dims(toset2, axis=2),np.expand_dims(blank, axis=2)), axis=2)
    return tosend
    
def SLMprocess(toset,wavelength, width):
    #tosend is the desired phase!
    if wavelength <480:
        raise ValueError('Wavelength outwith range')
    elif wavelength <= 550:
        #old versions are commented (I think incorrect!)
        #toset=(0.15+toset/(np.pi*22.3))*65535
        toset=(0.168+toset/(np.pi*22.3))*65535
    elif wavelength <= 600:
        #toset=(0.15+toset/(np.pi*19.9))*65535
        toset=(0.161+toset/(np.pi*19.9))*65535
    elif wavelength <= 650:
        #toset=(0.15+toset/(np.pi*17.5))*65535
        toset=(0.159+toset/(np.pi*17.5))*65535
    elif wavelength <= 700:
        #toset=(0.13+toset/(np.pi*15.5))*65535
        toset=(0.148+toset/(np.pi*15.5))*65535
    elif wavelength <= 750:
        #toset=(0.13+toset/(np.pi*14.6))*65535
        toset=(0.143+toset/(np.pi*14.6))*65535
    elif wavelength <= 800:
        #toset=(0.13+toset/(np.pi*13.3))*65535
        toset=(0.133+toset/(np.pi*13.3))*65535
    elif wavelength <= 810:
        toset=( (toset/np.pi)**5 * (-0.38455961) 
            + (toset/np.pi)**4 * (2.09392385) 
            + (toset/np.pi)**3 * (-1.24688497)
            + (toset/np.pi)**2 * (-10.03337945) 
            + (toset/np.pi) * (26.62612235))*65535*0.01
        #params: -0.38455961, 2.09392385, -1.24688497, -10.03337945, 26.62612235
    elif wavelength <= 820:
        toset=( (toset/np.pi)**5 * (-0.25994767) 
            + (toset/np.pi)**4 * (0.68186957) 
            + (toset/np.pi)**3 * (3.93830709)
            + (toset/np.pi)**2 * (-17.43210166) 
            + (toset/np.pi) * (30.11298442))*65535*0.01
        #params: -0.25994767, 0.68186957, 3.93830709, -17.43210166, 30.11298442
    elif wavelength <= 830:
        toset=( (toset/np.pi)**5 * (-0.53073512) 
            + (toset/np.pi)**4 * (1.96923007)
            + (toset/np.pi)**3 * (2.73958557)
            + (toset/np.pi)**2 * (-19.02779726) 
            + (toset/np.pi) * (32.28894872))*65535*0.01
        #params: -0.53073512, 1.96923007, 2.73958557, -19.02779726, 32.28894872
    elif wavelength <= 840:
        toset=( (toset/np.pi)**5 * (-0.61121743) 
            + (toset/np.pi)**4 * (2.76524192) 
            + (toset/np.pi)**3 * (0.35339914)
            + (toset/np.pi)**2 * (-16.59084706) 
            + (toset/np.pi) * (31.90949444))*65535*0.01
        #params: -0.61121743, 2.76524192, 0.35339914, -16.59084706, 31.90949444
    else:
        raise ValueError('Wavelength outwith range')
    toset=toset.astype('uint16')
    toset1= np.bitwise_and(toset,0xff)
    toset2=np.bitwise_and(np.right_shift(toset,8),0xff)
    blank=np.zeros(np.shape(toset1))
    tosend=np.concatenate((np.expand_dims(toset1, axis=2),np.expand_dims(toset2, axis=2),np.expand_dims(blank, axis=2)), axis=2)
    return tosend

def SLMinit(mtr):
    slm = slmpy.SLMdisplay(monitor = mtr)
    resX, resY = slm.getSize()
    slmstate=np.zeros((resY,resX,3))
    return (slm,slmstate)
    
def SLMsend(toset,slmstate,slm):
    slmstate[0:512,0:512,:]=toset
    slm.updateArray(slmstate.astype('uint8'))
    return

def useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix):
    separator=separate(slmpix,spfrq)
    holo=genhol(distribution,separator,inverter)
    tosend=SLMprocess(holo,wavelength)
    SLMsend(tosend,slmstate,slm)
    #time.sleep(0.5)
    return tosend