#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 10:25:44 2016

@author: Andrei
"""
import numpy as np
import slmpy
import time

#The purpose of this class is to control an SLM: useful functions are class initialisation, useSLM and closeslm.

class SLM():
    def __init__(self,slmp, wav, sep, mtr):
        #Sets up instance of class: slmp is pixel size of SLM, usually [512 512], wav is the wavelength used,
        #sep is the spatial frequency used to separate 1st order, mtr is the monitor on which the SLM can be found, usually 1.
        self.slmpix=slmp
        self.spfrq=sep
        self.wavelength=wav
        self.inverter=np.loadtxt('asinc.txt')
        self.slm = slmpy.SLMdisplay(monitor = mtr)
        resX, resY = self.slm.getSize()
        self.slmstate=np.zeros((resY,resX,3))
        return


    def separate(self):
        x=np.linspace(-(self.slmpix[0]-1)/2,(self.slmpix[0]-1)/2,self.slmpix[0])
        xmat=np.tile(x,(self.slmpix[1],1))
        y=np.linspace(-(self.slmpix[1]-1)/2,(self.slmpix[1]-1)/2,self.slmpix[1])
        ymat=np.tile(y,(self.slmpix[0],1))
        ymat=np.transpose(ymat)  
        return (2*np.pi*self.spfrq/np.sqrt(2)*(xmat-ymat))
    
    def genhol(self, distribution, separator):
        ampaim=np.absolute(distribution)
        phaseaim=np.angle(distribution)
        ampcorr=1-np.interp(ampaim,self.inverter[:,0],self.inverter[:,1])
        tosend=np.multiply(ampcorr,np.mod(separator+phaseaim,np.pi*2))
        return np.mod(tosend,2*np.pi)


    def SLMprocess(self, toset):
        if self.wavelength <500:
            raise Exception('Wavelength outwith range')
        elif self.wavelength <= 550:
            toset=(0.15+toset/(np.pi*22.3))*65535
        elif self.wavelength <= 600:
            toset=(0.15+toset/(np.pi*19.9))*65535
        elif self.wavelength <= 650:
            toset=(0.15+toset/(np.pi*17.5))*65535
        elif self.wavelength <= 700:
            toset=(0.13+toset/(np.pi*15.5))*65535
        elif self.wavelength <= 750:
            toset=(0.13+toset/(np.pi*14.6))*65535
        elif self.wavelength <= 800:
            toset=(0.13+toset/(np.pi*13.3))*65535
        else:
            raise Exception('Wavelength outwith range')
        toset=toset.astype('uint16')
        toset1= np.bitwise_and(toset,0xff)
        toset2=np.bitwise_and(np.right_shift(toset,8),0xff)
        blank=np.zeros(np.shape(toset1))
        tosend=np.concatenate((np.expand_dims(toset1, axis=2),np.expand_dims(toset2, axis=2),np.expand_dims(blank, axis=2)), axis=2)
        return tosend
    
    def SLMsend(self, toset):
        self.slmstate[0:self.slmpix[0],0:self.slmpix[1],:]=toset
        self.slm.updateArray(self.slmstate.astype('uint8'))
        return

    def useSLM(self, distribution):
        separator=self.separate()
        holo=self.genhol(distribution,separator)
        tosend=self.SLMprocess(holo)
        self.SLMsend(tosend)
        #time.sleep(0.5)
        return tosend
    
    def closeslm(self):
        self.slm.close()
        return