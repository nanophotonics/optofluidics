#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 10:25:44 2016

@author: Andrei
"""
from __future__ import division
import numpy as np
import slmpy
import time

#custom exceptions used in this class
class WavelengthError(Exception):
    pass

#The purpose of this class is to control an SLM: useful functions are class initialisation, useSLM and closeslm.

class SLM():
    def __init__(self,slmp, wav, sep, mtr=1, width=10):
        """
        Sets up instance of class: slmp is pixel size of SLM, usually [512 512], wav is thecentral wavelength used,
        sep is the spatial frequency used to separate 1st order, mtr is the monitor on which the SLM can be found, usually 1.
        width is the width (in nm) of incident wavelength (default: 10 nm)
        """
        self.slmpix=slmp
        self.spfrq=sep
        self.wavelength=wav
        self.width=width
        self.inverter=np.loadtxt('asinc.txt')
        self.slm = slmpy.SLMdisplay(monitor = mtr)
        resX, resY = self.slm.getSize()
        self.slmstate=np.zeros((resY,resX,3))
        self.calibData = np.load('SLMCalibrationParam.npy')
        #upper edges of calibration bands
        #this must be ipdated if calibration is updated!!!!!!!
        self.calibBandWidth=10
        self.calibBands = np.arange(500,841,self.calibBandWidth) # first band is 480-490 nm
        self.minAllowedWav = 490
        self.maxAllowedWav = 840
        return


    def separate(self):
        x=np.linspace(-(self.slmpix[0]-1)/2,(self.slmpix[0]-1)/2,self.slmpix[0])
        xmat=np.tile(x,(self.slmpix[1],1))
        y=np.linspace(-(self.slmpix[1]-1)/2,(self.slmpix[1]-1)/2,self.slmpix[1])
        ymat=np.tile(y,(self.slmpix[0],1))
        ymat=np.transpose(ymat)  
        #adjust spfrq based on length
        spfrq = self.spfrq / (self.wavelength/655)
        return (2*np.pi*spfrq/np.sqrt(2)*(xmat-ymat))
    
    def genhol(self, distribution, separator):
        ampaim=np.absolute(distribution)
        phaseaim=np.angle(distribution)
        ampcorr=1-np.interp(ampaim,self.inverter[:,0],self.inverter[:,1])
        tosend=np.multiply(ampcorr,np.mod(separator+phaseaim,np.pi*2))
        return np.mod(tosend,2*np.pi)


    def SLMprocess_old(self, toset):
        """Relates target phase-delay to SLM values"""
        if self.wavelength <500:
            raise WavelengthError('Wavelength outside range!\n'+
            'Got {} but must be in range 500-800 nm.'.format(toset))
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
            raise WavelengthError('Wavelength outside range!\n'+
            'Got {} nm but must be in range 480-840 nm.'.format(toset))
        toset=toset.astype('uint16')
        print(np.max(toset))
        toset1= np.bitwise_and(toset,0xff)
        toset2=np.bitwise_and(np.right_shift(toset,8),0xff)
        blank=np.zeros(np.shape(toset1))
        tosend=np.concatenate((np.expand_dims(toset1, axis=2),np.expand_dims(toset2, axis=2),
                               np.expand_dims(blank, axis=2)), axis=2)
        return tosend
    
    def SLMprocess(self, toset):
        """
        Relates target phase-delay to SLM values
        
        Calibration data was obtained for 10 nm wide bands
        Other bands are approximated using weighted averages of calibration data
        """
        minWav = self.wavelength - self.width/2
        #print(minWav)
        maxWav = self.wavelength + self.width/2
        
        if minWav <self.minAllowedWav or maxWav > self.maxAllowedWav:
            raise WavelengthError('Wavelength outside range!\n'+
            'Got {}-{} nm but must be in range 490-840 nm.'.format(
            minWav, maxWav))
        #print(minWav)
        #print(maxWav)
        #break up in to relevant calibration bands
        minInd, maxInd = np.searchsorted(self.calibBands, [minWav, maxWav])
        #print(minInd)
        #print(maxInd)
        wholeBands = self.calibBands[minInd+1:maxInd]
        #print(wholeBands)
        minWeight = (self.calibBands[minInd]-minWav) /10
        maxWeight = (maxWav + self.calibBandWidth 
                    - self.calibBands[maxInd] ) /10
        #print(minWeight)
        #print(maxWeight)
        #prefent error in self._applyCalib
        if minWav == 490:
            minWav += 1
        #if maxWav == 840:
        #    maxWav -= 1
        temp = self._applyCalib(toset,minWav)*minWeight
        #print(temp)
        temp += self._applyCalib(toset,maxWav)*maxWeight
        #print(temp)
        for band in wholeBands:
            temp += self._applyCalib(toset,band)
        #print(temp)
        temp = temp / (len(wholeBands)+minWeight+maxWeight)
        #print(temp)
        temp=temp.astype('uint16')
        #print(np.max(temp))
        toset1= np.bitwise_and(temp,0xff)
        toset2=np.bitwise_and(np.right_shift(temp,8),0xff)
        blank=np.zeros(np.shape(toset1))
        tosend=np.concatenate((np.expand_dims(toset1, axis=2),np.expand_dims(toset2, axis=2),
                               np.expand_dims(blank, axis=2)),
                                axis=2)
        return tosend
    
    def _applyCalib(self, toset, band):
        """
        Converts the toset value according to the specified band to an SLM value.
        
        band is expected as a wavelength /nm
        band should be within the chosen calibration band the lower band bound
        is excluded and the upper band bound included.
        """
        if band <=self.minAllowedWav or band > self.maxAllowedWav:
            raise WavelengthError('Wavelength outside range!\n'+
            'Got {} nm but must be in range 490-840 nm.'.format(band))
        
        bandInd = np.searchsorted(self.calibBands, band)
        #print(bandInd)
        #print(bandInd)
        toset=( (toset/np.pi)**5 * (self.calibData[bandInd,0]) 
            + (toset/np.pi)**4 * (self.calibData[bandInd,1]) 
            + (toset/np.pi)**3 * (self.calibData[bandInd,2])
            + (toset/np.pi)**2 * (self.calibData[bandInd,3]) 
            + (toset/np.pi) * (self.calibData[bandInd,4]) )*65535*0.01
        return toset
    
    def SLMsend(self, toset):
        self.slmstate[0:self.slmpix[0],0:self.slmpix[1],:]=toset
        self.slm.updateArray(self.slmstate.astype('uint8'))
        return

    def useSLM(self, distribution):
        """Disply distribution in first order reflection of SLM"""
        separator=self.separate()
        holo=self.genhol(distribution,separator)
        #print(np.max(holo))
        tosend=self.SLMprocess(holo)
        self.SLMsend(tosend)
        #time.sleep(0.5)
        return tosend
    
    def closeslm(self):
        self.slm.close()
        return