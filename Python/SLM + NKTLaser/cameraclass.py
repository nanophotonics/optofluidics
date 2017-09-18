# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 13:55:58 2016

@author: Hera
"""

import sys
import time
#from PySide.QtGui import *
from instrumental import instrument,gui, errors
import scipy.misc
import numpy as np
from autofunct import cropcurrent

class camera():
    def __init__(self,expos):
        self.maxexp=19.0
        self.pics=int(1)
        self.exp=float(0)
        self.cam = instrument('uc480') #uc480
        self.exposure(float(expos))
        self.exposure(float(expos))
        self.counter=0
        self.countmax=400
        print(self.cam._get_exposure())
        return

    def closecam(self):
        self.cam.stop_live_video()
        self.cam.close()
        return

    '''
    def takeimagecolour(self):
        b=np.zeros([1024,1280,3], dtype=int)
        
        for i in range(0, self.pics):
            frame_ready = self.cam.wait_for_frame()
            if frame_ready:
                a = self.cam.latest_frame()
                a=a[:,:,0:3]
            else:
                a=0
                raise NameError('Failed to get frame.')
            b=b+a
        b[b>255]=255
        #b=b.astype('float')
        #b=255*np.power((b/255),1/1.6)
        #b=b.astype('int')
        return b
        '''
        
    def takeimagebw(self):
        b=np.zeros([1024,1280], dtype=int)
        
        for i in range(0, self.pics):
            self.counter=self.counter+1
            frame_ready = False
            while (frame_ready==False):
                frame_ready = self.cam.wait_for_frame()
                if frame_ready:
                    a = self.cam.latest_frame()
                    #a=a[:,:,0:3]
                    #a=np.sum(a,2)
                else:
                    #a=0
                    #raise NameError('Failed to get frame.')
                    self.restartcam()
                #if (self.counter>=self.countmax):
                 #   self.restartcam()
            b=b+a
        b[b>255]=255
        return b#/3
        
        
    def restartcam(self):
        self.counter=0
        self.cam.stop_live_video()
        self.cam.close()
        self.cam = instrument('uc480')
        self.exposure(float(self.pics)*self.exp)
        print("Reset Camera")
        return
        
    '''
    def autoexposure(self):
        emax=60.0
        emin=0.0
        a=self.takeimagecolour()
        if np.sum(a==255)>10:
            emax=self.exp*self.pics
        else:
            emin=self.exp*self.pics
        while ((emax-emin)/(emax/2+emin/2)>=0.2):
            self.exposure((emax+emin)/2)
            a=self.takeimagecolour()
            if np.sum(a==255)>10:
                emax=(emax+emin)/2
            else:
                emin=(emax+emin)/2
        self.exposure(emin)
        return emin
        
    def autoexposurecut(self, cent, dist, angle):
        emax=60.0
        emin=0.0
        a=self.takeimagecolour()
        a=cropcurrent(a, cent, dist, angle)
        if np.sum(a==255)>10:
            emax=self.exp*self.pics
        else:
            emin=self.exp*self.pics
        while ((emax-emin)/(emax/2+emin/2)>=0.2):
            self.exposure((emax+emin)/2)
            a=self.takeimagecolour()
            a=cropcurrent(a, cent, dist, angle)
            if np.sum(a==255)>10:
                emax=(emax+emin)/2
            else:
                emin=(emax+emin)/2
        self.exposure(emin)
        return emin
    '''

    def autoexposure(self):
        emax=100.0
        emin=0.0
        a=self.takeimagebw()
        if np.sum(a==255)>10:
            emax=self.exp*self.pics
        else:
            emin=self.exp*self.pics
        while ((emax-emin)/(emax/2+emin/2)>=0.2):
            self.exposure((emax+emin)/2)
            a=self.takeimagebw()
            if np.sum(a==255)>10:
                emax=(emax+emin)/2
            else:
                emin=(emax+emin)/2
        self.exposure(emin)
        return emin
        
    def autoexposurecut(self, cent, dist, angle):
        emax=300.0
        emin=0.0
        a=self.takeimagebw()
        a=cropcurrent(a, cent, dist, angle)
        if np.sum(a==255)>10:
            emax=self.exp*self.pics
        else:
            emin=self.exp*self.pics
        while ((emax-emin)/(emax/2+emin/2)>=0.2):
            self.exposure((emax+emin)/2)
            a=self.takeimagebw()
            a=cropcurrent(a, cent, dist, angle)
            if np.sum(a==255)>10:
                emax=(emax+emin)/2
            else:
                emin=(emax+emin)/2
        self.exposure(emin)
        return emin
    
    def exposure(self, expo):
        self.pics=int(np.ceil(expo/self.maxexp))
        self.exp=expo/self.pics
        self.cam.start_live_video(exposure_time="{} millisecond".format(str(self.exp)))
        return
#print(np.shape(a))
#print(np.max(np.max(a)))