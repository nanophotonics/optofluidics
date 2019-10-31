# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:55:58 2019

@author: Takashi
"""

import sys
import time
#from PySide.QtGui import *
from instrumental import instrument,gui, errors
import scipy.misc
import numpy as np
from autofunct import cropcurrent
from autofunct import areainterest
from scipy.ndimage.interpolation import zoom, rotate
import matplotlib.image as mpimg
import scipy.ndimage.filters
import time
from datetime import datetime

class camera():
    def __init__(self,expos):
        """Initialises camera and prints current exposure"""
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
        """Close camera connetction
        
        Only one connetction may be open at a time!        
        """
        self.cam.stop_live_video()
        self.cam.close()
        return

    '''
    def takeimagecolour(self):
        b=np.empty([1024,1280,3], dtype=int)
        
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
        """Take a black & wight image"""
        b=np.empty([1024,1280], dtype=int)
        
        for i in range(0, self.pics):
            self.counter=self.counter+1
            frame_ready = False
            while (frame_ready==False):
                frame_ready = self.cam.wait_for_frame()
                if frame_ready:
                    a = self.cam.latest_frame()
                    #a=a[:,:,0:3]
                    #a=np.sum(a,2)
                    #a=a/3
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
        """Closes and restartys camera connection"""
        self.counter=0
        self.cam.stop_live_video()
        self.cam.close()
        self.cam = instrument('uc480')
        self.exposure(float(self.pics)*self.exp)
        print("Reset Camera")
        return

    def autoexposure(self):
        """Automatically adjust cexposure
        
        Returns exposure value"""
        emax=100.0
        emin=0.01
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
        """Autoexpose croped section of image?
        
        Returns exposure value"""
        emax=300.0
        emin=0.01
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
        """Sets camera exposure to expo
        
        Input: expo (scalar)"""
        self.pics=int(np.ceil(expo/self.maxexp))
        self.exp=expo/self.pics
        self.cam.start_live_video(exposure_time="{} millisecond".format(str(self.exp)))
        return

#setup data
cam=camera(10)
print(cam.autoexposure())
tint=5 #in min
number=6
impath='C://Users//Hera//Desktop//Camera Images//'

#camera loop
for i in range(0,number):
    #prints exposure value to console
    now = datetime.now()
    current_time = now.strftime("%H-%M-%S")
    print("Current Time =", current_time)
    print(cam.autoexposure())
    image=cam.takeimagebw()
    scipy.misc.imsave(impath+'Cam-{}.png'.format(current_time), image)
    time.sleep(tint*60)
cam.closecam()
