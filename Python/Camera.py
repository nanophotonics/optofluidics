# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 12:11:58 2017
"""

#The purpose of this class is to control an camera.
import sys
import time
import scipy.misc
import numpy as np
from autofunct import cropcurrent
from instrumental import instrument, gui, errors
 
class camera():
    def __init__(self,expos):
              #Initialise camera with an exposure: this can be changed later, maxexp is set to 19ms and
              #is the maximum camera exposure for a single image.
        self.maxexp=19.0
        self.pics=int(0)
        self.exp=float(0)
        self.cam = instrument('uc480') #uc480
        self.exposure(float(expos))
        self.exposure(float(expos))
        self.counter=0
        self.countmax=400
        print(self.cam._get_exposure())
        return
 
    def closecam(self):
              #Close the camera before program ends.
        self.cam.stop_live_video()
        self.cam.close()
        return
 
    def takeimagebw(self):
              #Take a black and white image: note that if the exposure is greater than maxexp then
              #multiple images will be taken and combined.
        b=np.empty([1024,1280], dtype=int)
       
        for i in range(0, self.pics):
            self.counter=self.counter+1
            frame_ready = False
            while (frame_ready==False):
                frame_ready = self.cam.wait_for_frame()
                if frame_ready:
                    a = self.cam.latest_frame()
                else:
                    self.restartcam()
            b=b+a
        b[b>255]=255
        return b
       
        
    def restartcam(self):
              #Sometimes the camera freezes, if this is the case then this function will restart it.
        self.counter=0
        self.cam.stop_live_video()
        self.cam.close()
        self.cam = instrument('uc480')
        self.exposure(float(self.pics)*self.exp)
        print("Reset Camera")
        return
 
    def autoexposure(self):
              #Perform an autoexposure: sets the exposure to the result.
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
              #Performs an autoexposure on a hexagonal region of the image specified by: cent (the centre of the hexagon),
              #dist (the radius of the hexagon) and angle (it's orientation).
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
              #Sets the camera exposure to the given value.
        self.pics=int(np.ceil(expo/self.maxexp))
        self.exp=expo/self.pics
        self.cam.start_live_video(exposure_time="{} millisecond".format(str(self.exp)))
        return