# -*- coding: utf-8 -*-
"""
Created on Fri Dec 02 14:35:49 2016

@author: Hera
"""
from camera import initcam, autoexposure, closecam, takeimagebw, takeimagecolour
from functions import SLMinit, useSLM
from beamshapes import LG, square, gausscomp, nothing, corrector
import numpy as np
def autoflatten(slm,cam,spfrq,wavelength,slmstate,slmpix):
    sigmax=400
    sigmin=100
    while (sigmax-sigmin>=10):
        useSLM(gausscomp((sigmax+sigmin)/2,slmpix),spfrq,wavelength,slm,slmstate,slmpix)
        a=takeimagebw(cam)
        a.astype(float)
        cent=((np.shape(a))[0]/2,(np.shape(a))[1]/2)
        #edge=np.mean([np.sum(a[0:20,0:20]),np.sum(a[0:20,-20:]),np.sum(a[-20:,0:20]),np.sum(a[-20:,-20:])])
        #edge=np.mean([np.sum(a[cent[0]-10:cent[0]+10,0:20]),np.sum(a[cent[0]-10:cent[0]+10,-20:]),np.sum(a[-20:,cent[1]-10:cent[1]+10]),np.sum(a[0:20,cent[1]-10:cent[1]+10])])
        #edge=np.mean([np.sum(a[0:20,cent[1]-10:cent[1]+10]),np.sum(a[-20:,cent[1]-10:cent[1]+10])])
        edge=np.mean([np.sum(a[cent[0]-10:cent[0]+10,0:20]),np.sum(a[cent[0]-10:cent[0]+10,-20:])])
        mid=np.sum(a[cent[0]-10:cent[0]+10,cent[1]-10:cent[1]+10])
        if mid>=edge:
            sigmax=(sigmax+sigmin)/2
        else:
           sigmin=(sigmax+sigmin)/2
    return (sigmin+sigmax)/2
   
slmpix=[512,512]
spfrq=0.17
wavelength=625
slm, slmstate =SLMinit(1)
    
cam=initcam(20)
autoexposure(cam)
sigma=autoflatten(slm,cam,spfrq,wavelength,slmstate,slmpix)
print(sigma)
useSLM(corrector(nothing(slmpix),sigma,slmpix),spfrq,wavelength,slm,slmstate,slmpix)
slm.close()
closecam(cam)
