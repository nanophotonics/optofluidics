#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 21:07:03 2017

@author: Andrei
"""
import numpy as np
from scipy.ndimage.interpolation import zoom, rotate
import matplotlib.image as mpimg
from autofunct import areainterest
import scipy.misc
import scipy.ndimage.filters
from beamshapes import square, nothing
from functions import SLMinit, useSLM
from cameraclass import camera
import time

cam=camera(10)
slm, slmstate =SLMinit(1)
slmpix=[512,512]
spfrq=0.15
wavelength=601
wavelength=625
squant=np.array([512,512])
squant=np.array([450,450])

distribution=square(squant[0], slmpix)
useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
time.sleep(1)
print(cam.autoexposure())
image=cam.takeimagebw()
print('Max pixel=',np.amax(image))
cent, dist, angle = areainterest(image, 4)
distribution=nothing(slmpix)
useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
print(cam.autoexposure())
image=cam.takeimagebw()
scipy.misc.imsave('before.png', image)   

mask=np.zeros(slmpix)
mask[(slmpix[0]-squant[0])/2:(slmpix[0]+squant[0])/2,(slmpix[1]-squant[1])/2:(slmpix[1]+squant[1])/2]=np.ones(squant)
oldist=np.ones(slmpix)
for i in range(0,10):
    print(cam.autoexposure())
    image=cam.takeimagebw()
    image=image[int(cent[1]-dist):int(cent[1]+dist), int(cent[0]-dist):int(cent[0]+dist)]    
    scipy.misc.imsave('{}.png'.format(i), image)
    image=rotate(image,angle/(2*np.pi)*360+45,reshape=False)
    size=np.shape(image)
    image=image[int(size[0]/2-dist/np.sqrt(2)):int(size[0]/2+dist/np.sqrt(2)), int(size[1]/2-dist/np.sqrt(2)):int(size[1]/2+dist/np.sqrt(2))]
    zoomf=squant.astype(float)/(np.array(np.shape(image))).astype(float)
    image= zoom(image, zoomf)
    #correct fo lr image inversion from  mirror?
    image = np.flipud(image)
    image=scipy.ndimage.filters.gaussian_filter(image.astype(float),(10,10))
    scipy.misc.imsave('{}c.png'.format(i), image)
    image=1/np.sqrt(image)
    #image=1/image
    image=image/(np.amax(image))
    intaim=np.ones(slmpix)
    intaim[(slmpix[0]-squant[0])/2:(slmpix[0]+squant[0])/2,(slmpix[1]-squant[1])/2:(slmpix[1]+squant[1])/2]=image
    intaim=np.transpose(intaim)    
    intaim=oldist*intaim
    intaim=intaim/np.amax(intaim*mask)*mask-(mask-1)
    oldist=intaim
    useSLM(intaim,spfrq,wavelength,slm,slmstate,slmpix)


useSLM(intaim,spfrq,wavelength,slm,slmstate,slmpix)
cam.autoexposure()
image=cam.takeimagebw()
scipy.misc.imsave('after.png', image)
np.save('calib.npy', intaim)
slm.close()
cam.closecam()
