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
from beamclass import beamshapes
from slmclass import SLM
from cameraclass import camera
import time

#setup data
cam=camera(10)
slmpix=[512,512]
spfrq=0.15
wavelength=651
squant=np.array([500,500])
slm=SLM(slmpix, wavelength, spfrq, 1)
beamgen=beamshapes(slmpix)

#setup
distribution=beamgen.square(squant[0])
slm.useSLM(distribution)
time.sleep(1)
print(cam.autoexposure())
#cam.exposure(10)
image=cam.takeimagebw()
print('Max pixel=',np.amax(image))
cent, dist, angle = areainterest(image, 4)
distribution=beamgen.nothing()
slm.useSLM(distribution)
print(cam.autoexposure())
image=cam.takeimagebw()
scipy.misc.imsave('before.png', image)   

mask=np.zeros(slmpix)
#make mask for SLM
mask[(slmpix[0]-squant[0])/2:(slmpix[0]+squant[0])/2,(slmpix[1]-squant[1])/2:(slmpix[1]+squant[1])/2]=np.ones(squant)
oldist=np.ones(slmpix)

#iterate 5 times to get flatter distribution
for i in range(0,10):
    #prints exposure value to console
    print(cam.autoexposure())
    image=cam.takeimagebw()
    #crop image to only contain region of interest
    image=image[int(cent[1]-dist):int(cent[1]+dist), int(cent[0]-dist):int(cent[0]+dist)]
    scipy.misc.imsave('{}.png'.format(i), image)
    #add 45 as angle is to corner?
    image=rotate(image,angle/(2*np.pi)*360+45,reshape=False)
    size=np.shape(image)
    image=image[int(size[0]/2-dist/np.sqrt(2)):int(size[0]/2+dist/np.sqrt(2)), int(size[1]/2-dist/np.sqrt(2)):int(size[1]/2+dist/np.sqrt(2))]
    zoomf=squant.astype(float)/(np.array(np.shape(image))).astype(float)
    image= zoom(image, zoomf)
    #flips image vertically
    image=(np.flipud(image))
    #flips image horizontally (added as flatening has lr problem, if doesn't work remove this and vertical flip?)
    #flips relate to transpose further down?!
    #image=(np.fliplr(image))
    #apply gaussian filter to smooth image
    image=scipy.ndimage.filters.gaussian_filter(image.astype(float),(10,10))
    scipy.misc.imsave('{}c.png'.format(i), image)
    #update target intensity of SLM
    image=1/np.sqrt(image)
    image=image/(np.amax(image))
    intaim=np.ones(slmpix)
    intaim[(slmpix[0]-squant[0])/2:(slmpix[0]+squant[0])/2,(slmpix[1]-squant[1])/2:(slmpix[1]+squant[1])/2]=image
    #transpose (change depending on +ve/-ve spfrq)
    intaim=np.transpose(intaim)    
    intaim=oldist*intaim
    intaim=intaim/np.amax(intaim*mask)*mask-(mask-1)
    oldist=intaim
    #put updated mask on to SLM
    slm.useSLM(intaim)


slm.useSLM(intaim)
cam.autoexposure()
image=cam.takeimagebw()
#save flattened distribution image & calibration
scipy.misc.imsave('after.png', image)
np.save('calib.npy', intaim)

#close connections
#slm.closeslm()
cam.closecam()
