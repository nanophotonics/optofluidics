# -*- coding: utf-8 -*-
"""
Created on Fri Dec 02 15:29:42 2016

@author: Hera
"""
import numpy as np
import scipy.misc
from camera import initcam, closecam, takeimagebw, autoexposure
from scipy.ndimage.interpolation import zoom
from numpy import genfromtxt
slmpix=[512,512]
phases = genfromtxt('phases.csv', delimiter=',')
amp = genfromtxt('amplitude.csv', delimiter=',')

phases= zoom(phases, 0.7)
L=np.shape(phases)[1]
phaseaim=np.zeros(slmpix)
phaseaim[(slmpix[0]-L)/2:(slmpix[0]+L)/2,(slmpix[1]-L)/2:(slmpix[1]+L)/2]=phases

amp= zoom(amp, 0.7)
L=np.shape(amp)[1]
ampaim=np.zeros(slmpix)
ampaim[(slmpix[0]-L)/2:(slmpix[0]+L)/2,(slmpix[1]-L)/2:(slmpix[1]+L)/2]=amp


scipy.misc.imsave('test.png', ampaim)
scipy.misc.imsave('test2.png', phaseaim)

