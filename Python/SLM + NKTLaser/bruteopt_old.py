#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 10:34:45 2017

@author: Andrei
"""
import matplotlib.image as mpimg
import scipy.misc
import scipy.signal
from autofunct import optsetup, compare, cropcurrent
from scipy.ndimage.interpolation import zoom, rotate
from beamshapes import test, LG, displace, angled, rotator
from functions import SLMinit, useSLM
from cameraclass import camera
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib import cm
import cv2
import scipy.optimize
import os

'''
Put identifier mode onto SLM.
Take a piture.
'''
savefolder='23rd March/optex'
masker='32e.png'
#masker='31a.png'
masker='10a.png'
maskrot=1
wavelength=651
l=1
g=0
size=81
#rotang=27
rotang=53

if not os.path.exists(savefolder):
    os.makedirs(savefolder)

intaim=np.load('calib.npy')
cam=camera(20)
slm, slmstate =SLMinit(1)
slmpix=[512,512]
spfrq=0.15


distribution=test(70, slmpix)*intaim
useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
time.sleep(1)
#cam.autoexposure()
img=cam.takeimagebw()
mask, cutmask, hexcor=optsetup(img, masker,'cutmask.bmp', maskrot)
fig=plt.imshow(mask)
plt.show()
distribution=rotator(np.real(LG(l, g, size, slmpix)),rotang)
#angle=1.4
#distribution=np.real((LG(3,2,160,slmpix))*np.exp(1.j*angle))
scipy.misc.imsave(savefolder+'/maskres.png',mask)
xposs=0
yposs=0
xangs=0
yangs=0


trydist=displace(angled(distribution, 0, 0, slmpix),0,0,slmpix)
useSLM(trydist*intaim,spfrq,wavelength,slm,slmstate,slmpix)
time.sleep(0.1)
img=cam.takeimagebw()
a=cropcurrent(img, hexcor[0], hexcor[1], hexcor[2])
cam.autoexposurecut(hexcor[0], hexcor[1], hexcor[2])
img=cam.takeimagebw()
val, currentimg =compare(img, mask,cutmask, hexcor)
scipy.misc.imsave(savefolder+'/test.png', currentimg)
print('val',val)




def tester(posvars):
    ypos=posvars[0]
    xpos=posvars[1]
    xang=xangs
    yang=yangs
    trydist=displace(angled(distribution, yang, xang, slmpix),ypos,xpos,slmpix)
    useSLM(trydist*intaim,spfrq,wavelength,slm,slmstate,slmpix)
    time.sleep(0.1)
    img=cam.takeimagebw()

    a=cropcurrent(img, hexcor[0], hexcor[1], hexcor[2])
    if np.sum(a==255)>10:
        cam.exposure(cam.autoexposurecut(hexcor[0], hexcor[1], hexcor[2])*0.9)
        img=cam.takeimagebw()
    elif (np.amax(img)<100):
        cam.autoexposurecut(hexcor[0], hexcor[1], hexcor[2])
        img=cam.takeimagebw()
    val, currentimg =compare(img, mask,cutmask, hexcor)
    currentimg=np.expand_dims(currentimg,2)
    currentimg=np.concatenate((currentimg,currentimg,currentimg),axis=2)    
    fig.set_data(currentimg.astype(float)/255)
    plt.draw()
    print(xpos, ypos, val)
    plt.pause(0.01)
    return val
    
def tester2(posvars):
    xpos=xposs
    ypos=yposs
    yang=posvars[0]
    xang=posvars[1]
    trydist=displace(angled(distribution, yang, xang, slmpix),ypos,xpos,slmpix)
    useSLM(trydist*intaim,spfrq,wavelength,slm,slmstate,slmpix)
    time.sleep(0.2)
    img=cam.takeimagebw()

    a=cropcurrent(img, hexcor[0], hexcor[1], hexcor[2])
    if np.sum(a==255)>10:
        cam.exposure(cam.autoexposurecut(hexcor[0], hexcor[1], hexcor[2])*0.9)
        img=cam.takeimagebw()
    elif (np.amax(img)<100):
        cam.autoexposurecut(hexcor[0], hexcor[1], hexcor[2])
        img=cam.takeimagebw()

    val, currentimg =compare(img, mask,cutmask, hexcor)
    currentimg=np.expand_dims(currentimg,2)
    currentimg=np.concatenate((currentimg,currentimg,currentimg),axis=2)    
    fig.set_data(currentimg.astype(float)/255)
    plt.draw()
    print(xang, yang, val)
    plt.pause(0.01)
    return val

#print(scipy.optimize.minimize(tester, np.array([5,-20]), bounds=[(-50,50),(-50,50)]))
s=time.time()

x0, fval, grid, Jout= scipy.optimize.brute(tester, ((-15,15),(-15,15)),Ns=9, finish=None, full_output=1)#, bounds=[(-50,50),(-50,50),(-6*np.pi,6*np.pi),(-6*np.pi,6*np.pi)],method='Nelder-Mead'))
np.savetxt(savefolder+"/bruteposc.csv", Jout, delimiter=",")
np.savetxt(savefolder+"/poscgrid0.csv", grid[0,:,:], delimiter=",")
np.savetxt(savefolder+"/poscgrid1.csv", grid[1,:,:], delimiter=",")
scipy.misc.imsave(savefolder+'/joutpc.png', Jout)
yposs=x0[0]
xposs=x0[1]
print('Best posn was', x0)
print(' ')


a0, fval, grid, Jout= scipy.optimize.brute(tester2, ((-10,10),(-10,10)),Ns=6, finish=None, full_output=1)#, bounds=[(-50,50),(-50,50),(-6*np.pi,6*np.pi),(-6*np.pi,6*np.pi)],method='Nelder-Mead'))
np.savetxt(savefolder+"/bruteangc.csv", Jout, delimiter=",")
np.savetxt(savefolder+"/angcgrid0.csv", grid[0,:,:], delimiter=",")
np.savetxt(savefolder+"/angcgrid1.csv", grid[1,:,:], delimiter=",")
scipy.misc.imsave(savefolder+'/joutac.png', Jout)
yangs=a0[0]
xangs=a0[1]
print('Best angle was', a0)
print(' ')

#Ns=5 and 6

x0, fval, grid, Jout= scipy.optimize.brute(tester, ((yposs-10,yposs+10),(xposs-10,xposs+10)),Ns=10, finish=None, full_output=1)#, bounds=[(-50,50),(-50,50),(-6*np.pi,6*np.pi),(-6*np.pi,6*np.pi)],method='Nelder-Mead'))
np.savetxt(savefolder+"/bruteposf.csv", Jout, delimiter=",")
np.savetxt(savefolder+"/posfgrid0.csv", grid[0,:,:], delimiter=",")
np.savetxt(savefolder+"/posfgrid1.csv", grid[1,:,:], delimiter=",")
scipy.misc.imsave(savefolder+'/joutpf.png', Jout)
yposs=x0[0]
xposs=x0[1]
print('Best posn was', x0)
print(' ')
a0, fval, grid, Jout= scipy.optimize.brute(tester2, ((yangs-4,yangs+4),(xangs-4,xangs+4)),Ns=11, finish=None, full_output=1)#, bounds=[(-50,50),(-50,50),(-6*np.pi,6*np.pi),(-6*np.pi,6*np.pi)],method='Nelder-Mead'))
np.savetxt(savefolder+"/bruteangf.csv", Jout, delimiter=",")
np.savetxt(savefolder+"/angfgrid0.csv", grid[0,:,:], delimiter=",")
np.savetxt(savefolder+"/angfgrid1.csv", grid[1,:,:], delimiter=",")
scipy.misc.imsave(savefolder+'/joutaf.png', Jout)
yangs=a0[0]
xangs=a0[1]
print('Best angle was', a0)
print(' ')
print('Time:', time.time()-s)


distribution=displace(angled(distribution, yangs, xangs, slmpix),yposs,xposs,slmpix)
useSLM(distribution*intaim,spfrq,wavelength,slm,slmstate,slmpix)
time.sleep(0.1)
cam.autoexposure()
img=cam.takeimagebw()
val, currentimg =compare(img, mask,cutmask, hexcor)
currentimg=currentimg.astype('float')

text_file = open(savefolder+"/paras.txt", "w")
text_file.write('L: {}\n'.format(l))
text_file.write('G: {}\n'.format(g))
text_file.write('Size: {}\n'.format(size))
text_file.write('Rotang: {}\n'.format(rotang))
text_file.write('Horizontal position: {}\n'.format(xposs))
text_file.write('Horizontal angle: {}\n'.format(xangs))
text_file.write('Vertical position: {}\n'.format(yposs))
text_file.write('Vertical angle: {}\n'.format(yangs))
text_file.write('Optimised parameter: {}'.format(val))
text_file.write('Wavelength: {}'.format(wavelength))
text_file.close()
scipy.misc.imsave(savefolder+'/fullimage.png',img)
scipy.misc.imsave(savefolder+'/cropped.png',currentimg)
currentimg=np.expand_dims(currentimg,2)
currentimg=np.concatenate((currentimg,currentimg,currentimg),axis=2)
scipy.misc.imsave(savefolder+'/croppedjet.png',cm.jet(np.sum(currentimg,2)/np.amax(np.sum(currentimg,2))))
scipy.misc.imsave(savefolder+'/amp.png', np.abs(distribution))
scipy.misc.imsave(savefolder+'/phase.png', np.angle(distribution))
#print(scipy.optimize.brute(tester2, ((-15,15),(-15,15)),Ns=30, finish=None))
#print(scipy.optimize.basinhopping(tester, np.array([0,0,0,0]), niter=50, stepsize=[5,5,0.5,0.5]))#, bounds=[(-50,50),(-50,50),(-6*np.pi,6*np.pi),(-6*np.pi,6*np.pi)],method='Nelder-Mead'))
      
slm.close()
cam.closecam()
