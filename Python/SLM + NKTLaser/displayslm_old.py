#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:06:14 2016

@author: Andrei
"""

import numpy as np
import scipy.misc
from functions import SLMinit, useSLM
from beamshapes import LG, square, gausscomp, nothing, pic, test, angled, displace, rotator
intaim=np.load('calib.npy')


slmpix=[512,512]
spfrq=0.15
wavelength=651
wavelength=625

#distribution=nothing(slmpix)

distribution=np.real(LG(9,2,100,slmpix))*intaim
#distribution=angled(distribution, 2.09, -6.28, slmpix)
#distribution=displace(distribution, -33, -17, slmpix)*intaim
angle=24
#distribution=rotator(np.real((LG(3,2,100,slmpix))),angle)*intaim

#distribution=(pic(352,0.37,slmpix))
#angle=-np.pi*0.65+np.pi/2
angle=27
angle=16
distribution=rotator(np.real((LG(0,0,200,slmpix))),angle)*intaim
#distribution=rotator(np.real((LG(3,2,70,slmpix))),angle)#x`*intaim
#distribution=intaim
#distribution=rotator(np.real((LG(6,3,15,slmpix))),angle)*intaim
#distribution=rotator(np.real((LG(3,2,20,slmpix))),angle)*intaim
#distribution=rotator(np.real((LG(3,2,100,slmpix))),angle)*intaim
#distribution=rotator(np.real((LG(0,0,160,slmpix))),angle)*intaim
#distribution=angled(distribution, -2, 2, slmpix)
#distribution=displace(distribution, -10, 6, slmpix)*intaim
#angle=-np.pi*0.2
#distribution=np.real((LG(3,2,70,slmpix))*np.exp(1.j*angle))
#angle=-np.pi*0.65
#distribution=np.real((LG(6,6,160,slmpix))*np.exp(1.j*angle))
#angle=-np.pi*0.5
#distribution=np.real((LG(2,0,80,slmpix))*np.exp(1.j*angle))
distribution=square(512,slmpix)*intaim
distribution=square(512,slmpix)
#distribution=square(50,slmpix)*intaim
#distribution=square(200,slmpix)
#distribution=rotator(pic('test.csv',0.8,slmpix),57)
#distribution=nothing(slmpix)
#distribution=test(70,slmpix)
#separator=separate(slmpix,spfrq)

#holo=genhol(distribution,separator,inverter)
#distribution=np.fliplr(pic2(0,1, slmpix))*intaim
#tosend=SLMprocess(holo,wavelength)
slm, slmstate =SLMinit(1)

#SLMsend(tosend,slmstate,slm)
#time.sleep(0.5)

tosend=useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
#slm.close()

tosend=tosend.astype('uint16')
scipy.misc.imsave('hologram.png', tosend[:,:,0]+np.left_shift(tosend[:,:,1],8))
scipy.misc.imsave('amp.png', np.abs(distribution))
scipy.misc.imsave('phase.png', np.angle(distribution))













