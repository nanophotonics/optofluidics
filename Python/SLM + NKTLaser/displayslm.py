# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 15:46:17 2017

@author: Hera
"""
from beamclass import beamshapes
from slmclass import SLM
import numpy as np
intcalib=np.load('calib.npy')

wavelength=495
slmpix=[512,512]
spfrq=0.15
beamgen=beamshapes(slmpix)
slm=SLM(slmpix, wavelength, spfrq, 1)

#different distributions below

distribution=np.real(beamgen.LG(3,2,70)*intcalib)
#distribution=beamgen.rotator(distribution,27)

distribution=np.real(beamgen.square(512))*intcalib
##distribution=np.real(beamgen.square(50))
#distribution=np.real(beamgen.nothing())
#print(np.max(distribution))
#pust to SLM
slm.useSLM(distribution)