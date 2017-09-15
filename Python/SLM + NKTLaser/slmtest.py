 # -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:28:06 2016

@author: Hera
"""

import slmpy
import numpy as np
slm = slmpy.SLMdisplay(monitor = 1)
resX, resY = slm.getSize()
print resX
print resY
slmout=np.zeros((resY,resX))

#slmactive=np.ones((512,512))*255
#slmblank=np.zeros((512,256))
#slmactive[0:512,0:256]=slmblank
X,Y = np.meshgrid(np.linspace(0,512,512),np.linspace(0,512,512))
testIMG = np.round((2**8-1)*(0.5+0.5*np.sin(2*np.pi*X/50))).astype('uint8')
slmout[0:512,0:512]=testIMG


slmout=slmout.astype('uint8')
#print slmout[0:512,0:512]
slm.updateArray(slmout)
#slm.close()