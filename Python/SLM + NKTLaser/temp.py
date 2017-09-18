# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from LG import LG
test=LG(1,1,100,512,512)
scipy.misc.imsave('outfile.png', np.abs(test))
