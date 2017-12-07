# -*- coding: utf-8 -*-
"""
Created on Mon Dec 04 15:00:20 2017

@author: Ana Andres
"""
import os
from scipy import misc
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np
from MatplotlibSettings import *
from lmfit.models import GaussianModel

    
def main():
    
    image_directory = 'R:/3-Temporary/sjd87/2017.11.30 - TiSa cuvette/test/'
    
    #def read_image(image_path)
    
    for index, image_name in enumerate(os.listdir(image_directory)): 
        if '.png' in image_name:
            print "Reading file: " + str(index+1)
            print image_name
            
            image_path = image_directory + image_name
            image_array = misc.imread(image_path)
            
            plt.figure()
            plt.imshow(image_array)
            plt.xlabel('x (px)')
            plt.ylabel('y (px)')
            plt.title(image_name)
            plt.colorbar()
            
            
            
            max_index = np.unravel_index(image_array.argmax(), image_array.shape)
            print max_index
    
            plt.figure()
            max_line_x = image_array[max_index[0],:]
            plt.plot(max_line_x, label='x')
            mean_x,std_x = norm.fit(max_line_x)
            x = range(image_array.shape[0])
            y = norm.pdf(x, mean_x, std_x)
            plt.plot(x, y, label='x Gaussian fit')
            
#            from numpy import loadtxt
#            from lmfit.models import GaussianModel
            
#            data = loadtxt('test_peak.dat')
            data = image_array
            x = data[:, 0]
            y = data[:, 1]
            
            mod = GaussianModel()
            
            pars = mod.guess(y, x=x)
            out  = mod.fit(y, pars, x=x)
            print(out.fit_report(min_correl=0.25))
            
            
            max_line_y = image_array[:,max_index[1]]
            
            plt.plot(max_line_y, label='y')
            plt.xlabel('pixel')
            plt.ylabel('max (grayscale)')
            plt.legend()
        
        
if __name__ == '__main__':
    main()