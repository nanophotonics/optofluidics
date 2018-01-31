# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 14:16:59 2018

@author: Stijn, Ana, Omid
"""
import h5py
import numpy as np
import scipy.optimize as opt
from scipy.optimize import curve_fit
import matplotlib.pylab as plt


###### define gaussian ######
def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

###### define 2d gaussian fitting ######
def twoD_gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
                                       + c * ((y - yo) ** 2)))
    return g.ravel()
    
    
if __name__ == "__main__":
    data_file = h5py.File('C:\Users\Hera\Desktop\scan_test.h5', 'r')
    
    if 'wavelength_scan_0' in data_file.keys():
        data_group = data_file['wavelength_scan_0']
        
    for i in range(len(data_group.items())):
        try:
            temp_convert = np.array(data_group["image_%d" % i])
            temp_convert = (temp_convert - temp_convert.mean()).clip(min=0)
    
            # Subtracts the mean from each pixel and then sets the min to 0
            # Presumably there is a more scientific way of removing noise
    
            ###### x-direction fit ######
            x_max_loc = np.unravel_index(temp_convert.argmax(), temp_convert.shape)[1]
            # Finds the  y-coordinate of the brightest pixel
    
            temp_convert_x = temp_convert[:,x_max_loc]
            # Selects the row with the brightest pixel in it
    
            x_fit_x = np.arange(-np.argmax(temp_convert_x), len(temp_convert_x) - np.argmax(temp_convert_x))
            # Defines the x-axis range such that the maximum is at x=0
    
            x_fit_params, x_fit_pcov = curve_fit(gaus, x_fit_x, temp_convert_x)
            # Fits the x-data along the row of the brightest pixel to a 1d gaussian
    
            x_min_lim = np.int_(x_max_loc-(3.5*x_fit_params[2]))
            x_max_lim = np.int_(x_max_loc+(3.5*x_fit_params[2]))
            cropped_x = temp_convert[:,x_min_lim:x_max_lim]
            # Crops range to within 3.5 sigma of centre (current unused)
    
            ###### Check fit ######
            #print(x_fit_params)
            #plt.figure()
            #plt.plot(x_fit_x, temp_convert_x)
            #plt.plot(x_fit_x, gaus(x_fit_x, *x_fit_params))
            #plt.show()
    
            ###### y-direction fit ######
            y_max_loc = np.unravel_index(temp_convert.argmax(), temp_convert.shape)[0]
            temp_convert_y = temp_convert[y_max_loc,:]
            y_fit_y = np.arange(-np.argmax(temp_convert_y), len(temp_convert_y) - np.argmax(temp_convert_y))
            y_fit_params, y_fit_pocv = curve_fit(gaus, y_fit_y, temp_convert_y)
            y_min_lim = np.int_(y_max_loc-(3.5*y_fit_params[2]))
            # Previously included a -1 in the above line
            y_max_lim = np.int_(y_max_loc+(3.5*y_fit_params[2]))
            #print(y_fit_params)
    
            # Have removed cropping in this version
            cropped_xy = temp_convert
            #plt.imshow(cropped_xy, origin='bottom')
            #plt.show()
    
            ###### set initial guess params for 2d gaussian
            amp_g = cropped_xy.max()
            xo_g = np.unravel_index(cropped_xy.argmax(), cropped_xy.shape)[1]
            yo_g = np.unravel_index(cropped_xy.argmax(), cropped_xy.shape)[0]
            sigma_x_g = x_fit_params[2]
            sigma_y_g = y_fit_params[2]
    
            ###### create x and y indices ######
            x_size = np.shape(cropped_xy)[1]
            y_size = np.shape(cropped_xy)[0]
            x = np.linspace(0, x_size, x_size)
            y = np.linspace(0, y_size, y_size)
            x, y = np.meshgrid(x, y)
    
    
    
    
            ###### create 2d gaussian data for fitting ######
            initial_guess = (amp_g, xo_g, yo_g, sigma_x_g, sigma_y_g, 0.2, 100)
            popt, pcov = opt.curve_fit(twoD_gaussian, (x, y), cropped_xy.reshape(x_size*y_size), p0=initial_guess)
            data_fitted = twoD_gaussian((x, y), *popt)
    
            ###### Plot results ######
            fig, ax = plt.subplots(1, 1)
            ax.imshow(cropped_xy, cmap=plt.cm.jet, origin='bottom')
            ax.contour(x, y, data_fitted.reshape(y_size, x_size), 5, colors='w')
            plt.show()
    
            pixel_size = 5.3e-6
            #print 'Frame', (frame.frame_no + 1), "/", len(frames)
            print 'amplitude', popt[0]
            print 'x0', popt[1]
            print 'y0', popt[2]
            print 'sigma_x', popt[3]
            print 'sigma_y', popt[4]
            print 'tilt', popt[5]
            print 'offset', popt[6]*pixel_size
    
            print 'x MFD: ', popt[4]*pixel_size*1000, 'um'
            print 'y MFD: ', popt[3]*pixel_size*1000, 'um'
    
        except RuntimeError:
            print "Couldn't fit a Gaussian to image", (frame.frame_no+1)
            
        

