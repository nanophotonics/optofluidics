# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 13:35:26 2018

@author: Ana Andres
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from MatplotlibSettings import *

class Ellipsoid():
    def __init__(self, a=40., b=40., c=40., num=15515, units='nm'):
        self.a = a # radius 1
        self.b = b # radius 2
        self.c = c  # radius 3
        self.num = int(num) # number of dipoles
        self.units =  units
        
        self.volume()
        self.set_shpar()        
        self.print_info()
    
    def volume(self):
        # volume
        self.V = 4./3.*np.pi*self.a*self.b*self.c
        # effective radius
        self.reff = (self.a*self.b*self.c)**(1./3.)
        
    def set_shpar(self):
        # dipole spacing
        self.s = (self.V/self.num)**(1./3.) 
        self.shpar1 = 2.*self.a/self.s
        self.shpar2 = 2.*self.b/self.s
        self.shpar3 = 2.*self.c/self.s
    
    def print_info(self):
        print "ELLIPSOID:"
        print "%f = a (%s)" % (self.a, self.units)
        print "%f = b (%s)" % (self.b, self.units)
        print "%f = c (%s)" % (self.c, self.units)
        print "%f = reff (%s)" % (self.reff, self.units)
        print "%f = s (%s^(-1))" % (self.s, self.units)
        print "%f = num dipoles" % self.num
        print "%f = shpar1" % self.shpar1
        print "%f = shpar2" % self.shpar2
        print "%f = shpar3" % self.shpar3
        print
    
        
class Cylndrcap():
    def __init__(self, L = 180, d = 45., num = 25000, units='nm'):
        self.L = L # length of the cylinder including end-caps
        self.d = d # nanorod diameter
        self.num = num # number of dipoles
        self.units = units

        self.set_l()        
        self.set_AR()
        
        self.volume()
        self.set_shpar()
        self.print_info()
    
    def volume(self):
        # volume
        self.V = 4./3.*np.pi*(self.d/2)**3 + np.pi*(self.l)*(self.d/2)**2
        # effective radius
        self.reff = (3.*self.V/4./np.pi)**(1./3.)
    
    def set_L(self):
        # length of the cylinder including end-caps
        self.L = self.l + self.d 
    
    def set_l(self):
        # length of the cylinder not including end-caps
        self.l = self.L - self.d 
        
    def set_AR(self):
        # aspect ratio
        self.AR = self.L/self.d
#        self.AR = self.l/self.d

    def set_shpar(self):
        # dipole spacing
        self.s = (self.V/self.num)**(1./3.) 
        self.shpar1 = self.l/self.s
        self.shpar2 = self.d/self.s
        
    def print_info(self):
        print "CYLNDRCAP:"
        print "%f = L (%s)" % (self.L, self.units)
        print "%f = d (%s)" % (self.d, self.units)
        print "%f = AR (L/R)" % self.AR
#        print "%f = AR (l/R)" % self.AR
        print "%f = reff (%s)" % (self.reff, self.units)
        print "%f = s (%s^(-1))" % (self.s, self.units)
        print "%f = num dipoles" % self.num
        print "%f = shpar1" % self.shpar1
        print "%f = shpar2" % self.shpar2
        print

class QTable():
    def __init__(self, file_path):
        self.file_path = file_path
        self.read_qtable(file_path)

    def read_qtable(self, file_path):
        file_object = open(file_path, 'r')
        
        line = file_object.readline() # DDSCAT
        self.ddscat_version = line.split('---')[1].strip()
        line = file_object.readline() # TARGET
        self.target = line.split()[2].strip(',')
        self.line = line
        line = file_object.readline() # method off fft solution
        self.fft_method = line.split()[0]
        line = file_object.readline() # prescription for dda polarisabilities
        self.dda_method = line.split()[0]
        line = file_object.readline() # shape
        self.shape = line.split()[0]
        
        line = file_object.readline() # number of dipoles
        self.num_dipoles = int(line.split()[0])
        line = file_object.readline() # refractive index file
        self.n = line.split()[0]
        
        line = file_object.readline() # beta
        self.beta_min = float(line.split()[0])
        self.beta_max = float(line.split()[1])
        self.nbeta = float(line.split()[-1])
        
        line = file_object.readline() # theta
        self.theta_min = float(line.split()[0])
        self.theta_max = float(line.split()[1])
        self.ntheta = float(line.split('=')[-1])
        
        line = file_object.readline() # phi
        self.phi_min = float(line.split()[0])
        self.phi_max = float(line.split()[1])
        self.nphi = float(line.split()[-1])        
        
        line = file_object.readline() # ETASCA parameter controlling # of scatt. dirs used to calculate <cos> etc.
        self.etasca = float(line.split()[0])
        
        line = file_object.readline() # target orientations
        self.num_target_orien = int(line.split()[3])
        line = file_object.readline() # indicent polarisations
        self.num_incident_pol = int(line.split()[1])
        
        line = file_object.readline() # table headers
        column_names = line.split()
        
        # table
        self.table = pd.read_table(file_path, skiprows=14, names=column_names, index_col=False, delim_whitespace=True)
    
    def print_info(self):
        print self.file_path
        print "%s = DDSCAT version" % self.ddscat_version
        print "%s = Target" % self.target
        print "%s = Method of FFT solution" % self.fft_method
        print "%s = Prescription for DDA polarisabilities" % self.dda_method
        print "%s = Shape" % self.shape
        print "%i = Number of dipoles" % self.num_dipoles
        print "%s = Target refractive index data" % self.n
        print "%f, %f, %i = beta_min, beta_max, nbeta" % (self.beta_min, self.beta_max, self.nbeta)
        print "%f, %f, %i = theta_min, theta_max, ntheta" % (self.theta_min, self.theta_max, self.ntheta)
        print "%f, %f, %i = phi_min, phi_max, nphi" % (self.phi_min, self.phi_max, self.nphi)
        print "%s = ETASCA param. controlling # of scatt. dirs used to calculate <cos> etc." % self.etasca
        print "%i = Number of target orientations" % self.num_target_orien
        print "%i = Number of incident polarisations" % self.num_incident_pol
        print
        
    def plot_q(self, ax):
        self.table.plot(x='wave', y=['Q_ext','Q_abs','Q_sca'], 
                        style='.-', markersize=20, ax=ax)
        title = "%s, %i dipoles" % (qtable.shape, qtable.num_dipoles)
        ax.set_title(title)
        
###############################################################################
if __name__ == '__main__':
        
#    e = Ellipsoid()
    c = Cylndrcap() 
    
    folder_paths = [
                    'C:/Users/Ana Andres/Desktop/DDSCAT/Au NR AR 2/',
                    'C:/Users/Ana Andres/Desktop/DDSCAT/Au NR AR 3/',
                    'C:/Users/Ana Andres/Desktop/DDSCAT/Au NR AR 4/',
                    ]

    file_name = 'qtable'                    
    fig = plt.figure()
    for i, folder_path in enumerate(folder_paths):
        file_path = folder_path + file_name
        qtable = QTable(file_path)
        qtable.print_info()
        subplot = int(str(len(folder_paths)) + '1' + str(i+1))
        ax = fig.add_subplot(subplot)
        qtable.plot_q(ax=ax)
        ax.set_title(folder_path)
        ax.set_xlim(0.3,0.9)
    
    
