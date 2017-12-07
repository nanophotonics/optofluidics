# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 16:11:09 2017

@author: Ana Andres
"""
import numpy as np
import matplotlib.pyplot as plt

# CANTILEVER FROM PAPER
ro = 1420.0
E = 3.1e9
L = 0.502
h = 0.89e-3
b = 1.7e-3
g = 2.
S = b*h
I = b*h**3/12
w_max = 200

# OPTICAL FIBRE
ro = 2200.0
E = 73.0e9
L = 10.2e-3
R = 125e-6
g = 200 # this one has a huge effect and it's hard to find in the literature
S = np.pi*R**2
I = np.pi*R**4/4
w_max = 20000


# CALCULATIONS
w = np.linspace(1,w_max,10000)

c = (E*I/ro/S)**(1./2.)
q = ((w/c)**2*(1-1j*g/w))**(1./4.)
Delta = 2*(1+np.cosh(q*L)*np.cos(q*L))
A = np.absolute(2/Delta*(np.cos(q*L)+np.cosh(q*L)))

plt.plot(w,A)