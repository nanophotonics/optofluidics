import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.integrate as integrate
################################# Define Equations here #################################
def ZR(w0,l):
	rl = np.pi*w0**2/(l)
	return rl

def wz(z,w0,l):
	w = w0*(np.sqrt(1+(z/ZR(w0,l))**2))
	return w

def R(z, w0, l): 
	r= z*(1 + ((np.pi*w0**2)/(l*z))**2)
	return r 

def NA(l,w0):
	n_aperture = nmed(l)*(l)/(np.pi*w0)
	return n_aperture

def n_particle(l):
	wavelength,n,ik = np.loadtxt("ri.txt", unpack=True, skiprows=1)
	fit_n = interp1d(wavelength,n, kind='cubic')
	fit_ik = interp1d(wavelength,ik, kind='cubic')
	ri_particle = fit_n(l*pow(10,9))+1j*fit_ik(l*pow(10,9))
	return ri_particle

def nmed(l):
	wl,n = np.loadtxt("water.txt", unpack=True, skiprows=1)
	fit_w = interp1d(wl,n, kind='cubic')
	water_ri = fit_w(l*pow(10,9))
	return water_ri

def alpha(d,l):
	e0= 8.85418782*pow(10,-12)
	a = 4*np.pi*pow(d,3)*(n_particle(l)**2 - nmed(l)**2)/(n_particle(l)**2 + 2*nmed(l)**2)
	al = (a+1j*((np.conjugate(a)*a)/(6*np.pi) * pow(2*np.pi*nmed(l)/l,3)))
	al_si=al*(pow(10,6)/(4*np.pi*e0))
	return al

def c_ext(d,l):
	e0= 8.85418782*pow(10,-12)
	a = 4*np.pi*pow(d,3)*(n_particle(l)**2 - nmed(l)**2)/(n_particle(l)**2 + 2*nmed(l)**2)
	al = ((a.imag*(2*np.pi*nmed(l)/l))+1j*((np.conjugate(a)*a)/(6*np.pi) * pow(2*np.pi*nmed(l)/l,4)))
	return al

def Efield(r,z,w0,l,p):
	e0 = 8.85418782*pow(10,-12)
	c = 2.99792458*pow(10,8)
	E00 = np.sqrt((4*p)/(np.pi*e0*nmed(l)*c*w0**2))
	efield = E00*(w0/wz(z,w0,l))*np.exp(-(r**2)/wz(z, w0, l)**2)*np.exp(-1j*(((2*np.pi*z)/l)-np.arctan(z/ZR(w0,l))+(np.pi*r**2)/(l*R(z,w0,l))))
	return efield

#def potential(r,z,d,w0,l,p,s):
#	c = 2.99792458*pow(10,8)
#	e0 = 8.85418782*pow(10,-12)
#	kt = 1.38064852*pow(10,-23)*(273+22)
#	etott = Efield(r,z,w0,l,p)+Efield(r,(s-z),w0,l,p)
#	econj = etott*np.conjugate(etott)
#	potent = -(((e0*nmed(l)**2)/4)*alpha(d,l).real*econj)/kt
#	return potent

def f_grad(r,z,d,w0,l,p):
	c = 2.99792458*pow(10,8)
	e0 = 8.85418782*pow(10,-12)
	E0 = (4*p)/(np.pi*e0*nmed(l)*c*w0**2)
	etot = Efield(r,z,w0,l,p)*np.conjugate(Efield(r,z,w0,l,p))
	gradr, gradz = np.gradient(etot,5.5*pow(10,-9))
	fgrad = ((e0*nmed(l)**2)/4)*alpha(d,l).real*gradr
	return fgrad

def f_scat(r,z,d,w0,l,p):
	e0 = 8.85418782*pow(10,-12)
	c = 2.99792458*pow(10,8)
	E0 = np.sqrt((4*p)/(np.pi*e0*nmed(l)*c*w0**2))
	econj = Efield(r,z,w0,l,p)*np.conjugate(Efield(r,z,w0,l,p))
	fscat = e0*((nmed(l)**2)/2)*c_ext(d,l)*econj
	return fscat

################################# Define Parameters  #################################
l0= 850*pow(10,-9)
l=l0/nmed(l0)
w0 = (5.3*pow(10,-6))/2
d=60*pow(10,-9)
s=ZR(w0,l)
radial= np.linspace(-4*w0,4*w0,1000,endpoint=True)
longitudinal = np.linspace(-5*s,s*5,1000, endpoint=True)
r,z = np.meshgrid(radial,longitudinal)
p=100*pow(10,-3)
print "Rayleigh Length: ", s
print "NA: ", NA(l,w0)


################################# Cumulative #################################
#gradr =  potential(r,z,d,w0,l,p,s)
#correrted = np.nan_to_num(gradr)
#uuuu = np.cumsum(correrted, axis=1)
#print "max: ", np.nanmax(f_grad(r,z,d,w0,l,p))


################################# Plot #################################
plt.contourf(r*pow(10,6), z*pow(10,6), f_grad(r,z,d,w0,l,p) ,100,cmap=cm.hot)
font = {'weight' : 'bold', 'size'   : 22}
plt.matplotlib.rc('font', **font)
plt.ticklabel_format(style='plain', scilimits=(0,0))
plt.colorbar().set_label('Force [fN]', rotation=270, labelpad=40)
plt.xlabel('r [${\mu}m$] (Fibre Core)')
plt.ylabel('z [${\mu}m$] (1 Rayleigh Length)')
plt.show()