#
# electron_equations
#
"""
This is a collection of functions that connects Gamma, Beta,
kinetic energy, total energy and wavelength of an accelerated electron

They should be all correct, I double checked ;)
"""

from pylab import *

e= 1.60217657e-19
m= 9.10938291e-31
c=2.99792458e8
h=6.62606957e-34
sigma2fwhm = 2.355
fwhm2sigma =1/2.355
keV = 1e3; um=1e-6;

# keep in mind, kinetic energy KE=qV, where q is the charge of particle
# and V is acceleration voltage
# for convenience, all inputs and outputs of KE
# should be in eV rather than joules

# beta = sqrt(1.0-1.0/gamma**2)
def gamma2beta(gamma):
    return sqrt(1.0-1.0/gamma**2)
# gamma = 1.0/sqrt(1.0-beta**2)
def beta2gamma(beta):
    return 1.0/sqrt(1.0-beta**2)

# gamma = 1.0+ KE*e/(m*c**2)
def KE2gamma(KE):
    return 1.0+(KE*e)/(m*c**2)
# KE = e*(gamma-1)*m*c**2
def gamma2KE(gamma):
    return e*(gamma-1)*m*c**2

def KE2beta(KE):
    return  gamma2beta(KE2gamma(KE))

# Total energy = e*gamma*m*c**2
def gamma2energy(gamma):
    return e*gamma*m*c**2

# wavelength = h*c/sqrt((KE*e+m*c**2)**2 -m**2 *c**4)
def KE2wavelength(KE):
    return h*c/sqrt((KE*e+m*c**2)**2 -m**2 *c**4)
# KE = sqrt((h*c/wavelength)**2 +(m*c**2)**2) -m*c**2
def wavelength2KE(wavelength):
    return sqrt((h*c/wavelength)**2 +(m*c**2)**2) -m*c**2

#-- stdz to pulse duration ---#
def stdz2duration(KE,L):
    v= KE2beta(KE)*c
    dt = L/v
    return dt

#-- electron charge to electron number
def coulomb2number(Q):
    return Q/e

#--- Wien filter equations ---#
def fieldlength (B, KE, lengthfactor):
    avgv= KE2beta(KE)*c
    omega_c = e*B/m
    Period = pi/omega_c 
    return Period*avgv*lengthfactor

#--- relativistic E and B change
def Ex_prime(Ex,By,KE):
    avgv = KE2beta(KE)*c
    G = KE2gamma(KE)
    return G*(Ex - avgv*By)

def By_prime(Ex,By,KE):
    avgv = KE2beta(KE)*c
    G = KE2gamma(KE)
    return G*(By - avgv*Ex/c**2)    

def meters2picosecond(length):
    return (length/c)*1e12

def picosecond2meters(ps):
    return ps*(1e-12)*c

def radius2theta(radius, detector_distance):  # assuming small angle and return in radians
    return radius/detector_distance

def radius2dhkl(radius, detector_distance, wavelength):
    return 2.0*wavelength*radius/detector_distance






