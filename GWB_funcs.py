
import numpy as np
import scipy.integrate as intg
import scipy.optimize as opti
import math as ma

c = 2.9979*10**(10)
G = 6.673*10**(-8)
Msun = 1.998*10**(33)
pc2cm = 3.08567758*10**(18)
yr2sec = 3600.*24.*365.25

mp = 1.6726219 * 10.**(-24) #proton mass in grams for LEdd
sigT = 6.6524587158*10**(-25) #cm^2
LEdd_Fac = 4.*ma.pi* G * mp*c/sigT 


###TRAP int
Ntrap_z = 30 #25
Ntrap_L = 30 #25
Lmx = 31.0#10.**30


Fbin(Mchrp, z, ...)

	return

def d2nOdzdMc():
	return d2ndVdM() * dVOdz * Fbin(Mchrp, z)*(1.+q)**(6./5.)/q**(3./5.)

def hc2ARG():
	return Mchrp**(5./3.) * fr**(2./3.) * d2nOdzdMc()


def h2ofF():
	return ma.pi**(2./3.)/3. * 4./ma.pi/frq^2 * np.trapz(np.trapz( hc2ARG() ,Mchrps ), zs)
