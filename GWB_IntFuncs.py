
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

#### INTEGRATION ERROR TOLS
###TRAP int
Ntrap_z = 41 #25
Ntrap_L = 41 #25

Ntrp_P = 41.
Ntrp_q = 41.

Lmx = 31.0


##################
### COSMOLOGY FUNCS
###################
### Ee defined 1/(sqrt(Om(1+z)^3 + OmL))
def OneoEe(z, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		return 1./ma.sqrt( (1.+z)**3 * Om  +  OL )
	else:
		res=[]
		for i in range(len(z)):
			res.append(1./ma.sqrt( (1.+z[i])**3 * Om  +  OL ))
		return np.array(res)


def Dang(z, h, Om, OL):
	H0 = 100.*h #kms Mpc^{-1}
	DH = c/(H0 *10.**5/(10.**6*pc2cm))
	tH = DH/c

	if (type(z) is float or type(z) is np.float64):
		return 	DH/(1.+ z) * intg.quad(OneoEe, 0., z, args=(Om,OL), epsabs=0.0)[0]
	else:
		res=[]
		for i in range(len(z)):
			res.append( DH/(1.+ z[i]) * intg.quad(OneoEe, 0., z[i], args=(Om,OL), epsabs=0.0)[0] )
	return np.array(res)
	
	

def dVdzdOm(z, h, Om, OL):
	H0 = 100.*h #kms Mpc^{-1}
	DH = c/(H0 *10.**5/(10.**6*pc2cm))
	tH = DH/c
	return  DH*(1.+z)*(1.+z)*Dang(z, h, Om, OL)*Dang(z, h, Om, OL)*OneoEe(z, Om, OL)



def dtrdz(z, h, Om, OL):
	H0 = 100.*h #kms Mpc^{-1}
	return OneoEe(z, Om, OL)/(H0 * (1.+z))


def fGW_pnt(P, M, qs):
	etaGW = 5./256. * 1./(2.*ma.pi)**(8./3.) * c**5/G**(5./3.)
	return etaGW * M**(-5./3.) * P**(8./3.) * qs**(-1.)



def fGas_pnt(qs, MdEff, eps):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd
	return qs/(4.*eps)

def t_tot(amax, M, qs, MdEff, eps):
	return np.minimum(fGas_pnt(qs, MdEff, eps), fGW_pnt(PmaxNPC(amax, M), M, qs))


def PbaseObs(Pbase,z):
	return Pbase/(1.+z)


def PmaxNPC(Npc, M):
	return 2.*ma.pi*(Npc)**(3./2.)/np.sqrt(G*M)


def PminRes(M, thmn, z, h, Om, OL):
	DA=Dang(z, h, Om, OL)
	thDA = (thmn * Dang(z, h, Om, OL))
	return 2.*ma.pi*(thDA)**(3./2.)/np.sqrt(G*M)

def qsofq(q):
	return 4.*q/(1. + q)**2

def qofqs(qs):
	return (2. - 2.*np.sqrt(1.-qs) - qs)/qs

def asep(P, M):
	return (P/(2.*ma.pi))**(2./3.) * (G*M)**(1./3.)







##################
### PTA FUNCS
###################
def fGW(P):
	return 2./P

def Mchirp(M, qs):
	return M*(qofqs(qs)/( 1. + qofqs(qs) )**2 )**(3./5.)

def hPTA(P,M,qs,z,h, Om, OL):
	DL = (1.+z)**2 * Dang(z, h, Om, OL)
	return G**(5./3.)/c**4 * 8./np.sqrt(10) * Mchirp(M, qs)/DL * ( ma.pi*Mchirp(M, qs)*fGW(P) )**(2./3.)
	#per freq
	#return 2./3./fGW(P) * G**(5./3.)/c**4 * 8./np.sqrt(10) * Mchirp(M, qs)/DL * ( ma.pi*Mchirp(M, qs)*fGW(P) )**(2./3.)


##################
### PTA FUNCS
###################







def fGW_int(P, qs, M):
	etaGW = 4.*5./256. * 1./(2.*ma.pi)**(8./3.) * c**5/G**(5./3.)
	return etaGW * M**(-5./3.) * P**(8./3.) * qs**(-1.) 


def fGas_int(qs, MdEff, eps):
	return qs/(4.*eps)

def tres_int(P, qs, M, MdEff, eps, tEdd, xi):
	return np.minimum( np.minimum( fGW_int(P, qs, M), fGas_int(qs, MdEff, eps) )/(xi*tEdd), 1.0)

def hc_int(P, qs, M, z, MdEff, eps, KQ, tEdd, h, Om, OL, xi):
	#The 6 is the max number of gav radii we allow BHs to get to eachother
	if (P < 2.*ma.pi*6.**(1.5)*G*M/c/c/c):
		#print "BHs are too close bra"
		return 0.0*hPTA(P,M,qs,z, h, Om, OL)
		##havent observed for 10 years yet...
	# elif (P>10.0*yr2sec):#(P>PmaxNPC(KQ*pc2cm, M)):
	#   	return np.exp(-P/(10.*yr2sec)) * np.minimum( np.minimum( fGW_int(P, qs, M), fGas_int(qs, MdEff, eps) )/(xi*tEdd), 1.0) * hPTA(P,M,qs,z, h, Om, OL)* hPTA(P,M,qs,z, h, Om, OL)
	else:
		return dtrdz(z, h, Om, OL) * np.minimum( np.minimum( fGW_int(P, qs, M), fGas_int(qs, MdEff, eps) )/(xi*tEdd), 1.0) * hPTA(P,M,qs,z, h, Om, OL)* hPTA(P,M,qs,z, h, Om, OL)





def GWBNum_f(z, M, fGW, thMn, qmin, eps, Pbase, MdEff, xi, KQ, h, Om, OL):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd

	Porb = 2./(fGW * (1.+z))  ##Rest frame period from observed GW freq
	qss = np.linspace(qmin, 1.0, Ntrp_q)
	dq = (1.0-qmin)/Ntrp_q
	return np.trapz(  hc_int(Porb, qss, M, z, MdEff, eps, KQ, tEdd, h, Om, OL, xi), dx=dq)
	#return  hc_int(Porb, 0.1, M, z, MdEff, eps, KQ, tEdd, h, Om, OL)

def GWBDen_f(z, M, fGW, thMn, qmin, eps, KQ, MdEff, xi, h, Om, OL):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd

	Porb = 2./(fGW * (1.+z))
	qss = np.linspace(qmin, 1.0, Ntrp_q)
	dq = (1.0-qmin)/Ntrp_q
	return  np.trapz(  tres_int(Porb, qss, M, MdEff, eps, tEdd, xi), dx=dq)




def GWB_GWgas_f(z, M, fGW, thMn, qmin_EHT, qmin_POP, eps_CBD, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	
	FF = GWBNum_f(z, M, fGW, thMn, qmin_POP, eps_CBD, Pbase, MdEff, xi, KQ, h, Om, OL)/(1.0-qmin_POP)
	#Numr = GWBNum_f(z, M, fGW, thMn, qmin_EHT, eps_CBD, Pbase, MdEff, xi, KQ, h, Om, OL)
	#Dnmr = GWBDen_f(z, M, fGW, thMn, qmin_POP, eps_CBD, KQ, MdEff, xi, h, Om, OL)

	## Should only happen if picked an inconsistent mass no? (numerical enforcement of integration bounds)
	# if (Dnmr == 0.0):
	# 	FF = 0.0
	# else:
	# 	FF = Numr/Dnmr

	# if (FF>1.0):
	# 	if (FF>1.01):
	# 		print "frac>1!" # - DOES THIS HAPPEN 
	# 		print FF
	# 	FF = 1.0
	#if (FF<0):
	#	FF=0.0

	return fbin*FF#np.maximum(fbin *FF, 1.e-14)








#########################
### INTEGRATE OVER f
########################

def GWBNum_nmr(z, M, thMn, qmin, eps, Pbase, MdEff, xi, KQ, h, Om, OL):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd
	thmn = thMn


	Npc = KQ*pc2cm
	#Pbase = np.minimum(Pbase, 2.*ma.pi*(Npc)**(3./2.)/np.sqrt(G*M)*(1.+z))


	#PMin = PminRes(M, thmn, z, h, Om, OL)
	#Pbase = np.maximum(Pbase, PMin*(1.+z))	


	P_PTAmx = 2./(10.**(-9)*(1.+z))
	P_PTAmn = np.minimum(2./(10.**(-7)*(1.+z)), PmaxNPC(Npc, M))
	if (P_PTAmx<=P_PTAmn or qmin>=1.0):
		return 0.0
	else:
	
	
		Ps  = np.linspace(P_PTAmn, P_PTAmx, Ntrp_P)
		qss = np.linspace(qmin, 1.0, Ntrp_q)

		Ivar = np.meshgrid(Ps, qss) 

		#dP = ( Pbase/(1.+z) - PMin)/Ntrp_P
		
		dP = ( P_PTAmx - P_PTAmn)/Ntrp_P
		dq = (1.0-qmin)/Ntrp_q

		fmid = 2.0 * 2.0/(P_PTAmx + P_PTAmn)
		

		return  np.trapz(  np.trapz(  hc_int(Ivar[0], Ivar[1], M, z, MdEff, eps, tEdd, h, Om, OL), dx=dP, axis=0), dx=dq, axis=0)










def GWBDen_nmr(z, M, thMn, qmin, eps, KQ, MdEff, xi, h, Om, OL):

	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd
	thmn = thMn
	DA=Dang(z, h, Om, OL)

	thDA = (thMn * Dang(z, h, Om, OL))

	PMax = PmaxNPC(KQ*pc2cm, M)

	Ps  = np.linspace(0.0, PMax, Ntrp_P)
	qss = np.linspace(qmin, 1.0, Ntrp_q)

	Ivar = np.meshgrid(Ps, qss) 

	dP = PMax/Ntrp_P
	dq = (1.0-qmin)/Ntrp_q
	#return np.trapz(  np.trapz(  tres_int(Ivar[0], Ivar[1], M, MdEff, eps, tEdd), dx=dP, axis=0), dx=dq, axis=0)
	return np.trapz(  np.trapz(  tres_int(Ivar[0], Ivar[1], M, MdEff, eps, tEdd), dx=dP, axis=0), dx=dq, axis=0)




def GWB_GWgas(z, M, thMn, qmin_EHT, qmin_POP, eps_CBD, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	# Numr = FNum(z, M, thMn, qmin, eps_CBD, Pbase, MdEff, xi, KQ, h, Om, OL)
	# Dnmr = FDen(z, M, thMn, qmin, eps_CBD, KQ, MdEff, xi, h, Om, OL)

	Numr = GWBNum_nmr(z, M, thMn, qmin_EHT, eps_CBD, Pbase, MdEff, xi, KQ, h, Om, OL)
	Dnmr = GWBDen_nmr(z, M, thMn, qmin_POP, eps_CBD, KQ, MdEff, xi, h, Om, OL)

	## Should only happen if picked an inconsistent mass no? (numerical enforcement of integration bounds)
	if (Dnmr == 0.0):
		FF = 0.0
	else:
		FF = Numr/Dnmr

	if (FF>1.0):
		if (FF>1.01):
			print "frac>1!" # - DOES THIS HAPPEN 
			print FF
		FF = 1.0
	#if (FF<0):
	#	FF=0.0

	return FF#np.maximum(fbin *FF, 1.e-14)
#########################
### ^INTEGRATE OVER f
########################








def Lmm2Mbn(Lmm, Mmx, f_Edd):
	BCUV = 4.2 
	nu14 = 1.4e9
	numm = c/(0.1)
	L14 = 10.**(Lmm)*1.e7 /( (3.e11/(1.4e9))**(-0.1) )
	Lbol = BCUV * 10.**(1.5 * np.log10(nu14 * L14) - 19.0 )
	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )
	return np.minimum(Mmx, Mbn)






def GWBofLmm(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	BCUV = 4.2 ## Lbol = BC lambda L_lambda From 1+2011 at 145 nm Runnoe+2012 Table 2 https://arxiv.org/pdf/1201.5155v1.pdf 
	nu14 = 1.4e9
	numm = c/(0.1)
	#L14 = (Lmm)*1.e7/chi  #in cgs, chi converts specific Lums, not nuL_nu
	L14 = 10.**(Lmm)*1.e7 /( (3.e11/(1.4e9))**(-0.1) )##CAREUFL WE ARE CONVERTING TO L14 here, not L408 liek in Lum Func #We still want observed lumm because BCs are in terms of observred values!!!! 10** if int log L
	Lbol = BCUV * 10.**(1.5 * np.log10(nu14 * L14) - 19.0 )  ## FROM TABLE 6 in from Baldi+2014 https://arxiv.org/pdf/1405.1711v1.pdf

	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )# in units of Msun


	Mbn = np.maximum( np.minimum(Mmx, Mbn), 10.**5)  ## we integrate L to large values, but cutoff M in F
	#Mbn = np.maximum( np.minimum(Mmx*10., Mbn), 10.**0)

	return GWB_GWgas(z, Mbn*Msun, thMn, qmin_EHT, qmin_POP, eps, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)



def GWBofLmm_f(Lmm, z, fGW, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	BCUV = 4.2 ## Lbol = BC lambda L_lambda From 1+2011 at 145 nm Runnoe+2012 Table 2 https://arxiv.org/pdf/1201.5155v1.pdf 
	nu14 = 1.4e9
	numm = c/(0.1)
	#L14 = (Lmm)*1.e7/chi  #in cgs, chi converts specific Lums, not nuL_nu
	L14 = 10.**(Lmm)*1.e7 /( (numm/(nu14))**(-0.1) )##CAREUFL WE ARE CONVERTING TO L14 here, not L408 liek in Lum Func #We still want observed lumm because BCs are in terms of observred values!!!! 10** if int log L
	Lbol = BCUV * 10.**(1.5 * np.log10(nu14 * L14) - 19.0 )  ## FROM TABLE 6 in from Baldi+2014 https://arxiv.org/pdf/1405.1711v1.pdf

	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )# in units of Msun


	Mbn = np.maximum( np.minimum(Mmx, Mbn), 10.**5)  ## we integrate L to large values, but cutoff M in F


	return GWB_GWgas_f(z, Mbn*Msun, fGW, thMn, qmin_EHT, qmin_POP, eps, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)






##See Yaun, Z Wang J + 2016 Spectral index paper and https://arxiv.org/pdf/1602.04298.pdf

def smLF_Flat(Lmm, z, chi):
	#Lmm = Lmm/(10.**7) #put in SI for LF

	zsig = 0.32
	z0 = 0.91 #0.915
	mm = -0.75
	zGEQz0 = 0.5*(np.sign(z-z0) + 1.0)
	zLEQz0 = 0.5*(np.sign(z0-z) + 1.0)
	# if (z <= z0):
	# 	ee1 = z**mm * np.exp(-0.5 * ((z - z0)/zsig)**2)
	# else:
	# 	ee1 = z**mm

	ee1 = zGEQz0 * z**mm * np.exp(-0.5 * ((z - z0)/zsig)**2) + zLEQz0 * z**mm

	Lstr = 10.**(26.40)#10.**(25.86)
	Aa = 1.76
	Bb = 0.53
	#phi0 = 10.**(-5.308)
	phi1 = 10.**(-5.58)

	k1 = -0.11#-0.107
	k2 = 0.79#0.796
	ee2 = 10.**(k1*z + k2*z*z)

	#Lum1p4 = Lmm/chi
	Lum408 = 10.**Lmm/chi

	p0 = phi1 / ( (Lum408/ee2/Lstr)**Aa + (Lum408/ee2/Lstr)**Bb ) 

	#return ee1 * p0 / Lum1p4/ee2 #div by L1p4 to put it into dN/dL not dN/DlogL
	return ee1 * p0 /np.log(10.) ## in log L  is this Log or log10!!!! will be extra factor of ln10 !!

def smLF_Steep(Lmm, z, chi):
	#Lmm = Lmm/(10.**7) #put in SI for LF

	zsig = 0.320
	z0 = 0.915
	mm = -0.73
	zGEQz0 = 0.5*(np.sign(z-z0) + 1.0)
	zLEQz0 = 0.5*(np.sign(z0-z) + 1.0)
	# if (z <= z0):
	# 	ee1 = z**mm * np.exp(-0.5 * ((z - z0)/zsig)**2)
	# else:
	# 	ee1 = z**mm

	ee1 = zGEQz0 * z**mm * np.exp(-0.5 * ((z - z0)/zsig)**2) + zLEQz0 * z**mm

	Lstr = 10.**(25.86)
	Aa = 1.76
	Bb = 0.53
	phi1 = 10.**(-5.3)#10.**(-4.44) ##https://arxiv.org/pdf/1602.04298v2.pdf has -5.3 instead of -4.44

	k1 = -0.107
	k2 = 0.796
	ee2 = 10.**(k1*z + k2*z*z)

	#Lum1p4 = Lmm/chi
	Lum408 = 10.**Lmm/chi
	p0 = phi1 / ( (Lum408/ee2/Lstr)**Aa + (Lum408/ee2/Lstr)**Bb ) 

	## div by Lum becuase above is in in d/dnp.logL
	#return ee1 * p0 / Lum1p4/ee2
	return ee1 * p0 /np.log(10.) ## in lnL in Eq 8 of Yuan Wang +2016 so convert to per log10L

def smLF(Lmm, z, chi):
	return smLF_Flat(Lmm, z, chi) + smLF_Steep(Lmm, z, chi)





def Lmin(z, h, Om, OL, Fmin):
	DL = (1.+z)**2 * Dang(z, h, Om, OL)
	#alp=1.0 ##alp=1 matches the 0.75 mJy -> L_nu curve in Figure 6 of Baldi, Capetti + 2014 (BECAUSE THEY ARE PLOTTING Lrest = Lobs(1+Z- but we wan Lobs))
	#LmmOVERLnuem = 1.0#(1.0+z)**alp ## shouldn't be large for small z here but should gen be <=1 so 1 is conservative
	## K correction from Lnu ~ nu^(-0.1) and so L_nu/L_{nu(1+z)} = (1+z)^(0.1)
	return np.log10( (1.+z)**(0.1) * 4.*ma.pi * DL*DL* Fmin/(1.+z)/(1.e7) ) #in cgs, LF converts to SI there only




def GWB_Integrand_GWgas(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	return dVdzdOm(z, h, Om, OL) * smLF(Lmm, z, chi)/(10.**6 * pc2cm)**3 * GWBofLmm(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)	


def GWB_Integrand_GWgas_f(Lmm, z, fGW, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	return dVdzdOm(z, h, Om, OL) * smLF(Lmm, z, chi)/(10.**6 * pc2cm)**3 * GWBofLmm_f(Lmm, z, fGW, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)	



def GWB_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
	MdEff = 0.1
	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	return (Lmx - Lmn)/(2.*Ntrap_L) * (2.0 * np.sum([GWB_Integrand_GWgas(L, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) for L in Ls]) - GWB_Integrand_GWgas(Lmn, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) - GWB_Integrand_GWgas(Lmx, z,  Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )


def GWB_Trap_dL_f(z, fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
	MdEff = 0.1
	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	return (Lmx - Lmn)/(2.*Ntrap_L) * (2.0 * np.sum([GWB_Integrand_GWgas_f(L, z, fGW, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) for L in Ls]) - GWB_Integrand_GWgas_f(Lmn, z, fGW, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) - GWB_Integrand_GWgas_f(Lmx, z, fGW, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )


# def GWB_Trap_dL_f(z, fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
# 	MdEff = 0.1
# 	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
# 	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
# 	return np.trapz(GWB_Integrand_GWgas_f(Ls, z, fGW, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL), Ls)



def GWB_Trap_dz(zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
	zs = np.linspace(0.000001, zmax, Ntrap_z)
	#return (zmax - 0.000001)/(2.*Ntrap_z) * (2.0 * np.sum([GWB_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) for z in zs]) - GWB_Trap_dL(zs[0], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) - GWB_Trap_dL(zs[Ntrap_z-1],  Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) )
	return 4.*ma.pi * (zmax - 0.000001)/(2.*Ntrap_z) * (2.0 * np.sum([GWB_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) for z in zs]) - GWB_Trap_dL(zs[0], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) - GWB_Trap_dL(zs[Ntrap_z-1],  Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) )


def GWB_Trap_dz_f(fGW, zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
	zs = np.linspace(0.000001, zmax, Ntrap_z)
	return 4.*ma.pi * (zmax - 0.000001)/(2.*Ntrap_z) * (2.0 * np.sum([GWB_Trap_dL_f(z, fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) for z in zs]) - GWB_Trap_dL_f(zs[0], fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) - GWB_Trap_dL_f(zs[Ntrap_z-1], fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) )
	#return (zmax - 0.000001)/(2.*Ntrap_z) * (2.0 * np.sum([GWB_Trap_dL_f(z, fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) for z in zs]) - GWB_Trap_dL_f(zs[0], fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) - GWB_Trap_dL_f(zs[Ntrap_z-1], fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) )


# def GWB_Trap_dz_f(fGW, zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
# 	zs = np.linspace(0.000001, zmax, Ntrap_z)
# 	return 4.*ma.pi * np.trapz(GWB_Trap_dL_f(zs, fGW, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL), zs)



def IntzZ_Trap_GWB(p, zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL):
	print p
	eps, KQ = p
	hGWB = GWB_Trap_dz(zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL)
	print "hGWB=%g" %hGWB
	return -hGWB



def IntzZ_Trap_GWB_f(p, fGW, zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL):
	print p
	eps, KQ = p
	hGWB = np.sqrt(GWB_Trap_dz_f(fGW, zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL))
	print "hGWB=%g" %hGWB
	return -hGWB
















