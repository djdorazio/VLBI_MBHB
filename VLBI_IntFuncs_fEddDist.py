
import numpy as np
import scipy.integrate as intg
import scipy.optimize as opti
import scipy.special as spc
import math as ma

c = 2.9979*10**(10)
G = 6.673*10**(-8)
Msun = 1.998*10**(33)
pc2cm = 3.08567758*10**(18)
yr2sec = 3600.*24.*365.25

mp = 1.6726219 * 10.**(-24) #proton mass in grams for LEdd
sigT = 6.6524587158*10**(-25) #cm^2
LEdd_Fac = 4.*ma.pi* G * mp*c/sigT 

#JET CONSTANTS
nummGHz = c/0.1/1.e9 ##1mm in GHz
gamj = 10.0 ##typical, maybe high
thobs = 1./gamj#ma.pi/4. ##Eval at thobs such that nu_SSA is max, most conservative?
ke = 1.0 #constant order unity
Delc = np.log(1.e5) ##COul log of rmax/rmin
Lam = np.log(1.e5) ##COul log of gam_e max/gam_e min


xminDraw = -5.5 
aDraw = -0.3 
x0Draw = -0.6
sigDraw = 0.3

#### INTEGRATION ERROR TOLS
###TRAP int
Ntrap_z = 41 #25
Ntrap_L = 41 #25

Ntrp_P = 41
Ntrp_q = 31

Lmx = 30.0#10.*30
#Lmx = 25.0 ##LLAGN

#Lmx = np.log(10.**28) #(per nu, per ) ## any higher doesn't change answer much, also >~LEdd for 10^10 Msun
#QUADPACK
myrel = 1.49e-6#0.001
myabs = 1.49e-6
reclim = 1
limlst = 1
maxp1 = 1
fo = 1

##COSMOLOGY
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



def asep(P, M):
	return (P/(2.*ma.pi))**(2./3.) * (G*M)**(1./3.)

def MresBase(Pbase, thmn, z, h, Om, OL):
	thDA = (thmn * Dang(z, h, Om, OL))
	return 4.*ma.pi**2 * thDA**3 * ((1.+z)/Pbase)**2/G


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


def PTrans(M, qs, MdEff, eps):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd
	return 2. * ma.pi * (16./5.)**(3./8.) * G**(5./8.)/c**(15./8.) * qs**(3./4.) * (M)**(5./8.)*(eps)**(-3./8.)

def aTrans(M, qs, MdEff, eps):
	## eps is dimensionaless if write in terms of cst pc already
	#MEdd = LEdd_Fac/(MdEff*c*c)
	#tEdd = 1./MEdd
	#eps = eps/tEdd
	return 8.3 * 10.**(-3) * qs**(1./2.) * (M/(10.**8*Msun))**(3./4.)* eps**(-1./4.) *  pc2cm


def FbinPM(P, M, qs, MdEff, eps, xi):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd
	if (P <= PTrans(M, qs, MdEff, eps)):
		tres = (fGW_pnt(P, M, qs))/(xi*tEdd)
	else:
		tres = (fGas_pnt(qs, MdEff, eps))/(xi*tEdd)
	return min(tres,1.0)


def PbaseObs(Pbase,z):
	return Pbase/(1.+z)


def PmaxNPC(Npc, M):
	return 2.*ma.pi*(Npc)**(3./2.)/np.sqrt(G*M)


def PminRes(M, thmn, z, h, Om, OL):
	DA=Dang(z, h, Om, OL)
	thDA = (thmn * Dang(z, h, Om, OL))
	return 2.*ma.pi*(thDA)**(3./2.)/np.sqrt(G*M)


def PISCO(M):
	return 2.*np.pi * (6.*G*M/c/c)**(1.5) / np.sqrt(G*M)

##################
### Doppler FUNCS
###################
def qsofq(q):
	return 4.*q/(1. + q)**2

def qofqs(qs):
	return (2. -2.*np.sqrt(1.-qs) - qs)/qs

def asep(P, M):
	return (P/(2.*ma.pi))**(2./3.) * (G*M)**(1./3.)

def vsec(P,M,qs):
	return 1./(1.+qofqs(qs)) * np.sqrt(G*M/asep(P,M))

def gams(P, M,qs):
	return 1./np.sqrt(1.-(vsec(P,M,qs)/c)**2)

def DopMax(P,M,qs, alpha):
	return (1./(gams(P,M,qs)*(1. - vsec(P,M,qs)/c )))**(3.-alpha)
##################
### Doppler FUNCS
###################




##################
### PTA FUNCS
###################
def fGW(P):
	return 2./P

def Mchirp(M, qs):
	return M*(qofqs(qs)/(1. + qofqs(qs) )**2 )**(3./5.)

def hPTA(P,M,qs,z,h, Om, OL):
	DL = (1.+z)**2 * Dang(z, h, Om, OL)
	return G**(5./3.)/c**4 * 8./np.sqrt(10) * Mchirp(M, qs)/DL * ( ma.pi*Mchirp(M, qs)*fGW(P) )**(2./3.)
##################
### PTA FUNCS
###################





##################
### Self Grav amax
###################
def RSGff(eps, M, MdEff):
	return 30.99 * 10.**3 * (2*G*M/c**2) * (M/(1.e7*Msun))**(52./45.)  * (eps/0.1)**(-22./45.) * (MdEff/0.1)**(22./45.)

def RSGes(eps, M, MdEff):
	return 12.60 * 10.**3 * (2*G*M/c**2) * (M/(1.e7*Msun))**(-26./27.) * (eps/0.1)**(-8./27.)  * (MdEff/0.1)**(8./27.)


def RSG2(amax, eps, M, MdEff):
	Resff = 4.10  * 10.**3 * (2*G*M/c**2)
	RSGes = 12.60 * 10.**3 * (2*G*M/c**2) * (M/(1.e7*Msun))**(-26./27.) * (eps/0.1)**(-8./27.)  * (MdEff/0.1)**(8./27.)
	RSGff = 30.99 * 10.**3 * (2*G*M/c**2) * (M/(1.e7*Msun))**(52./45.)  * (eps/0.1)**(-22./45.) * (MdEff/0.1)**(22./45.)

	if (amax<Resff):
		RSG = RSGes
	else:
		RSG = RSGff



	return RSG


def RSG(eps, M, MdEff):
	Resff = 4.10  * 10.**3 * (2*G*M/c**2)
	RSGes = 12.60 * 10.**3 * (2*G*M/c**2) * (M/(1.e7*Msun))**(-26./27.) * (eps/0.1)**(-8./27.)  * (MdEff/0.1)**(8./27.)
	RSGff = 30.99 * 10.**3 * (2*G*M/c**2) * (M/(1.e7*Msun))**(52./45.)  * (eps/0.1)**(-22./45.) * (MdEff/0.1)**(22./45.)


	if (RSGes < Resff):
		RSG = RSGes
	else:
		RSG = RSGff

	return RSG

##################
### Self Grav amax
###################






##################
### Synchrotron crit frquencies for computing emission region size limits.
###################
def betj(gamj):
	return np.sqrt(1.-1./gamj/gamj)

def Dopj(gamj, thobs):
	return 1./(gamj*(1. - betj(gamj)*np.cos(thobs)))

def L44(fEdd, M, ke, Delc, Lam):
	return fEdd*LEdd_Fac*M/10.**(44) * ( Delc*(1. + 2./3.*ke*Lam) )



def L44_L(Lmm, ke, Delc, Lam):
	BCUV = 50.0#4.2 
	#nu14 = 1.4e9
	numm = c/(0.1)
	#L14 = 10.**(Lmm)*1.e7 /( (3.e11/(1.4e9))**(-0.1) )
	#Lbol = BCUV * 10.**(1.5 * np.log10(nu14 * L14) - 19.0 )
	Lbol = BCUV * 10.**(Lmm)*1.e7 * numm  ## for LLAGN, UV nu_UV lum_UV is nearly the same as nu_mm Lmm so use UV correction 
	return Lbol/10.**(44.) * ( Delc*(1. + 2./3.*ke*Lam) )

#rob in pc
def nu_SSA(z, rob, fEdd, M, thobs, gamj, ke, Delc, Lam):
	# From mathemaetica notebook "mmemissregion_critera"
	#return (3.*ke**0.3333333333333333*(Dopj(gamj,thobs))**0.6666666666666666*(L44(fEdd,M,ke,Delc,Lam))**0.6666666666666666*(np.sin(thobs))**0.6666666666666666)/(gamj**0.3333333333333333*(Delc*(1. + (2.*ke*Lam)/3.))**0.6666666666666666*rob*(1. + z)* (betj(gamj))**0.6666666666666666)  
	return 3./(1.+z) * ke**(1./3.) * (Delc*(1.+2./3.*ke*Lam))**(-2./3.) * gamj**(-4./3.) * betj(gamj)**(-2./3.) * Dopj(gamj,thobs)**(2./3.) * ((1./gamj)/np.sin(thobs) )**(-1.) * np.sin(thobs)**(-1./3.) * (L44(fEdd,M,ke,Delc,Lam))**(2./3.) / rob

#rob in pc
def BofL(L44, rob, fEdd, M, thobs, gamj, ke, Delc, Lam): ##all in cgs -> Gauss
	# From mathemaetica notebook "mmemissregion_critera"
	#return (20000000000000000000000.*np.sqrt(L44))/(np.sqrt(c*np.sqrt(1. - gamj**(-2.)))*np.sqrt(Delc*(1. + (2.*ke*Lam)/3.))*pc2cm*rob)
	#return (2.*np.sqrt(L44))/(np.sqrt(c*np.sqrt(1. - gamj**(-2.)))*np.sqrt(Delc*(1. + (2.*ke*Lam)/3.))*pc2cm*rob)
	return 2. * (Delc*(1.+2./3.*ke*Lam))**(-1./2.) * 1./gamj * (c*betj(gamj) )**(-1./2.) * gamj/(pc2cm*rob) * np.sqrt(L44*10.**44)


def nu_loss(z, rob, fEdd, M, thobs, gamj, ke, Delc, Lam):
	# From mathemaetica notebook "mmemissregion_critera"
	#return 0.07/(1. + z)* gamj**(2.)*(betj(gamj))**(2.)*(Dopj(gamj, thobs))/np.sin(thobs)/(BofL(L44(fEdd, M, ke, Delc, Lam), 1, fEdd, M, thobs, gamj, ke, Delc, Lam))**3 * rob   ##this is B at 1pc as per Eqn 21 BK79
	return 0.07/(1.+z) * gamj**2. * betj(gamj)**2. * Dopj(gamj, thobs)/np.sin(thobs) * (BofL(L44(fEdd, M, ke, Delc, Lam), 1.0, fEdd, M, thobs, gamj, ke, Delc, Lam))**(-3.) * rob


def nu_SSA_L(z, rob, Lmm, thobs, gamj, ke, Delc, Lam):
	# From mathemaetica notebook "mmemissregion_critera"
	#return (3.*ke**0.3333333333333333*(Dopj(gamj,thobs))**0.6666666666666666*(L44(fEdd,M,ke,Delc,Lam))**0.6666666666666666*(np.sin(thobs))**0.6666666666666666)/(gamj**0.3333333333333333*(Delc*(1. + (2.*ke*Lam)/3.))**0.6666666666666666*rob*(1. + z)* (betj(gamj))**0.6666666666666666)  
	return 3./(1.+z) * ke**(1./3.) * (Delc*(1.+2./3.*ke*Lam))**(-2./3.) * gamj**(-4./3.) * betj(gamj)**(-2./3.) * Dopj(gamj,thobs)**(2./3.) * ((1./gamj)/np.sin(thobs) )**(-1.) * np.sin(thobs)**(-1./3.) * (L44_L(Lmm,ke,Delc,Lam))**(2./3.) / rob

#rob in pc
def BofL_L(L44, rob, thobs, gamj, ke, Delc, Lam): ##all in cgs -> Gauss
	# From mathemaetica notebook "mmemissregion_critera"
	#return (20000000000000000000000.*np.sqrt(L44))/(np.sqrt(c*np.sqrt(1. - gamj**(-2.)))*np.sqrt(Delc*(1. + (2.*ke*Lam)/3.))*pc2cm*rob)
	#return (2.*np.sqrt(L44))/(np.sqrt(c*np.sqrt(1. - gamj**(-2.)))*np.sqrt(Delc*(1. + (2.*ke*Lam)/3.))*pc2cm*rob)
	return 2. * (Delc*(1.+2./3.*ke*Lam))**(-1./2.) * 1./gamj * (c*betj(gamj) )**(-1./2.) * gamj/(pc2cm*rob) * np.sqrt(L44*10.**(44))


def nu_loss_L(z, rob, Lmm, thobs, gamj, ke, Delc, Lam):
	# From mathemaetica notebook "mmemissregion_critera"
	#return 0.07/(1. + z)* gamj**(2.)*(betj(gamj))**(2.)*(Dopj(gamj, thobs))/np.sin(thobs)/(BofL(L44(fEdd, M, ke, Delc, Lam), 1, fEdd, M, thobs, gamj, ke, Delc, Lam))**3 * rob   ##this is B at 1pc as per Eqn 21 BK79
	return 0.07/(1.+z) * gamj**2. * betj(gamj)**2. * Dopj(gamj, thobs)/np.sin(thobs) * (BofL_L(L44_L(Lmm, ke, Delc, Lam), 1.0, thobs, gamj, ke, Delc, Lam))**(-3.) * rob




def fGW_int(P, qs, M):
	etaGW = 4.*5./256. * 1./(2.*ma.pi)**(8./3.) * c**5/G**(5./3.)
	return etaGW * M**(-5./3.) * P**(8./3.) * qs**(-1.) 


def fGas_int(qs, eps):
	return qs/(4.*eps)  ## eps passed here is eps/tEdd

# def tres_int(P, qs, M, MdEff, eps, fEdd, tEdd, z):
# 	##Check if emission region is larger than binary orbit, if so make tres=0
# 	##DO WE WANT REST FRAME ASEP??
# 	##CANT DO IF STATEMENT need MINS?
# 	#res=[]
# 	for i in range(Ntrp_P):
# 		for j in range(Ntrp_q):
# 	#	for j in range(Ntrap_L):
# 		#print asep(P[i],M)/2.
# 		#if (nu_SSA(z, pc2cm*0.001, fEdd, M, thobs, gamj, ke, Delc, Lam)>nummGHz or nu_loss(z, pc2cm*0.001, fEdd, M, thobs, gamj, ke, Delc, Lam)<nummGHz):
# 			if (nu_SSA(z, asep(P[i][j],M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)>nummGHz or nu_loss(z, asep(P[i][j],M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)<nummGHz):
# 				#res.append(0.0)
# 				return 0.0*P[i][j]
# 				print "mm-emission region too big mang"
# 			else:
# 				return np.minimum( fGW_int(P[i][j], qs, M), fGas_int(qs, eps))/tEdd 
# 				#TESTING
# 				#return fGas_int(qs, eps)/tEdd 
# 				#return fGW_int(P[i][j], qs, M)/tEdd 



# 				#res.append( np.minimum( fGW_int(P, qs, M), fGas_int(qs, eps))/tEdd )
				
# 			# return np.array(res)
		

def tres_intnu(P, qs, M, MdEff, eps, fEdd, tEdd, z):
	##Check if emission region is larger than binary orbit, if so make tres=0
	## settning rest frame asep = robs max
	Tres = P*0.0
	#fEdd = 0.001
	for i in range(Ntrp_q):
		for j in range(Ntrp_P):
			if (nu_SSA(z, asep(P[i][j],M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)<=nummGHz and nu_loss(z, asep(P[i][j],M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)>=nummGHz):
				Tres[i][j] = np.minimum( fGW_int(P[i][j], qs[i][j], M), fGas_int(qs[i][j], eps))/tEdd 
			#else:
			#	print "mm too big"
	return Tres
	
	#return np.minimum( fGW_int(P, qs, M), fGas_int(qs, eps))/tEdd 



def tres_intnu_L(P, qs, M, MdEff, eps, fEdd, tEdd, z, Lmm):
	##Check if emission region is larger than binary orbit, if so make tres=0
	##DO WE WANT REST FRAME ASEP??

	Tres = P*0.0
	#fEdd = 0.001
	for i in range(Ntrp_q):
		for j in range(Ntrp_P):
			if (nu_SSA_L(z, asep(P[i][j],M)/pc2cm/2.0, Lmm, thobs, gamj, ke, Delc, Lam)<=nummGHz and nu_loss_L(z, asep(P[i][j],M)/pc2cm/2.0, Lmm, thobs, gamj, ke, Delc, Lam)>=nummGHz):
				Tres[i][j] = np.minimum( fGW_int(P[i][j], qs[i][j], M), fGas_int(qs[i][j], eps))/tEdd 
				#print "mm right size"
			# else:
			# 	Tres[i][j] = 0.0
			#	print "mm too big"
	return Tres

	#return np.minimum( fGW_int(P, qs, M), fGas_int(qs, eps))/tEdd 


def tres_int_ALL(P, qs, M, MdEff, eps, tEdd):
	##Cdont cut out those with too large an emission region (for denominator of prob)
	return np.minimum( fGW_int(P, qs, M), fGas_int(qs, eps))/tEdd



def FNum_nmr(z, M, thMn, qmin, eps, fEdd, Pbase, MdEff, xi, KQ, h, Om, OL):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd
	thmn = thMn


	Npc = KQ*pc2cm
	PMax = PmaxNPC(KQ*pc2cm, M)
	Pbase = np.minimum(Pbase, PMax*(1.+z))
	PMin = np.maximum( PminRes(M, thmn, z, h, Om, OL), PISCO(M) ) ##rest frame


	if (Pbase<=PMin*(1.+z) or qmin>=1.0):
		return 0.0
	else:
		
		Ps  = np.linspace(PMin, Pbase/(1.+z), Ntrp_P) ##integrate in rest frame
		qss = np.linspace(qmin, 1.0, Ntrp_q)

		Ivar = np.meshgrid(Ps, qss) 

		dP = ( Pbase/(1.+z) - PMin)/Ntrp_P
		dq = (1.0-qmin)/Ntrp_P
		
	 	#return intg.dblquad(tres_int, qmin, 1.0, lambda qs: PMin, lambda nu: Pbase/(1.+z),  args =(M, MdEff, eps, tEdd), epsabs=myabs, epsrel=myrel )[0]
		return np.trapz(  np.trapz(  tres_intnu(Ivar[0], Ivar[1], M, MdEff, eps, fEdd, tEdd, z), dx=dP, axis=0), dx=dq, axis=0)




def FNum_nmr_L(z, M, thmn, qmin, eps, fEdd, Pbase, MdEff, xi, KQ, h, Om, OL, Lmm):
	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd


	Npc = KQ*pc2cm
	PMax = PmaxNPC(KQ*pc2cm, M)
	Pbase = np.minimum(Pbase, PMax*(1.+z))
	PMin = np.maximum( PminRes(M, thmn, z, h, Om, OL), PISCO(M) ) ##rest frame


	if (Pbase<=PMin*(1.+z) or qmin>=1.0):
		return 0.0
	else:
		
		Ps  = np.linspace(PMin, Pbase/(1.+z), Ntrp_P) ##integrate in rest frame
		qss = np.linspace(qmin, 1.0, Ntrp_q)

		Ivar = np.meshgrid(Ps, qss) 

		dP = ( Pbase/(1.+z) - PMin)/Ntrp_P
		dq = (1.0-qmin)/Ntrp_P
		
	 	#return intg.dblquad(tres_int, qmin, 1.0, lambda qs: PMin, lambda nu: Pbase/(1.+z),  args =(M, MdEff, eps, tEdd), epsabs=myabs, epsrel=myrel )[0]
		return np.trapz(  np.trapz(  tres_intnu_L(Ivar[0], Ivar[1], M, MdEff, eps, fEdd, tEdd, z, Lmm), dx=dP, axis=0), dx=dq, axis=0)





def FDen_nmr(z, M, thMn, qmin, eps, KQ, MdEff, xi, h, Om, OL):

	MEdd = LEdd_Fac/(MdEff*c*c)
	tEdd = 1./MEdd

	eps = eps/tEdd


	PMax = PmaxNPC(KQ*pc2cm, M)
	PIsco = PISCO(M)
	if (PMax<=PIsco):
		return 0.0
	else:

		Ps  = np.linspace(PIsco, PMax, Ntrp_P)
		qss = np.linspace(qmin, 1.0, Ntrp_q)

		Ivar = np.meshgrid(Ps, qss) 

		dP = (PMax-PIsco)/Ntrp_P
		dq = (1.0-qmin)/Ntrp_q
		
		
		return np.trapz(  np.trapz(  tres_int_ALL(Ivar[0], Ivar[1], M,  MdEff, eps, tEdd), dx=dP, axis=0), dx=dq, axis=0)

def fbin_GWgas(z, M, thMn, qmin_EHT, qmin_POP, eps_CBD, fEdd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):

	Numr = FNum_nmr(z, M, thMn, qmin_EHT, eps_CBD, fEdd, Pbase, MdEff, xi, KQ, h, Om, OL)
	Dnmr = FDen_nmr(z, M, thMn, qmin_POP, eps_CBD, KQ, MdEff, xi, h, Om, OL)

	if (Dnmr == 0.0):
		FF = 0.0
	else:
		FF = Numr/Dnmr

	if (FF>1.0):
		#print "frac>1!" # - DOES THIS HAPPEN - no only if 1.000000001 roundoff
		#print FF
		FF = 1.0
	#if (FF<0):
	#	FF=0.0

	return np.maximum(fbin * FF, 1.e-14)



def fbin_GWgas_L(z, M, thMn, qmin_EHT, qmin_POP, eps_CBD, fEdd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL, Lmm):

	Numr = FNum_nmr_L(z, M, thMn, qmin_EHT, eps_CBD, fEdd, Pbase, MdEff, xi, KQ, h, Om, OL, Lmm)
	Dnmr = FDen_nmr(z, M, thMn, qmin_POP, eps_CBD, KQ, MdEff, xi, h, Om, OL)

	if (Dnmr == 0.0):
		FF = 0.0
	else:
		FF = Numr/Dnmr

	if (FF>1.0):
		#print "frac>1!" # - DOES THIS HAPPEN - no only if 1.000000001 roundoff
		#print FF
		FF = 1.0
	#if (FF<0):
	#	FF=0.0

	return np.maximum(fbin * FF, 1.e-14)




def Mbn2Lmm(Mbn):
	f_Edd = 10.**draw_fEdd(xminDraw, aDraw, x0Draw, sigDraw)
	BCUV = 50.0#4.2 ## Lbol = BC lambda L_lambda From 1+2011 at 145 nm Runnoe+2012 Table 2 https://arxiv.org/pdf/1201.5155v1.pdf 
	numm = c/(0.1)
	Lbol = Mbn*(f_Edd * LEdd_Fac ) #* Msun
	Lmm = Lbol/BCUV/1.e7/numm ## for LLAGN, UV nu_UV lum_UV is nearly the same as nu_mm Lmm so use UV correction 
	return Lmm



def Mbn2Lmm_fxd(Mbn, f_Edd):
	BCUV = 50.0#4.2 ## Lbol = BC lambda L_lambda From 1+2011 at 145 nm Runnoe+2012 Table 2 https://arxiv.org/pdf/1201.5155v1.pdf 
	numm = c/(0.1)
	Lbol = Mbn*(f_Edd * LEdd_Fac ) #* Msun
	Lmm = Lbol/BCUV/1.e7/numm ## for LLAGN, UV nu_UV lum_UV is nearly the same as nu_mm Lmm so use UV correction 
	return Lmm




def Lmm2Mbn(Lmm, Mmx, f_Edd):
	BCUV = 50.0#4.2 
	#nu14 = 1.4e9
	numm = c/(0.1)
	#L14 = 10.**(Lmm)*1.e7 /( (3.e11/(1.4e9))**(-0.1) )
	#Lbol = BCUV * 10.**(1.5 * np.log10(nu14 * L14) - 19.0 )
	Lbol = BCUV * 10.**(Lmm)*1.e7 * numm  ## for LLAGN, UV nu_UV lum_UV is nearly the same as nu_mm Lmm so use UV correction 
	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )
	#return Mbn
	return np.maximum(np.minimum(Mmx*1., Mbn), 0.0)



def Lmm2Mbn_draw(Lmm):
	f_Edd = 10.**draw_fEdd(xminDraw, aDraw, x0Draw, sigDraw)
	BCUV = 50.0#4.2 
	#nu14 = 1.4e9
	numm = c/(0.1)
	#L14 = 10.**(Lmm)*1.e7 /( (3.e11/(1.4e9))**(-0.1) )
	#Lbol = BCUV * 10.**(1.5 * np.log10(nu14 * L14) - 19.0 )
	Lbol = BCUV * 10.**(Lmm)*1.e7 * numm  ## for LLAGN, UV nu_UV lum_UV is nearly the same as nu_mm Lmm so use UV correction 
	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )
	return Mbn
	#return np.maximum(np.minimum(5.*10.**10, Mbn), 0.0)




def step(x):
    return 1 * (x > 0)

def pdf_fEdd(x, xmin, ap, x0, sig):
	#(x is logf)
 	#return (2.*(1. + a)*(np.exp(-(-x + x0)**2/(2.*sig**2)) + (-x)**a)*(1. - step(x))* step(x - xmin))/((2.*(-xmin)**a*xmin + (1. + a)*np.sqrt(2*np.pi)*sig*(spc.erf((0.7071067811865475*x0)/sig) - spc.erf((0.7071067811865475*(x0 - xmin))/sig)))*step(-xmin)*(-1. + step(xmin)))
 	return ((10.**x)**ap + 1./(np.exp( (x - x0)**2/(2.*sig**2) )*np.sqrt(2.*np.pi)*sig))/((0.43429448190325176*(1. - 10.**(ap*xmin)))/ap + 0.5*(-spc.erf(x0/(np.sqrt(2.)*sig)) + spc.erf((x0 - xmin)/(np.sqrt(2.)*sig))))


def CDF_fEdd(xcum, xmin, a, x0, sig):
	Nftrap = 41
	xs = np.linspace(xmin, xcum, Nftrap)
	return np.trapz( pdf_fEdd(xs, xmin, a, x0, sig), xs  )

def fEdd_slv_func(xcum, xmin, a, x0, sig, Uu):
	return (CDF_fEdd(xcum, xmin, a, x0, sig) - Uu)

def fEdd_slv_func2(xcum, xmin, a, x0, sig, Uu):
	return (CDF_fEdd(xcum, xmin, a, x0, sig) - Uu)*(CDF_fEdd(xcum, xmin, a, x0, sig) - Uu)

def iCDF_fEdd(xmin, a, x0, sig, Uu):
	#return opti.brentq(fEdd_slv_func, xmin, 0.0, args=(xmin, a, x0, sig, Uu))
	return opti.fmin(fEdd_slv_func2, -4.0, args=(xmin, a, x0, sig, Uu), xtol=0.01, ftol = 0.01, disp=0)[0]

def draw_fEdd(xmin, a, x0, sig):
	#np.random.seed(2)
	Uu = np.random.rand()
	return iCDF_fEdd(xmin, a, x0, sig, Uu)

def FbinofLmm(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	# draw fEdd
	f_exp = draw_fEdd(xminDraw, aDraw, x0Draw, sigDraw)
	# if (f_exp>=-3.0):
	# 	return 0.0
	#else:
	BCUV = 50.0#4.2 ## Lbol = BC lambda L_lambda From 1+2011 at 145 nm Runnoe+2012 Table 2 https://arxiv.org/pdf/1201.5155v1.pdf 
	#nu14 = 1.4e9
	numm = c/(0.1)
	
	#L14 = 10.**(Lmm)*1.e7 /( (3.e11/(1.4e9))**(-0.1) )##CAREUFL WE ARE CONVERTING TO L14 here, not L408 liek in Lum Func #We still want observed lumm because BCs are in terms of observred values!!!! 10** if int log L
	#Lbol = BCUV * 10.**(1.5 * np.log10(nu14 * L14) - 19.0 )  ## FROM TABLE 6 in from Baldi+2014 https://arxiv.org/pdf/1405.1711v1.pdf

	Lbol = BCUV * 10.**(Lmm)*1.e7 * numm  ## for LLAGN, UV nu_UV lum_UV is nearly the same as nu_mm Lmm so use UV correction 


	f_Edd = 10.**f_exp
	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )# in units of Msun
	#NOTE THAT we could allow f_Edd to be eps above, but the frac LEdd that L is doesnt have to be linked to CBD accretion rat which drives migration
	# Lmin comes in W so correct to erg/s here


	#Mbn = np.maximum( np.minimum(Mmx, Mbn), 0.)  ## we integrate L to large values, but cutoff M in F - shouldnt lumfunc take care of this?

	# if (Mbn > Mmx):
	# 	return 0.0#np.exp(-Mbn/Mx) * fbin_GWgas_L(z, Mbn*Msun, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL, Lmm)
	# else:
	# 	return fbin_GWgas(z, Mbn*Msun, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	

	return np.exp(-(Mbn/Mmx)**20) * fbin_GWgas_L(z, Mbn*Msun, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL, Lmm)

	#return  fbin_GWgas_L(z, Mbn*Msun, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL, Lmm)



def FbinofLopt(Lopt, z, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	BCopt = 10.0 #Richards+2006
	nuVbnd = c/(5.45*10**(-5))
	Lbol = BCopt * Lopt*nuVbnd*1.e7  #input from Lmin is in SI and specific Lum
	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )# in units of Msun
	return fbin_GWgas(z, Mbn*Msun, thMn, qmin, eps, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)






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
	#return smLF_Flat(Lmm, z, chi) ## FLat Spectrum are core dominated rather then jet dominated (https://ned.ipac.caltech.edu/level5/Cambridge/Cambridge1_3_1.html)
	#return smLF_Steep(Lmm, z, chi)

# vecsmLF_Flat = np.vectorize(smLF_Flat)
# vecsmLF_Steep = np.vectorize(smLF_Steep)
# def vecsmLF(Lmm, z, chi):
# 	return vecsmLF_Flat(Lmm, z, chi) + vecsmLF_Steep(Lmm, z, chi)
# 	#return vecsmLF_Flat(Lmm, z, chi)  ## FLat Spectrum are core dominated rather then jet dominated (https://ned.ipac.caltech.edu/level5/Cambridge/Cambridge1_3_1.html)




def smLF(Lmm, z, chi):
	#https://arxiv.org/pdf/1708.03087.pdf
	Lstar = 10.**(24.79)#10.**(24.79)
	phi1 = 10.**(-4.72)#10.**(-4.72) ##this is phi0?
	p1 = -1.29 ##SIGNS OF p1 nad p2 switched in YUAN PAPER!
	p2 = 6.80
	bet = 0.45 
	gam = 0.31
	zc= 0.78


	p0 = (1./(1.+zc))**p1  +  (1./(1.+zc))**p2
	e1 = p0 / ( ((1.+z)/(1.+zc))**p1 + ((1.+z)/(1.+zc))**p2 )

	k2 = -0.16  ##k1 and k2 ARE BACKWARDS IN THE YAUN PAPER
	k1 = 1.44
	e2 = 10.**(k1*z + k2*z*z)


	Lum408 = 10.**Lmm/chi
	## do in d/dLogL for integration
	L = Lum408/e2
	rho = phi1/np.log(10.) * (L/Lstar)**(-bet) * np.exp(-(L/Lstar)**gam) 
	return e1*rho











# def Lmin(z, h, Om, OL, Fmin):
# 	DL = (1.+z)**2 * Dang(z, h, Om, OL) 
# 	return 4.*ma.pi * DL*DL* Fmin/(1.+z)/(1.e7) #in cgs, LF converts to SI there only

def Lmin(z, h, Om, OL, Fmin):
	DL = (1.+z)**2 * Dang(z, h, Om, OL)
	#alp=1.0 ##alp=1 matches the 0.75 mJy -> L_nu curve in Figure 6 of Baldi, Capetti + 2014 (BECAUSE THEY ARE PLOTTING Lrest = Lobs(1+Z- but we wan Lobs))
	#LmmOVERLnuem = 1.0#(1.0+z)**alp ## shouldn't be large for small z here but should gen be <=1 so 1 is conservative
	## K correction from Lnu ~ nu^(-0.1) and so L_nu/L_{nu(1+z)} = (1+z)^(0.1)
	return np.log10( (1.+z)**(0.1) * 4.*ma.pi * DL*DL* Fmin/(1.+z)/(1.e7) ) #in cgs, LF converts to SI there only

def FInt_smLF(z, chi, Fmin, h, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		Lmn = Lmin(z, h, Om, OL, Fmin)
		return 	intg.quad(smLF, Lmn, Lmx, args=(z, chi), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] 
	else:
		res=[]
		#zup = z/z*20.0
		for i in range(len(z)):
			Lmn = Lmin(z[i], h, Om, OL, Fmin)
			res.append(  intg.quad(smLF, Lmn, Lmx, args=(z[i], chi), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
		return np.array(res)


def FInt_OLF(z, chi, Fmin_opt, h, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		Lmn = Lmin(z, h, Om, OL, Fmin)
		return 	intg.quad(OLF, Lmn, Lmx, args=(z, chi), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] 
	else:
		res=[]
		#zup = z/z*20.0
		for i in range(len(z)):
			Lmn = Lmin(z[i], h, Om, OL, Fmin)
			res.append(  intg.quad(OLF,  Lmn, Lmx, args=(z[i], chi), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
		return np.array(res)




def Fbin_NOLF_Integrand_GWgas(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	return FbinofLmm(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)/Lmm ##Not mult bu dN/DL/dV here so need to div by Lmm (just a diagnostic)	



def Fbin_Integrand_GWgas(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	return dVdzdOm(z, h, Om, OL) * smLF(Lmm, z, chi)/(10.**6 * pc2cm)**3 * FbinofLmm(Lmm, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)	


def Fbin_Integrand_GWgas_Mbn(Mbn, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	return dVdzdOm(z, h, Om, OL) * smLF(np.log10(Mbn2Lmm(Mbn)), z, chi)/(10.**6 * pc2cm)**3 * fbin_GWgas(z, Mbn, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)

#vecdVdzdOm = np.vectorize(dVdzdOm)
#vecFbinofLmm = np.vectorize(FbinofLmm)
#def vecFbin_Integrand_GWgas(Lmm, z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
#	return vecdVdzdOm(z, h, Om, OL) * vecsmLF(Lmm, z, chi)/(10.**6 * pc2cm)**3 * vecFbinofLmm(Lmm, z, Mmx, chi#, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
#vecFbin_Integrand_GWgas = np.vectorize(vecFbin_Integrand_GWgas)

def Fbin_OPTICAL_Integrand_GWgas(Lopt, z, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	return dVdzdOm(z, h, Om, OL) * OLF(Lopt, z, chi)/(10.**6 * pc2cm)**3 * FbinofLopt(Lopt, z, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)	


def NtotDZ_Integrand(Lmm, z, chi, h, Om, OL):
	return dVdzdOm(z, h, Om, OL) * smLF(Lmm, z, chi)/(10.**6 * pc2cm)**3 


# def NtotDZ_Integrand(Lmm, z, chi, h, Om, OL):
# 	return vecdVdzdOm(z, h, Om, OL) * vecsmLF(Lmm, z, chi)/(10.**6 * pc2cm)**3 
# #vecNtotDZ_Integrand = np.vectorize(NtotDZ_Integrand)


def NtotDZ_OPTICAL_Integrand(Lopt, z, chi, h, Om, OL):
	return dVdzdOm(z, h, Om, OL) * OLF(Lopt, z, chi)/(10.**6 * pc2cm)**3 


# def NtotDZ_OPTICAL_Integrand(Lopt, z, chi, h, Om, OL):
# 	return vecdVdzdOm(z, h, Om, OL) * vecOLF(Lopt, z, chi)/(10.**6 * pc2cm)**3 
# #vecNtotDZ_OPTICAL_Integrand = np.vectorize(NtotDZ_OPTICAL_Integrand)


def FtotDZ_GWgas(z, Mmx, DZ, Fmin, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		#print "here"
		Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
		return  intg.quad(Fbin_NOLF_Integrand_GWgas, Lmn, Lmx,  args=(z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#elif(type(z) is np.ndarray):
	else:
		#print "there"
		res=[]
		#zup = z/z*40.0
		print "Why are you here?"
		for i in range(len(z)):
			Lmn = np.minimum( Lmin(z[i], h, Om, OL, Fmin), Lmx)
			res.append( intg.quad(Fbin_NOLF_Integrand_GWgas, Lmn, Lmx,  args=(z[i], Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
		return np.array(res)


def NtotDZ_GWgas(z, Mmx, DZ, Fmin, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
		return 4.*ma.pi*z * intg.quad(Fbin_Integrand_GWgas, Lmn, Lmx,  args=(z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	else:
		res=[]
		#zup = z/z*40.0
		for i in range(len(z)):
			Lmn = np.minimum( Lmin(z[i], h, Om, OL, Fmin), Lmx)
			res.append( 4.*ma.pi*z[i] * intg.quad(Fbin_Integrand_GWgas, Lmn, Lmx,  args=(z[i], Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
		return np.array(res)



def NtotDZ_OPTICAL_GWgas(z, DZ, Fmin_opt, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
		return 4.*ma.pi*z * intg.quad(Fbin_OPTICAL_Integrand_GWgas, Lmn, Lmx,  args=(z, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	else:
		res=[]
		#zup = z/z*40.0
		for i in range(len(z)):
			Lmn = np.minimum( Lmin(z[i], h, Om, OL, Fmin), Lmx)
			res.append( 4.*ma.pi*z[i]* intg.quad(Fbin_OPTICAL_Integrand_GWgas, Lmn, Lmx,  args=(z[i], chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
		return np.array(res)


def NtotDZ_RLF(z, Fmin, chi, h, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
		return 4.*ma.pi*z * intg.quad(NtotDZ_Integrand, Lmn, Lmx,  args=(z, chi, h, Om, OL) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	else:
		res=[]
		#zup = z/z*40.0
		for i in range(len(z)):
			Lmn = np.minimum( Lmin(z[i], h, Om, OL, Fmin), Lmx)
			res.append( 4.*ma.pi*z[i] * intg.quad(NtotDZ_Integrand,  Lmn, Lmx,  args=(z, chi, h, Om, OL) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
		return np.array(res)


def NtotDL_RLF(z, Fmin, chi, h, Om, OL):
	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
	return intg.quad(NtotDZ_Integrand, Lmn, Lmx,  args=(z, chi, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]


def NtotDL_Trap_RLF(z, Fmin, chi, h, Om, OL):
	#MdEff = 0.1
	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	return (Lmx - Lmn)/(2.*Ntrap_L) * (2.0 * np.sum([NtotDZ_Integrand(L, z, chi, h, Om, OL)  for L in Ls]) - NtotDZ_Integrand(Lmn, z, chi, h, Om, OL) - NtotDZ_Integrand(Lmx, z, chi, h, Om, OL))

	# if (type(z) is float or type(z) is np.float64):
	# 	Lmn = Lmin(z, h, Om, OL, Fmin)
	# 	#Lmx =  Lmin(40.0, h, Om, OL, Fmin)
	# 	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	# 	Trap_sub = vecNtotDZ_Integrand(Ls[0], z, chi, h, Om, OL) + vecNtotDZ_Integrand(Ls[Ntrap_L-1], z, chi, h, Om, OL)
	# 	Trap_int = -Trap_sub + np.sum(2.* vecNtotDZ_Integrand(Ls, z, chi, h, Om, OL) )
	# 	return (Lmx-Lmn)/(2.*Ntrap_L) * (Trap_int)
	# else:
	# 	res = []
	# 	for i in range(len(z)):
	# 		Lmn = Lmin(z[i], h, Om, OL, Fmin)
	# 		#Lmx = 1.e5 
	# 		#Lmx =  Lmin(40.0, h, Om, OL, Fmin)
	# 		Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	# 		Trap_sub = vecNtotDZ_Integrand(Ls[0], z[i], chi, h, Om, OL) + vecNtotDZ_Integrand(Ls[Ntrap_L-1], z[i], chi, h, Om, OL)
	# 		Trap_int = -Trap_sub + np.sum(2.* vecNtotDZ_Integrand(Ls, z[i], chi, h, Om, OL) )
	# 		res.append((Lmx-Lmn)/(2.*Ntrap_L) * (Trap_int))
	# 	return np.array(res)		



def NtotDL_Optical_Trap_RLF(z, Fmin, chi, h, Om, OL):
	#MdEff = 0.1
	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	return (Lmx - Lmn)/(2.*Ntrap_L) * (2.0 * np.sum([NtotDZ_Optical_Integrand(L, z, chi, h, Om, OL)  for L in Ls]) - NtotDZ_Optical_Integrand(Lmn, z, chi, h, Om, OL) - NtotDZ_Optical_Integrand(Lmx, z, chi, h, Om, OL))
	# if (type(z) is float or type(z) is np.float64):
	# 	Lmn = Lmin(z, h, Om, OL, Fmin)
	# 	#Lmx =  Lmin(40.0, h, Om, OL, Fmin)
	# 	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	# 	Trap_sub = vecNtotDZ_Optical_Integrand(Ls[0], z, chi, h, Om, OL) + vecNtotDZ_Optical_Integrand(Ls[Ntrap_L-1], z, chi, h, Om, OL)
	# 	Trap_int = -Trap_sub + np.sum(2.* vecNtotDZ_Optical_Integrand(Ls, z, chi, h, Om, OL) )
	# 	return (Lmx-Lmn)/(2.*Ntrap_L) * (Trap_int)
	# else:
	# 	res = []
	# 	for i in range(len(z)):
	# 		Lmn = Lmin(z[i], h, Om, OL, Fmin)
	# 		#Lmx = 1.e5 
	# 		#Lmx =  Lmin(40.0, h, Om, OL, Fmin)
	# 		Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	# 		Trap_sub = vecNtotDZ_Optical_Integrand(Ls[0], z[i], chi, h, Om, OL) + vecNtotDZ_Optical_Integrand(Ls[Ntrap_L-1], z[i], chi, h, Om, OL)
	# 		Trap_int = -Trap_sub + np.sum(2.* vecNtotDZ_Optical_Integrand(Ls, z[i], chi, h, Om, OL) )
	# 		res.append((Lmx-Lmn)/(2.*Ntrap_L) * (Trap_int))
	# 	return np.array(res)	








def Ntot_RLF(zmax, Fmin, chi, h, Om, OL):
	return intg.quad(NtotDL_RLF, 0.0, zmax, args=(Fmin, chi, h, Om, OL), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] 

def Ntot_Optical_RLF(zmax, Fmin, chi, h, Om, OL):
	return intg.quad(NtotDL_Optical_RLF, 0.0, zmax, args=(Fmin, chi, h, Om, OL), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] 


def Ntot_Trap_RLF(zmax, Fmin, chi, h, Om, OL):
	zs = np.linspace(0.000001, zmax, Ntrap_z)
	return 4.*ma.pi * (zmax-0.000001)/(2.*Ntrap_z) * (  2.0 * np.sum([NtotDL_Trap_RLF(z, Fmin, chi, h, Om, OL) for z in zs]) - NtotDL_Trap_RLF(zs[0], Fmin, chi, h, Om, OL) -  NtotDL_Trap_RLF(zs[Ntrap_z-1], Fmin, chi, h, Om, OL)  )


	# Trap_sub = NtotDL_Trap_RLF(zs[0], Fmin, chi, h, Om, OL) + NtotDL_Trap_RLF(zs[Ntrap_z-1], Fmin, chi, h, Om, OL)
	# Trap_int = -Trap_sub + np.sum(2.* NtotDL_Trap_RLF(zs, Fmin, chi, h, Om, OL))
	# return 4.*ma.pi * (zmax-0.000001)/(2.*Ntrap_z) * (Trap_int)



def Ntot_Optical_Trap_RLF(zmax, Fmin, chi, h, Om, OL):
	zs = np.linspace(0.000001, zmax, Ntrap_z)
	return 4.*ma.pi * (zmax-0.000001)/(2.*Ntrap_z) * (  2.0 * np.sum([NtotDL_Optical_Trap_RLF(z, Fmin, chi, h, Om, OL) for z in zs]) - NtotDL_Optical_Trap_RLF(zs[0], Fmin, chi, h, Om, OL) -  NtotDL_Optical_Trap_RLF(zs[Ntrap_z-1], Fmin, chi, h, Om, OL)  )


	# Trap_sub = NtotDL_Optical_Trap_RLF(zs[0], Fmin, chi, h, Om, OL) + NtotDL_Optical_Trap_RLF(zs[Ntrap_z-1], Fmin, chi, h, Om, OL)
	# Trap_int = -Trap_sub + np.sum(2.* NtotDL_Optical_Trap_RLF(zs, Fmin, chi, h, Om, OL))
	# return 4.*ma.pi * (zmax-0.000001)/(2.*Ntrap_z) * (Trap_int)





def NtotDZ_OLF(z, DZ, Fmin_opt, chi, h, Om, OL):
	if (type(z) is float or type(z) is np.float64):
		Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
		return 4.*ma.pi*z * intg.quad(NtotDZ_OPTICAL_Integrand, Lmn, Lmx,  args=(z, chi, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	else:
		res=[]
		#zup = z/z*40.0
		for i in range(len(z)):
			Lmn = np.minimum( Lmin(z[i], h, Om, OL, Fmin), Lmx)
			res.append( 4.*ma.pi*z * intg.quad(NtotDZ_OPTICAL_Integrand,  Lmn, Lmx,  args=(z, chi, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
		return np.array(res)




def OptNEHT(p, Mmx, DZ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL):
	print p
	z, eps, KQ = p
	#z, eps, KQ, MdEff = p
	MdEff = 0.1
	#f_Edd = eps
	if (z<0.0 or eps<0.0 or KQ<=0.0 or (MdEff > 1.0 or MdEff < 0.0)):
		NEHT = -10000000000.
	else:
		Lmn = Lmin(z, h, Om, OL, Fmin)
		NEHT = 4.*ma.pi*z * intg.quad(Fbin_Integrand_GWgas, Lmn, Lmx,  args=(z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	print "NEHT=%g" %NEHT
	return -NEHT


def OptNEHT_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL):
	MdEff = 0.1
	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
	NEHT = intg.quad(Fbin_Integrand_GWgas, Lmn, Lmx,  args=(z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	return NEHT

def OptNEHT_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
	MdEff = 0.1
	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	return (Lmx - Lmn)/(2.*Ntrap_L) * (2.0 * np.sum([Fbin_Integrand_GWgas(L, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) for L in Ls]) - Fbin_Integrand_GWgas(Lmn, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) - Fbin_Integrand_GWgas(Lmx, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )


# #vecLmin = np.vectorize(Lmin)
# def OptNEHT_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
# 	MdEff = 0.1
# 	Lmn = np.minimum( Lmin(z, h, Om, OL, Fmin), Lmx)
# 	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
# 	return (Lmx - Lmn)/(2.*Ntrap_L) * (2.0 * np.sum([Fbin_Integrand_GWgas(L, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) for L in Ls]) - Fbin_Integrand_GWgas(Lmn, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) - Fbin_Integrand_GWgas(Lmx, z, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )

	# if (type(z) is float or type(z) is np.float64):
	# 	Lmn = Lmin(z, h, Om, OL, Fmin)
	# 	#Lmx =  Lmin(40.0, h, Om, OL, Fmin)
	# 	Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	# 	Trap_sub = vecFbin_Integrand_GWgas(Lmn, z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) + vecFbin_Integrand_GWgas(Lmx, z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	# 	Trap_int = -Trap_sub + sum(2.* vecFbin_Integrand_GWgas(Ls, z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )
	# 	# Trap_int = 0.0
	# 	# for i in range(Ntrap_L):
	# 	# 	Trap_int += 2.* vecFbin_Integrand_GWgas(Ls[i], z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) 

	# 	# Trap_int -= Trap_sub 
	# 	return (Lmx-Lmn)/(2.*Ntrap_L) * (Trap_int)
	# else:
	# 	res = []
	# 	for i in range(len(z)):
	# 		Lmn = Lmin(z[i], h, Om, OL, Fmin)
	# 		#Lmx =  Lmin(40.0, h, Om, OL, Fmin)
	# 		Ls = np.linspace(Lmn,Lmx, Ntrap_L)
	# 		Trap_sub = vecFbin_Integrand_GWgas(Lmn, z[i], Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) + vecFbin_Integrand_GWgas(Lmx, z[i], Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	# 		Trap_int = -Trap_sub + sum(2.* vecFbin_Integrand_GWgas(Ls, z[i], Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )
	# 		# Trap_int = 0.0
	# 		# for j in range(Ntrap_L):
	# 		# 	Trap_int += 2.* vecFbin_Integrand_GWgas(Ls[j], z[i], Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) 
	# 		# Trap_int -= Trap_sub

	# 		res.append((Lmx-Lmn)/(2.*Ntrap_L) * (Trap_int))
	# 	return np.array(res)

# Ls = np.zeros([Ntrap_z, Ntrap_L])
# for j in range(len(zs)):
# 	Ls[:][j] = np.linspace( Lmin(zs[j], h, Om, OL, Fmin), Lmin(20.0, h, Om, OL, Fmin), Ntrap_L)

#Ls[i] sums over diff L's
#T.Ls indexies diff zs


##NOTE working, but mayeb faster:
# def OptNEHT_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL):
# 	MdEff = 0.1
# 	if (type(z) is float or type(z) is np.float64):
# 		Lmn = Lmin(z, h, Om, OL, Fmin)
# 		Lmx =  Lmin(20.0, h, Om, OL, Fmin)
# 		Ls = np.linspace(Lmn,Lmx, Ntrap_L)
# 		Trap_sub = vecFbin_Integrand_GWgas(Ls[0], z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) + vecFbin_Integrand_GWgas(Ls[Ntrap_L-1], z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
# 		Trap_int = -Trap_sub + np.sum(2.* vecFbin_Integrand_GWgas(Ls, z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )
# 		return (Lmx-Lmn)/(2.*Ntrap_L) * (Trap_int)
# 	else:
# 		Ls = np.zeros([Ntrap_z, Ntrap_L])
# 		for j in range(len(z)):
# 			Ls[:][j] = np.linspace( Lmin(z[j], h, Om, OL, Fmin), Lmin(20.0, h, Om, OL, Fmin), Ntrap_L)

# 		Trap_sub = vecFbin_Integrand_GWgas(Ls[0], z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) + vecFbin_Integrand_GWgas(Ls[Ntrap_L-1], z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
# 		Trap_int = -Trap_sub + np.sum(2.* vecFbin_Integrand_GWgas(Ls, z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL), axis=1 )
# 		return (Lmin(20.0, h, Om, OL, Fmin)-np.transpose(Ls)[0])/(2.*Ntrap_L) * (Trap_int)




def OptNEHT_dz(zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL):
	#f_Edd = eps
	return 4.*ma.pi * intg.quad(OptNEHT_dL, 0.0, zmax,  args=(Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]



# def OptNEHT_Trap_dz(zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
# 	#cant mak zminany smaller or canta invert CDF for fEdd_draw
# 	zs = np.linspace(0.000001, zmax, Ntrap_z)
# 	return 4.*ma.pi * (zmax - 0.000001)/(2.*Ntrap_z) * (2.0 * np.sum([OptNEHT_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) for z in zs]) - OptNEHT_Trap_dL(zs[0], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) - OptNEHT_Trap_dL(zs[Ntrap_z-1], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) )


def OptNEHT_Trap_dz(zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL):
	zs = np.linspace(0.000001, zmax, Ntrap_z)
	return 4.*ma.pi * (zmax - 0.000001)/(2.*Ntrap_z) * (2.0 * np.sum([OptNEHT_Trap_dL(z, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) for z in zs]) - OptNEHT_Trap_dL(zs[0], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) - OptNEHT_Trap_dL(zs[Ntrap_z-1], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL) )


	# Trap_sub = OptNEHT_Trap_dL(zs[0], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL) + OptNEHT_Trap_dL(zs[Ntrap_z-1], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL)
	# Trap_int = -Trap_sub + sum(2.* OptNEHT_Trap_dL(zs, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL) )
	# # Trap_int = 0.0
	# # for i in range(Ntrap_z):
	# # 	Trap_int += 2.* OptNEHT_Trap_dL(zs[i], Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL) 

	# # Trap_int -= Trap_sub
	# return 4.*ma.pi * (zmax-0.000001)/(2.*Ntrap_z) * Trap_int
	


def IntzZ_OptNEHT(p, zmax, Mmx, Fmin, chi, thMn, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL):
	print p
	eps, KQ = p
	#z, eps, KQ, MdEff = p
	#MdEff = 0.1
	#f_Edd = eps
	#NEHT = OptNEHT_dz(zmax, eps, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL)
	if (KQ<0.0 or eps<0.0):
	 	NEHT = -100000000.
	else:
		NEHT = OptNEHT_dz(zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL)

	print "NEHT=%g" %NEHT
	return -NEHT



def IntzZ_Trap_OptNEHT(p, zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL):
	print p
	eps, KQ = p
	#z, eps, KQ, MdEff = p
	#MdEff = 0.1
	if (KQ<0.0 or eps<0.0):
	 	NEHT = -100000000.
	else:
		NEHT = OptNEHT_Trap_dz(zmax, Mmx, eps, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL)
	print "NEHT=%g" %NEHT
	return -NEHT

























# def FtotDZ(DZ, z, h, Om, OL, Mmax, Mmin, thMn, qmin, eps, Pbase, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim):
# 	arg =  4.*np.pi *DZ* dVdzdOm(z, h, Om, OL) * fbin_GWgas(z, Mmax, Mmin, thMn, qmin, eps, Pbase, KQ, MdEff, xi, fbin, h, Om, OL) * FInt_smLF(z, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim, h, Om, OL) 	
# 	return arg/(10.**6 * pc2cm)**3  ##LF in Mpc^(-3)	





# def Ftot(zmax, h, Om, OL, alpha, beta, gamma, Mmax, Mmin, thMx, thMn, Pbase, PNyq, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim, tQuasar, xii):
# 	if (type(zmax) is float or type(zmax) is np.float64):
# 		return intg.quad(Fbin_Integrand, 0., zz,  args=(h, Om, OL, alpha, beta, gamma, Mmax, Mmin, thMx, thMn, Pbase, PNyq, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim, tQuasar, xii), epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
# 	else:
# 		for i in range(len(zmax)):
# 			res.append( intg.quad(Fbin_Integrand, 0., zz[i],  args=(h, Om, OL, alpha, beta, gamma, Mmax, Mmin, thMx, thMn, Pbase, PNyq, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim ,tQuasar, xii), epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
# 		return np.array(res)


# def Ftot_GWgas(zmax, h, Om, OL, Mmax, Mmin,thMn, qmin, eps, Pbase, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim, tQuasar, xii):
# 	if (type(zmax) is float or type(zmax) is np.float64):
# 		return intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=(h, Om, OL, Mmax, MMin, thMn, qmin, eps, Pbase, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
# 	else:
# 		for i in range(len(zmax)):
# 			res.append( intg.quad(Fbin_Integrand_GWgas, 0., zmax[i],  args=(h, Om, OL, Mmax, Mmin, thMn, qmin, eps, Pbase, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim), epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
# 		return np.array(res)



# def Ntot_Integrand(z, h, Om, OL, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim):
# 	return dVdzdOm(z, h, Om, OL) * FInt_smLF(z, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim, h, Om, OL) /(10.**6 * pc2cm)**3   ##LF in Mpc^(-3)	

# def Ntot(zmax, h, Om, OL, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim):
# 	if (type(zmax) is float or type(zmax) is np.float64):
# 		return 	intg.quad( Ntot_Integrand, 0, zmax, args=(h, Om, OL, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] 
# 	else:
# 		res=[]
# 		for i in range(len(zmax)):
# 			res.append(  intg.quad(Ntot_Integrand, 0, zmax[i], args=(h, Om, OL, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]  )
# 		return np.array(res)




# def NtotDZ(DZ, z, h, Om, OL, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim):
# 	return 4.*np.pi *DZ* dVdzdOm(z, h, Om, OL)/(10.**6 * pc2cm)**3  * FInt_smLF(z, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim, h, Om, OL)













