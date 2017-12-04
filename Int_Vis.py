import numpy as np


import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 18})
import matplotlib.pyplot as plt


import scipy.integrate as intg
import math as ma

import GWB_IntFuncs as hGWB 
from GWB_IntFuncs import *


## method optioncs
fEdd_Dist = True
textlabels = True

ztst=0.1
f_Edd = 0.1#10**(-4.1)  

if (fEdd_Dist):
	Ng = 100
else:
	Ng = 100
sumz = 200
NgHi = 200


if (fEdd_Dist):
	import VLBI_IntFuncs_fEddDist as IFs 
	from VLBI_IntFuncs_fEddDist import *
else:
	import VLBI_IntFuncs_V2 as IFs 
	from VLBI_IntFuncs_V2 import *


##Constants
c = 2.9979*10**(10)
G = 6.673*10**(-8)
Msun = 1.998*10**(33)
pc2cm = 3.08567758*10**(18)
yr2sec = 3600.*24.*365.25

mu_as2rad = ma.pi/180. / 3600. /10.**6
mJy2cgs = 10.**(-26)





## ALL THE PARAMETERS!
###MEASURED PARAMETERS















#Accretion Params
KQ = 10.**(-1.0) ## sets number of pc at which RL turns on
eps = 1.0#10**(-3.75)  ## sets migration (accretion rate in CBD pushing two together)

# if (Lmx==24.0):	
# 	f_Edd = 0.0001
# 	fbin = 1.00
# else:
# 	f_Edd = 0.1  ## sets L to M connection (accretion rate onto shining BH)

fbin = 1.0
MdEff = 0.1

## binary pop params
qmin_EHT = 0.01   ### qminof EHT sample
qmin_POP = np.minimum(qmin_EHT, 0.01)  ### qmin of all MBHBS 

#zeval = 0.5  #eval at this z
zmax = 3.0 ### integrateo out to zmax=5.0

##Instrument params
Fmin = 1.0 * mJy2cgs
thMn = 1.0 * mu_as2rad 
Pbase = 10.0*yr2sec


###Cosmology
h=0.7
Om = 0.3
OL=0.7


###
Mmx = 1.*10.**10 ## This is for LLAGN - choose L, and randomly draw f_edd can give very large inferred masses - so limit by obs knowledge of BH mass
Mmax = 1.*10.**10*Msun
Mmin= 0.0*10.**5*Msun ## no really necessary, but in case wan tot limit this


###SED scaling
chi  = 0.5 #Elvis 1994 Radio Loud AGN - nulnu propto nu^(0.9) and #0.1  #chi L_(408MHz) = L_mm (specifc fluxes)
xi=1.0


#JET CONSTANTS
nummGHz = c/0.1/1.e9 ##1mm in GHz
gamj = 10.0 ##typical, maybe high
thobs = 1./gamj#ma.pi/4. ##Eval at thobs such that nu_SSA is max, most conservative?
ke = 1.0 #constant order unity
Delc = np.log(1.e5) ##COul log of rmax/rmin
Lam = np.log(1.e5) ##COul log of gam_e max/gam_e min


#### Shade region that is bright enough
def Flxnu(z, Mbn, f_Edd, h, Om, OL, Fmin, fEdd_Dist):
	if (fEdd_Dist):
		Fmm = Mbn2Lmm(Mbn)*1.e7/4./np.pi/( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2
	else:
		Fmm = Mbn2Lmm(Mbn, f_Edd)*1.e7/4./np.pi/( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2
	if (Fmm>=Fmin):
		return 0.0
	else:
		return 1.0


def FlxnuL(z, Lmm, f_Edd, h, Om, OL, Fmin, fEdd_Dist):
	if (fEdd_Dist):
		Fmm = 10.**Lmm*1.e7/4./np.pi/( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2
	else:
		Fmm = 10.**Lmm*1.e7/4./np.pi/( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2
	if (Fmm>=Fmin):
		return 0.0
	else:
		return 1.0


def asep_bnd(P, M, z, thmn, h, Om, OL, Pbase, kq):

	# ##P passed is obs frame
	Npc = kq*pc2cm
	PMax = PmaxNPC(KQ*pc2cm, M)
	Pbase = np.minimum(Pbase, PMax*(1.+z))
	PMin = np.maximum( PminRes(M, thmn, z, h, Om, OL), PISCO(M) ) 

	#if (asep(P, M) >= thmn*Dang(z, h, Om, OL) and Pbase>PMin*(1.+z) and P*(1.+z)<=Pbase):
	if (P>=PMin*(1.+z) and Pbase>PMin*(1.+z) and P<=Pbase):
	#if (P>=PMin*(1.+z) and P<=Pbase):
		return 0.0
	else:
		return 1.0

def nubnd(P, M, z, fEdd, thobs, gamj, ke, Delc, Lam):
	P = P/(1.+z)
	if (nu_SSA(z, asep(P,M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)<=nummGHz and nu_loss(z, asep(P,M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)>=nummGHz):
		return 0.0
	else:
		return 1.0

#nubnd_L(10**Pbasez[j]*yr2sec, Lmm2Mbn_draw(Lmms[i])*Msun, ztst, Lmms[i], thobs, gamj, ke, Delc, Lam)
def nubnd_L(P, M, z, Lmm, thobs, gamj, ke, Delc, Lam):
	P = P/(1.+z)
	if (nu_SSA_L(z, asep(P,M)/pc2cm/2., Lmm, thobs, gamj, ke, Delc, Lam)<=nummGHz and nu_loss_L(z, asep(P,M)/pc2cm/2., Lmm, thobs, gamj, ke, Delc, Lam)>=nummGHz):
		return 0.0
	else:
		return 1.0


def FlxfrmL(Lmm, z, h, Om, OL):
	return Lmm*1.e7/4./ma.pi/( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2


def FlxfrmM(Mbn, f_Edd, z, h, Om, OL):
	BCXray = 50.0#4.2 ## Lbol = BC lambda L_lambda From 1+2011 at 145 nm Runnoe+2012 Table 2 https://arxiv.org/pdf/1201.5155v1.pdf 
	#nu14 = 1.4e9
	numm = c/(0.1)
	Lbol = Mbn*(f_Edd * LEdd_Fac ) #* Msun
	Lmm = Lbol/BCXray/numm
	return Lmm/4./ma.pi/( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2




def Flx2M(Flx, f_Edd, z, h, Om, OL):
	BCXray = 50.0#4.2 
	numm = c/(0.1)
	Lmm = (Flx * 4.*ma.pi*( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2)
	Lbol = BXray*Lmm*numm
	Mbn = Lbol  /(f_Edd * LEdd_Fac * Msun )
	return Mbn
	


PbasezHi = np.linspace(-1.0, np.log10(2.*Pbase/yr2sec), NgHi)
MbnzHi = np.linspace(5.0, 11.5, NgHi)

Pbasez = np.linspace(-1.0, np.log10(2.*Pbase/yr2sec), Ng)
Mbnz = np.linspace(5.0, 11.5, Ng)

Mavgs = np.zeros(Ng)
MavgsHi = np.zeros(NgHi)
Lmms = np.linspace(22.0, 26.0, Ng)
LmmsHi = np.linspace(22.0, 26.0, NgHi)
FlxsL = FlxfrmL(10.**Lmms, ztst, h, Om, OL)/mJy2cgs
FlxsLHi = FlxfrmL(10.**LmmsHi, ztst, h, Om, OL)/mJy2cgs
FlxsM = FlxfrmM(10.**Mbnz*Msun, f_Edd, ztst, h, Om, OL)/mJy2cgs
FlxsMHi = FlxfrmM(10.**MbnzHi*Msun, f_Edd, ztst, h, Om, OL)/mJy2cgs


Int_grid = np.zeros([Ng,Ng])
nubnds_avg = np.zeros([Ng,Ng])

FlxmnsLo = np.zeros([Ng,Ng])
abndsLo = np.zeros([Ng,Ng])
nubndsLo = np.zeros([Ng,Ng])

Flxmns = np.zeros([NgHi,NgHi])
abnds = np.zeros([NgHi,NgHi])
nubnds = np.zeros([NgHi,NgHi])
nubnds_avgHi = np.zeros([NgHi,NgHi])





#Int_grid_avg = np.zeros([Ng,Ng])


#### Contour the integrand
###ContourPlot DN/DL * F
#F depends on M, but integreates over Ps, so can only vary Pbase
if (fEdd_Dist):
	for k in range(sumz):
		for i in range(0,Ng):
			Mbndraw = Lmm2Mbn_draw(Lmms[i])
			for j in range(0,Ng):
				#Int_grid[j][i] =  Int_grid[j][i] +  Fbin_Integrand_GWgas(np.log10(Mbn2Lmm_fxd(10**Mbnz[i]*Msun, 0.001)), ztst, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, 10**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL)
				Int_grid[j][i]   = Int_grid[j][i]   + Fbin_Integrand_GWgas(Lmms[i], ztst, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, 10**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL)
				nubnds_avg[j][i] = nubnds_avg[j][i] + nubnd_L(10**Pbasez[j]*yr2sec, Mbndraw*Msun, ztst, Lmms[i], thobs, gamj, ke, Delc, Lam)
				#Flxmns[j][i]   =  Flxmns[j][i] + Flxnu(ztst, 10**Mbnz[i]*Msun, f_Edd, h, Om, OL, Fmin,fEdd_Dist)
			Mavgs[i] = Mavgs[i] + Mbndraw * np.exp(-(Lmm2Mbn_draw(Lmms[i])/Mmx)**4)
		print "step %g/%g" %(k+1, sumz)
	Int_grid = Int_grid/sumz
	#Flxmns = Flxmns/sumz
	Int_grid = 4.*np.pi*ztst * Int_grid * np.log10(Mbn2Lmm(10**Mbnz[Ng/2]*Msun))
	Mavgs = np.log10(Mavgs/sumz)
	nubnds_avg = nubnds_avg/sumz
	#Mavgs = Mavgs/sumz
else:		
	for i in range(0,Ng):
		for j in range(0,Ng):
			Int_grid[j][i] = Fbin_Integrand_GWgas(np.log10(Mbn2Lmm(10**Mbnz[i]*Msun, f_Edd)), ztst, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, 10**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL)
			
			

	Int_grid = np.maximum(1.e-14, np.nan_to_num( 4.*np.pi*ztst * Int_grid * np.log10(Mbn2Lmm(10**Mbnz[Ng/2]*Msun, f_Edd)) ) )
 


if (fEdd_Dist):
	for i in range(0, Ng):
		for j in range(0,Ng):
			FlxmnsLo[j][i]  = FlxnuL(ztst, Lmms[i], f_Edd, h, Om, OL, Fmin,fEdd_Dist)
			abndsLo[j][i]   = asep_bnd(10**Pbasez[j]*yr2sec, 10**Mbnz[i]*Msun, ztst, thMn, h, Om, OL, Pbase, KQ)
			nubndsLo[j][i]  = nubnd(10**Pbasez[j]*yr2sec, 10**Mbnz[i]*Msun, ztst, f_Edd, thobs, gamj, ke, Delc, Lam)

	# for k in range(sumz):
	# 	for i in range(0, NgHi):
	# 		for j in range(0,NgHi):
	# 			nubnds_avgHi[j][i]   = nubnd_L(10**PbasezHi[j]*yr2sec, Lmm2Mbn_draw(LmmsHi[i])*Msun, ztst, LmmsHi[i], thobs, gamj, ke, Delc, Lam)
	# 		#MavgsHi[i] = MavgsHi[i] + Lmm2Mbn_draw(LmmsHi[i])
	# nubnds_avg = nubnds_avg/sumz
	# #MavgsHi = np.log10(MavgsHi/sumz)
else:	
	for i in range(0, NgHi):
		for j in range(0,NgHi):
			Flxmns[j][i]   = Flxnu(ztst, 10**MbnzHi[i]*Msun, f_Edd, h, Om, OL, Fmin,fEdd_Dist)
			abnds[j][i]    = asep_bnd(10**PbasezHi[j]*yr2sec, 10**MbnzHi[i]*Msun, ztst, thMn, h, Om, OL, Pbase, KQ)
			nubnds[j][i]   = nubnd(10**PbasezHi[j]*yr2sec, 10**MbnzHi[i]*Msun, ztst, f_Edd, thobs, gamj, ke, Delc, Lam)


# abnds   = asep_bnd(10**Pbasez*yr2sec, 10**Mbnz*Msun, ztst, thMn, h, Om, OL, Pbase)
# nubnds  = nubnd(10**Pbasez*yr2sec, 10**Mbnz*Msun, ztst, f_Edd, thobs, gamj, ke, Delc, Lam)






	

fig=plt.figure(figsize=[7.5,6.6])
ax = fig.add_subplot(111)



# if (fEdd_Dist):
# 	plt.title(r"$\frac{dN^2}{dL dz}\Delta L\Delta z \mathcal{F}$, z=%g$" %(ztst))
# else:
# 	plt.title(r"$\frac{dN^2}{dL dz}\Delta L\Delta z \mathcal{F}$, z=%g, $f_{\rm{Edd}} = 10^{%g}}$" %(ztst,np.log10(f_Edd)), fontsize=14)
# plt.title(r"$\frac{d^2N}{dL dz} L z \mathcal{F}$", fontsize=15)

# integrand contours
#cnt = plt.contourf(Mbnz, Pbasez, np.log10(Int_grid), 200, cmap = "viridis", zorder=0)





if (fEdd_Dist):
	ax.contourf(np.log10(FlxsL), Pbasez, np.log10(Int_grid), 200, cmap = "viridis", zorder=0)
	ax2 = ax.twiny()

	
	ax2.contourf(Mavgs, Pbasez, np.log10(Int_grid), 200, cmap = "viridis", zorder=0)
#
	#cbar = plt.colorbar(cnt)
	lms = ax.contour(np.log10(FlxsL), Pbasez, np.log10(Int_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],zorder=5)
	ax.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	
#
	lms = ax2.contour(Mavgs, Pbasez, np.log10(Int_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],zorder=5)
	ax2.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	
else:
	ax.contourf(np.log10(FlxsM), Pbasez, np.log10(Int_grid), 200, cmap = "viridis", zorder=5)
	ax2=ax.twiny()
#
#
	ax2.contourf(Mbnz, Pbasez, np.log10(Int_grid), 200, zorder=-10)
	#cbar = plt.colorbar(cnt)
	lms = ax.contour(np.log10(FlxsM), Pbasez, np.log10(Int_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],zorder=5)
	ax.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	
#
	lms2 = ax2.contour(Mbnz, Pbasez, np.log10(Int_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],zorder=5)
	ax2.clabel(lms2, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	



# ax2 = fig.add_subplot(111)
# ax2.contourf(np.log10(FlxsL), Pbasez, np.log10(Int_grid), 200, cmap = "viridis", zorder=0)
# lms = ax2.contour(np.log10(FlxsL), Pbasez, np.log10(Int_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],zorder=5)
# ax2.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	



##asep bounds
# ax2.contourf(MbnzHi, PbasezHi, np.log10(abnds), colors="red", alpha=0.3, zorder=10)
# ax2.contour(MbnzHi, PbasezHi, abnds, colors="red", linewidths=[3], linestyle="--", levels = [0.5], zorder=10)




##Fmin bounds
if (fEdd_Dist):
	# ax.contourf(np.log10(FlxsLHi), PbasezHi, np.log10(Flxmns-0.5), colors="black", alpha=0.3, zorder=10)
	# ax.contour(np.log10(FlxsLHi), PbasezHi, Flxmns, colors="black", linewidths=[3], levels = [0.5], zorder=10)
	# ##nu bounds
	# ax.contourf(np.log10(FlxsL), Pbasez, np.log10(nubnds_avg-0.5), colors="yellow", alpha=0.3, zorder=10)
	# ax.contour(np.log10(FlxsL), Pbasez, nubnds_avg, colors="yellow", linewidths=[3], linestyle=":", levels = [0.5], zorder=10)
	# ##asep bounds
	# ax.contourf(np.log10(FlxsLHi), PbasezHi, np.log10(abnds), colors="red", alpha=0.3, zorder=10)
	# ax.contour(np.log10(FlxsLHi), PbasezHi, abnds, colors="red", linewidths=[3], linestyle="--", levels = [0.5], zorder=10)

	ax2.contourf(Mavgs, Pbasez, np.log10(FlxmnsLo-0.5), colors="black", alpha=0.3, zorder=10)
	ax2.contour(Mavgs, Pbasez, FlxmnsLo, colors="black", linewidths=[3], levels = [0.5], zorder=10)
#
	##nu bounds
	ax2.contourf(Mavgs, Pbasez, np.log10(nubnds_avg-0.5), colors="yellow", alpha=0.3, zorder=10)
	ax2.contour(Mavgs, Pbasez, nubnds_avg, colors="yellow", linewidths=[3], linestyle=":", levels = [0.5], zorder=10)
#
	##asep bounds
	ax2.contourf(Mavgs, Pbasez, np.log10(abndsLo), colors="red", alpha=0.3, zorder=10)
	ax2.contour(Mavgs, Pbasez, abndsLo, colors="red", linewidths=[3], linestyle="--", levels = [0.5], zorder=10)
else:
	# ax.contourf(np.log10(FlxsMHi), PbasezHi, np.log10(Flxmns), colors="black", alpha=0.3, zorder=10)
	# ax.contour(np.log10(FlxsMHi), PbasezHi, Flxmns, colors="black", linewidths=[3], levels = [0.5], zorder=10)
	# ##nu bounds
	# ax.contourf(np.log10(FlxsMHi), PbasezHi, np.log10(nubnds), colors="yellow", alpha=0.3, zorder=10)
	# ax.contour(np.log10(FlxsMHi), PbasezHi, nubnds, colors="yellow", linewidths=[3], linestyle=":", levels = [0.5], zorder=10)
#
	# ##asep bounds
	# ax.contourf(np.log10(FlxsMHi), PbasezHi, np.log10(abnds), colors="red", alpha=0.3, zorder=10)
	# ax.contour(np.log10(FlxsMHi), PbasezHi, abnds, colors="red", linewidths=[3], linestyle="--", levels = [0.5], zorder=10)
#	
	ax2.contourf(MbnzHi,  PbasezHi, np.log10(Flxmns), colors="black", alpha=0.3, zorder=10)
	ax2.contour(MbnzHi,  PbasezHi, Flxmns, colors="black", linewidths=[3], levels = [0.5], zorder=10)
	##nu bounds
	ax2.contourf(MbnzHi, PbasezHi, np.log10(nubnds), colors="yellow", alpha=0.3, zorder=10)
	ax2.contour(MbnzHi, PbasezHi, nubnds, colors="yellow", linewidths=[3], linestyle=":", levels = [0.5], zorder=10)
	##asep bounds
	ax2.contourf(MbnzHi, PbasezHi, np.log10(abnds), colors="red", alpha=0.3, zorder=10)
	ax2.contour(MbnzHi, PbasezHi, abnds, colors="red", linewidths=[3], linestyle="--", levels = [0.5], zorder=10)


# ax2.contourf(MbnzHi, PbasezHi, np.log10(nubnds), colors="yellow", alpha=0.3, zorder=10)
# ax2.contour(MbnzHi, PbasezHi, nubnds, colors="yellow", linewidths=[3], linestyle=":", levels = [0.5], zorder=10)



#Diags
#plt.plot(Mbnz, np.log10(PminRes(10.**Mbnz*Msun, thMn, ztst, h, Om, OL)/yr2sec ) )
#plt.plot(Mbnz, np.log10(PISCO(10.**Mbnz*Msun)/yr2sec ) )


FminSv = Fmin/mJy2cgs/1000.
thMnSv = thMn/mu_as2rad 
PbaseSv = Pbase/yr2sec
Lmx_cgs = Lmx + 7.0

if (textlabels):
	if (fEdd_Dist):
		plt.figtext(0.19,0.47, r"$z=%g$" %ztst, color='yellow', fontsize=15)
	else:
		plt.figtext(0.19,0.52, r"$f_{\rm{Edd}}=10^{%g}$" %np.log10(f_Edd), color='yellow', fontsize=15)
		plt.figtext(0.19,0.47, r"$z=%g$" %ztst, color='yellow', fontsize=15)

	plt.figtext(0.19,0.42, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
	plt.figtext(0.19,0.37, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
	if (FminSv*1000.0 <= 100):
		FminSv = FminSv*1000.0
		plt.figtext(0.19,0.32, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
	else:
		plt.figtext(0.19,0.32, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
	plt.figtext(0.19,0.27, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
	plt.figtext(0.19,0.22, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ), color='yellow', fontsize=15)
	plt.figtext(0.19,0.17, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)


if (fEdd_Dist):
	ax.set_ylabel(r'$\rm{log}_{10}[P_{\rm{base}}/yr]$')
	ax2.set_xlabel(r'$\rm{log}_{10}\left<M/M_{\odot}\right>$')
	ax.set_xlabel(r'$\rm{log}_{10}[F_{mm}/\rm{mJy}]$')
	#ax2.set_xlim(Mavgs[0], Mavgs[Ng-2])
else:
	ax.set_ylabel(r'$\rm{log}_{10}[P_{\rm{base}}/yr]$')
	ax2.set_xlabel(r'$\rm{log}_{10}[M/M_{\odot}]$')
	ax.set_xlabel(r'$\rm{log}_{10}[F_{mm}/\rm{mJy}]$')



plt.tight_layout()

Savename = "diffN_fEddDist%g_zeval%g_fEdd%g_amx%g_avg%g.png" %(fEdd_Dist, ztst, f_Edd, KQ, sumz)

Savename = Savename.replace('.', 'p')
Savename = Savename.replace('ppng', '.png')
plt.savefig(Savename)


#plt.show()

