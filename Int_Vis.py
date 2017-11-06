import numpy as np


import matplotlib
#matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt


import scipy.integrate as intg
import math as ma

import GWB_IntFuncs as hGWB 
from GWB_IntFuncs import *


## method optioncs
fEdd_Dist = True


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
f_Edd = 10**(-3.5)  
MdEff = 0.1

## binary pop params
qmin_EHT = 0.01   ### qminof EHT sample
qmin_POP = np.minimum(qmin_EHT, 0.01)  ### qmin of all MBHBS 

zeval = 0.5  #eval at this z
zmax = 5.0 ### integrateo out to zmax=5.0

##Instrument params
Fmin = 10.0 * mJy2cgs
thMn = 1.0 * mu_as2rad 
Pbase = 10.0*yr2sec


###Cosmology
h=0.7
Om = 0.3
OL=0.7


###DEFUNCT
Mmx = 2000000.*10.**10 ## jsut to not limit lum function - doesnt change anyhting when set at 2*10^10 number
Mmax = 2000000.*10.**10*Msun
Mmin= 0.0*10.**5*Msun 

###SED scaling
chi  = 0.5 #Elvis 1994 Radio Loud AGN - nulnu propto nu^(0.9) and #0.1  #chi L_(408MHz) = L_mm (specifc fluxes)
xi=1.0


#JET CONSTANTS
nummGHz = c/0.1/1.e9 ##1mm in GHz
thobs = 0.1 ##Eval at thobs such that nu_SSA is max, most conservative?
gamj = 10. ##typical, maybe high
ke = 1.0 #constant order unity
Delc = np.log(1.e5) ##COul log of rmax/rmin
Lam = np.log(1.e5) ##COul log of gam_e max/gam_e min



#### Shade region that is bright enough
def Flxnu(z, Mbn, f_Edd, h, Om, OL, Fmin):
	Fmm = Mbn2Lmm(Mbn, f_Edd)/4./np.pi/( (1.+ztst)**2 * Dang(z, h, Om, OL) )**2
	if (Fmm>=Fmin):
		return 0.0
	else:
		return 1.0


def asep_bnd(P, M, z, thmn, h, Om, OL, Pbase):
	if (asep(P, M) >= thmn*Dang(z, h, Om, OL)and P<Pbase):
		return 0.0
	else:
		return 1.0

def nubnd(P, M, z, fEdd, thobs, gamj, ke, Delc, Lam):
	if (nu_SSA(z, asep(P,M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)<=nummGHz and nu_loss(z, asep(P,M)/pc2cm, fEdd, M, thobs, gamj, ke, Delc, Lam)>=nummGHz):
		return 0.0
	else:
		return 1.0


ztst=0.2

if (fEdd_Dist):
	Ng = 20
else:
	Ng = 100
sumz = 200
Pbasez = np.linspace(-1.0, np.log10(2.*Pbase/yr2sec), Ng)
Mbnz = np.linspace(5.0, 11.5, Ng)


Int_grid = np.zeros([Ng,Ng])
Flxmns = np.zeros([Ng,Ng])
abnds = np.zeros([Ng,Ng])
nubnds = np.zeros([Ng,Ng])

#Int_grid_avg = np.zeros([Ng,Ng])


#### Contour the integrand
###ContourPlot DN/DL * F
#F depends on M, but integreates over Ps, so can only vary Pbase
if (fEdd_Dist):
	for k in range(sumz):
		for i in range(0,Ng):
			for j in range(0,Ng):
				Int_grid[j][i] =  Int_grid[j][i] + Fbin_Integrand_GWgas(np.log10(Mbn2Lmm(10**Mbnz[i]*Msun))-7., ztst, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, 10**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL)

	Int_grid = Int_grid/sumz
else:		
	for i in range(0,Ng):
		for j in range(0,Ng):
			Int_grid[j][i] = Fbin_Integrand_GWgas(np.log10(Mbn2Lmm(10**Mbnz[i]*Msun, f_Edd))-7., ztst, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps, f_Edd, 10**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL)
			
Int_grid = 4.*np.pi*ztst*np.log10(1.e44/(c/0.1)) * Int_grid


for i in range(0,Ng):
		for j in range(0,Ng):
			Flxmns[j][i]   = Flxnu(ztst, 10**Mbnz[i]*Msun, f_Edd, h, Om, OL, Fmin)
			abnds[j][i]    = asep_bnd(10**Pbasez[j]*yr2sec, 10**Mbnz[i]*Msun, ztst, thMn, h, Om, OL, Pbase)
			nubnds[j][i]   = nubnd(10**Pbasez[j]*yr2sec, 10**Mbnz[i]*Msun, ztst, f_Edd, thobs, gamj, ke, Delc, Lam)


### Shade resolvable separations



### shade regio



	

plt.figure(figsize=[7.5,6.1])

# if (fEdd_Dist):
# 	plt.title(r"$\frac{dN^2}{dL dz}\Delta L\Delta z \mathcal{F}$, z=%g$" %(ztst))
# else:
# 	plt.title(r"$\frac{dN^2}{dL dz}\Delta L\Delta z \mathcal{F}$, z=%g, $f_{\rm{Edd}} = 10^{%g}}$" %(ztst,np.log10(f_Edd)), fontsize=14)
plt.title(r"$\frac{dN^2}{dL dz}\Delta L\Delta z \mathcal{F}$")

### integrand contours
cnt = plt.contourf(Mbnz, Pbasez, np.log10(Int_grid), 200, cmap = "viridis", zorder=0)

##Fmin bounds
plt.contourf(Mbnz, Pbasez, np.log10(Flxmns), colors="black", alpha=0.3, zorder=10)
plt.contour(Mbnz, Pbasez, Flxmns, colors="black", linewidths=[3], levels = [0.5], zorder=10)

##asep bounds
plt.contourf(Mbnz, Pbasez, np.log10(abnds), colors="red", alpha=0.3, zorder=10)
plt.contour(Mbnz, Pbasez, abnds, colors="red", linewidths=[3], linestyle="--", levels = [0.5], zorder=10)


##nu bounds
plt.contourf(Mbnz, Pbasez, np.log10(nubnds), colors="yellow", alpha=0.3, zorder=10)
plt.contour(Mbnz, Pbasez, nubnds, colors="yellow", linewidths=[3], linestyle=":", levels = [0.5], zorder=10)


#cbar = plt.colorbar(cnt)
lms = plt.contour(Mbnz, Pbasez, np.log10(Int_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	



FminSv = Fmin/mJy2cgs/1000.
thMnSv = thMn/mu_as2rad 
PbaseSv = Pbase/yr2sec
Lmx_cgs = Lmx + 7.0

if (fEdd_Dist):
	plt.figtext(0.17,0.47, r"$z=%g$" %ztst, color='yellow', fontsize=14)
else:
	plt.figtext(0.17,0.52, r"$f_{\rm{Edd}}=10^{%g}$" %np.log10(f_Edd), color='yellow', fontsize=14)
	plt.figtext(0.17,0.47, r"$z=%g$" %ztst, color='yellow', fontsize=14)

plt.figtext(0.17,0.42, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=14)
plt.figtext(0.17,0.37, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=14)
plt.figtext(0.17,0.32, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=14)
plt.figtext(0.17,0.27, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=14)
plt.figtext(0.17,0.22, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ), color='yellow', fontsize=14)
plt.figtext(0.17,0.17, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=14)

plt.ylabel(r'$\rm{log}_{10}[P_{\rm{base}}/yr]$')
plt.xlabel(r'$\rm{log}_{10}[M/M_{\odot}]$')

plt.tight_layout()

Savename = "diffN_fEddDist%g_zeval%g_fEdd%g.png" %(fEdd_Dist, ztst, f_Edd)

Savename = Savename.replace('.', 'p')
Savename = Savename.replace('ppng', '.png')
plt.savefig(Savename)


#plt.show()

