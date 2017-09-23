import numpy as np


import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt

import scipy as sc
import scipy.integrate as intg
from scipy import optimize
from scipy.optimize import fmin

import math as ma
# import VLBI_IntFuncs_V2 as IFs 
# from VLBI_IntFuncs_V2 import *


## method optioncs
fEdd_Dist = True

if (fEdd_Dist):
	import VLBI_IntFuncs_fEddDist as IFs 
	from VLBI_IntFuncs_fEddDist import *
else:
	import VLBI_IntFuncs_V2 as IFs 
	from VLBI_IntFuncs_V2 import *



#OPTIONS
Fit = False
logN = True
fit_MdEff = False
TrapInt = True


##Constants
c = 2.9979*10**(10)
G = 6.673*10**(-8)
Msun = 1.998*10**(33)
pc2cm = 3.08567758*10**(18)
yr2sec = 3600.*24.*365.25

mu_as2rad = ma.pi/180. / 3600. /10.**6
mJy2cgs = 10.**(-26) ###1Jy is 10^(-23) erg/s/cm^2/Hz 










##IMPORTANT PARAMS HERE
fbin = 1.0
Fmin = 10.0 * mJy2cgs
qmin_EHT = 0.01   ### qminof EHT sample



Pbase = 10.0*yr2sec
KQ = 10.**(-1.0)  #0.01 ## sets number of pc at which RL turns on
qmin_POP = np.minimum(qmin_EHT, 0.01)  ### qmin of all MBHBS 

eps = 1.0  ## sets migration (accretion rate in CBD pushing two together)
f_Edd = 0.1 

# if (Lmx==24.0):
# 	f_Edd = 0.0001
# 	fbin = 0.001
# else:
# 	f_Edd = 1.0  ## sets L to M connection (accretion rate onto shining BH)
# 	fbin = 0.1













## ALL THE PARAMETERS!
###MEASURED PARAMETERS
h=0.7
Om = 0.3
OL=0.7


Mmx = 20000.*10.**10 ## jsut to not limit lum function - doesnt change anyhting when set at 2*10^10 number
Mmax = 20000.*10.**10*Msun
Mmin= 0.**5*Msun 





zmax = 10.0

##FREE (and KEY) PARAMETERS
chi  = 0.5 #From power law estimate of Elvis 1994  Radio Loud AGN (3e11/0.408e9)^(0.9-1.0)
#chi = 1./0.12  #chi L_14 = L_mm (specific Lums) 1.0 for blazars, 1./0.12 for Sgr A*
# for Sgr A*, nuLnu \propto nu^(1.4) (Loeb Waxman 2007) 


MdEff = 0.1
xi = 1.0
zeval = 0.5
DZ = 1.0
thMn = 1.0 * mu_as2rad 







#tsts = 400
#NEHT_tst = np.zeros(tsts)
#p0 = [zeval, eps, KQ, MdEff]

p0 = [zeval, eps, KQ]

#if not fitting:
pnoFit = [eps, KQ]

thmin_max = 20.0
Nz = 5
Ng = int(1.*thmin_max)/2
popt = np.zeros([Ng,len(p0)])
NEHT =np.zeros(Ng)
NEHT_Z1 =np.zeros(Ng)
NEHT_Z2 =np.zeros(Ng)
NEHT_Z3 =np.zeros(Ng)
NEHT_Z4 =np.zeros(Ng)
NEHT_Z5 =np.zeros(Ng)
NEHT_Z6 =np.zeros(Ng)
NEHT_Z7 =np.zeros(Ng)

thMns = np.linspace(1.0, thmin_max, Ng) 
#Zs = [0.02, 0.05, 0.1, 0.5, 1.0]#, 2.0]
Zs = [0.01, 0.1, 1.0, 2.0, 3.0]
#Zs = [0.01, 0.05, 0.1, 1.0, 1.5, 2.0, 2.5]



if (Fit):
	for i in range(0, Ng):
		print "i=%g" %i
		print r"$\theta_{\min} = %g$" %thMns[i]


		#initial guess
		#for j in range(tsts)
		NEHT_tst = 0.0
		# NEHT_tstZ1 = 0.0
		# NEHT_tstZ2 = 0.0
		# NEHT_tstZ3 = 0.0
		j=0
		k=0
		Nlim = 1.0
		print "Generating ICs"
		while ( (NEHT_tst<=Nlim) or j>10.**(len(p0))):
		#while ( (NEHT_tstZ1<=Nlim and NEHT_tstZ2<=Nlim and NEHT_tstZ3<=Nlim) or j>10.**(len(p0))):
			zeval = np.random.rand()
			eps = 10.**( -4.0*np.random.rand() )
			KQ = 10.**( -3.0*np.random.rand() )
			#MdEff = 10.**( -2.0*np.random.rand() )
			#p0 = [zeval, eps, KQ, MdEff]
			p0 = [zeval, eps, KQ]
			#p01 = [eps, KQ]
			
			# eps = 10.**( -4.0*np.random.rand() )
			# KQ = 10.**( -3.0*np.random.rand() )
			# p02 = [eps, KQ]

			# eps = 10.**( -4.0*np.random.rand() )
			# KQ = 10.**( -3.0*np.random.rand() )
			# p03 = [eps, KQ]
			NEHT_tst = -OptNEHT(p0, Mmx, DZ, Fmin, chi, thMns[i]* mu_as2rad , qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
			#NEHT_tstZ1 = -IntzZ_OptNEHT(p01, Zs[0], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL)
			#NEHT_tstZ2 = -IntzZ_OptNEHT(p02, Zs[1], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL)
			#NEHT_tstZ3 = -IntzZ_OptNEHT(p03, Zs[2], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL)

			j = j+1
			if (k==0 and j>10.**(len(p0))):
				k=1
				j=0
				Nlim = 0.0  ##try again with a weaker constraint
			# if (NEHT_tst>0.0):
			# 	pmax = p0
			# if (j==0):
			# 	pmax = p0
			# elif(  NEHT_tst[j]  >=  np.max(NEHT_tst) ):
			# 	pmax = p0
		
		# print "Optimizing Z1"
		# poptZ1[i] = sc.optimize.fmin(IntzZ_OptNEHT, p01, args=(Zs[0], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL), full_output=1, disp=False, ftol=1.e-2)[0]
		# NEHT_Z1 = -IntzZ_OptNEHT(poptZ1[i], Zs[0], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL)
		
		# print "Optimizing Z2"
		# poptZ2[i] = sc.optimize.fmin(IntzZ_OptNEHT, p02, args=(Zs[1], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL), full_output=1, disp=False, ftol=1.e-2)[0]
		# NEHT_Z2 = -IntzZ_OptNEHT(poptZ2[i], Zs[1], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL)
		
		# print "Optimizing Z3"
		# poptZ3[i] = sc.optimize.fmin(IntzZ_OptNEHT, p03, args=(Zs[2], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL), full_output=1, disp=False, ftol=1.e-2)[0]
		# NEHT_Z3 = -IntzZ_OptNEHT(poptZ3[i], Zs[2], Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, xi, fbin, h, Om, OL)

		popt[i]  = sc.optimize.fmin(OptNEHT,  p0, args=(Mmx, DZ, Fmin, chi, thMns[i]* mu_as2rad , qmin, Pbase, xi, fbin, h, Om, OL), full_output=1, disp=False, ftol=1.e-8)[0]
		NEHT[i] = -OptNEHT(popt[i], Mmx, DZ, Fmin, chi, thMns[i]* mu_as2rad , qmin, Pbase, xi, fbin, h, Om, OL)

		pdz = [popt[i][1], popt[i][2]]
		NEHT_Z1[i] = -IntzZ_OptNEHT(pdz, Zs[0], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
		NEHT_Z2[i] = -IntzZ_OptNEHT(pdz, Zs[1], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
		NEHT_Z3[i] = -IntzZ_OptNEHT(pdz, Zs[2], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
else:
	if (TrapInt):
		for i in range(0, Ng):
			#OptNEHT_dz(zmax, Mmx, eps, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL)
		

			NEHT_Z1[i] = -IntzZ_Trap_OptNEHT(pnoFit, Zs[0], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL)
		
			NEHT_Z2[i] = -IntzZ_Trap_OptNEHT(pnoFit, Zs[1], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL)
		
			NEHT_Z3[i] = -IntzZ_Trap_OptNEHT(pnoFit, Zs[2], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL)
			
			NEHT_Z4[i] = -IntzZ_Trap_OptNEHT(pnoFit, Zs[3], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL)

			NEHT_Z5[i] = -IntzZ_Trap_OptNEHT(pnoFit, Zs[4], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL)
			
			# NEHT_Z6[i] = -IntzZ_Trap_OptNEHT(pnoFit, Zs[5], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL)

			# NEHT_Z7[i] = -IntzZ_Trap_OptNEHT(pnoFit, Zs[6], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL)


		Ntot_Z1 = Ntot_Trap_RLF(Zs[0], Fmin, chi, h, Om, OL)
		Ntot_Z2 = Ntot_Trap_RLF(Zs[1], Fmin, chi, h, Om, OL)
		Ntot_Z3 = Ntot_Trap_RLF(Zs[2], Fmin, chi, h, Om, OL)
		Ntot_Z4 = Ntot_Trap_RLF(Zs[3], Fmin, chi, h, Om, OL)
		Ntot_Z5 = Ntot_Trap_RLF(Zs[4], Fmin, chi, h, Om, OL)
		# Ntot_Z6 = Ntot_Trap_RLF(Zs[5], Fmin, chi, h, Om, OL)
		# Ntot_Z7 = Ntot_Trap_RLF(Zs[6], Fmin, chi, h, Om, OL)

	else:
		for i in range(0, Ng):
			#OptNEHT_dz(zmax, Mmx, eps, KQ, Fmin, chi, thMn, qmin, Pbase, xi, fbin, h, Om, OL)
		

			NEHT_Z1[i] = -IntzZ_OptNEHT(pnoFit, Zs[0], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
		
			NEHT_Z2[i] = -IntzZ_OptNEHT(pnoFit, Zs[1], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
		
			NEHT_Z3[i] = -IntzZ_OptNEHT(pnoFit, Zs[2], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
			
			NEHT_Z4[i] = -IntzZ_OptNEHT(pnoFit, Zs[3], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)

			NEHT_Z5[i] = -IntzZ_OptNEHT(pnoFit, Zs[4], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
			
			# NEHT_Z6[i] = -IntzZ_OptNEHT(pnoFit, Zs[5], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
			
			# NEHT_Z7[i] = -IntzZ_OptNEHT(pnoFit, Zs[6], Mmx, Fmin, chi, thMns[i]* mu_as2rad, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)



		Ntot_Z1 = Ntot_RLF(Zs[0], Fmin, chi, h, Om, OL)
		Ntot_Z2 = Ntot_RLF(Zs[1], Fmin, chi, h, Om, OL)
		Ntot_Z3 = Ntot_RLF(Zs[2], Fmin, chi, h, Om, OL)
		Ntot_Z4 = Ntot_RLF(Zs[3], Fmin, chi, h, Om, OL)
		Ntot_Z5 = Ntot_RLF(Zs[4], Fmin, chi, h, Om, OL)
		# Ntot_Z6 = Ntot_RLF(Zs[5], Fmin, chi, h, Om, OL)
		# Ntot_Z7 = Ntot_RLF(Zs[6], Fmin, chi, h, Om, OL)











opac = 0.6
plt.figure(figsize=[8,6])
if (logN):
	# plt.scatter(thMns, np.log10(NEHT_Z1+1.0), color='#1b9e77')
	# plt.plot(thMns, np.log10(NEHT_Z1+1.0), color='#1b9e77', linewidth=3)

	# plt.scatter(thMns, np.log10(NEHT_Z2+1.0), color='#d95f02')
	# plt.plot(thMns, np.log10(NEHT_Z2+1.0), color='#d95f02', linewidth=3)

	# plt.scatter(thMns, np.log10(NEHT_Z3+1.0), color='#7570b3')
	# plt.plot(thMns, np.log10(NEHT_Z3+1.0), color='#7570b3', linewidth=3)

	# plt.scatter(thMns, np.log10(NEHT_Z4+1.0), color='#e7298a')
	# plt.plot(thMns, np.log10(NEHT_Z4+1.0), color='#e7298a', linewidth=3)

	# plt.scatter(thMns, np.log10(NEHT_Z5+1.0), color='gray')
	# plt.plot(thMns, np.log10(NEHT_Z5+1.0), color='gray', linewidth=3, linestyle='--')



	plt.scatter(thMns, np.log10(NEHT_Z1), color='#1b9e77')
	p1 = plt.plot(thMns, np.log10(NEHT_Z1), color='#1b9e77', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, np.log10(NEHT_Z2), color='#d95f02')
	p2 = plt.plot(thMns, np.log10(NEHT_Z2), color='#d95f02', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, np.log10(NEHT_Z3), color='#7570b3')
	p3 = plt.plot(thMns, np.log10(NEHT_Z3), color='#7570b3', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, np.log10(NEHT_Z4), color='#e7298a')
	p4 = plt.plot(thMns, np.log10(NEHT_Z4), color='#e7298a', linewidth=3, linestyle='--',alpha=opac)

	plt.scatter(thMns, np.log10(NEHT_Z5), color='#66a61e')
	p5 = plt.plot(thMns, np.log10(NEHT_Z5), color='#66a61e', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, np.log10(NEHT_Z6), color='#e6ab02')
	p6 = plt.plot(thMns, np.log10(NEHT_Z6), color='#e6ab02', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, np.log10(NEHT_Z7), color='#a6761d')
	p7 = plt.plot(thMns, np.log10(NEHT_Z7), color='#a6761d', linewidth=3, linestyle='--', alpha=opac)

	plt.ylim(-3.0, 3.0)






	#plt.scatter(thMns, np.log10(NEHT_Z6+1.0), color='black')
	#plt.plot(thMns, np.log10(NEHT_Z6+1.0), color='black', linewidth=3)

else:
	plt.scatter(thMns, NEHT_Z1, color='#1b9e77')
	p1 = plt.plot(thMns, NEHT_Z1, color='#1b9e77', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, NEHT_Z2, color='#d95f02')
	p2 = plt.plot(thMns, NEHT_Z2, color='#d95f02', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, NEHT_Z3, color='#7570b3')
	p3 = plt.plot(thMns, NEHT_Z3, color='#7570b3', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, NEHT_Z4, color='#e7298a')
	p4 = plt.plot(thMns, NEHT_Z4, color='#e7298a', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, NEHT_Z5, color='#66a61e')
	p5 = plt.plot(thMns, NEHT_Z5, color='#66a61e', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, NEHT_Z6, color='#e6ab02')
	p6 = plt.plot(thMns, NEHT_Z6, color='#e6ab02', linewidth=3, linestyle='--', alpha=opac)

	plt.scatter(thMns, NEHT_Z7, color='#a6761d')
	p7 = plt.plot(thMns, NEHT_Z7, color='#a6761d', linewidth=3, linestyle='--', alpha=opac)

	#plt.ylim(-3.0, 3.0)

plt.axhline(0.0, linestyle=":", color="black")

FminSv = Fmin/mJy2cgs/1000.
thMnSv = thMn/mu_as2rad 
PbaseSv = Pbase/yr2sec
Lmx_cgs = Lmx + 7.0

plt.figtext(0.57,0.89, r"$q^{\rm{EHT}}_{\rm{min}}=%g$" %qmin_EHT, color='black', fontsize=14)
plt.figtext(0.57,0.84, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='black', fontsize=14)
plt.figtext(0.57,0.79, r"$\dot{\mathcal{M}}=%g$" %eps, color='gray', fontsize=14)
plt.figtext(0.57,0.74, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='gray', fontsize=14)
plt.figtext(0.57,0.69, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ), color='gray', fontsize=14)
plt.figtext(0.57,0.64, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='gray', fontsize=14)
#plt.figtext(0.7,0.55, r"$L^{\rm{max}}_{mm}=10^{%g}$ erg s$^{-1}$" %Lmx_cgs, color='black', fontsize=15)


plt.figlegend([ p1[0], p2[0], p3[0], p4[0], p5[0]], ("z=%g"%Zs[0], "z=%g"%Zs[1],"z=%g"%Zs[2],"z=%g"%Zs[3],"z=%g"%Zs[4]), (0.76,0.68), fontsize=14)

#plt.figlegend([ p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0] ], ("z=%g"%Zs[0], "z=%g"%Zs[1],"z=%g"%Zs[2],"z=%g"%Zs[3],"z=%g"%Zs[4],"z=%g"%Zs[5],"z=%g"%Zs[6]), (0.76,0.58), fontsize=14)

# plt.figtext(0.2,0.26, r"$\eta^{\rm{bf}}=%g$" %popt[3], color='black', fontsize=14)

plt.xlabel(r"$\theta_{\rm{min}}$ [$\mu$as]")
if (logN):
	plt.ylabel(r"$\log_{10}[N_{\rm{EHT}}]$")
	#plt.ylabel(r"$\log_{10}[\rm{max}\left\{N_{\rm{EHT}},0.1\right\}]$")
else:
	plt.ylabel(r"$N_{\rm{EHT}}$")


plt.xlim(thMns[0], thMns[Ng-1])


plt.tight_layout()


if (TrapInt):
	if (Lmx==24.0):
		plt.title("LLAGN")
		if (logN):
			Savename = "log10NEHTmax_vs_thMn_LLAGN_TrapN%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc_Lmx%g.png" %(Ntrap_z, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ, Lmx_cgs)
		else:
			Savename = "NEHTmax_vs_thMn_LLAGN_TrapN%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc_Lmx%g.png" %(Ntrap_z, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ, Lmx_cgs)
	else:
		if (logN):
			Savename = "log10NEHTmax_vs_thMn_TrapN%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc_Lmx%g.png" %(Ntrap_z, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ, Lmx_cgs)
		else:
			Savename = "NEHTmax_vs_thMn_TrapN%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc_Lmx%g.png" %(Ntrap_z, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ, Lmx_cgs)
else:
	if (logN):
		Savename = "log10NEHTmax_vs_thMn_reclim%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc_Lmx%g.png" %(reclim, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ, Lmx_cgs)
	else:
		Savename = "NEHTmax_vs_thMn_reclim%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc_Lmx%g.png" %(reclim, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ, Lmx_cgs)

if (fEdd_Dist):
	Savename = Savename+ "_LLAGN_fDist_"

Savename = Savename.replace('.', 'p')
Savename = Savename.replace('ppng', '.png')
plt.savefig(Savename)











#if (Plt_Frcs):

plt.figure(figsize=[8,6])
plt.scatter(thMns, np.log10(NEHT_Z1/Ntot_Z1), color='#1b9e77')
plt.plot(thMns, np.log10(NEHT_Z1/Ntot_Z1), color='#1b9e77', linewidth=3)

plt.scatter(thMns, np.log10(NEHT_Z2/Ntot_Z2), color='#d95f02')
plt.plot(thMns, np.log10(NEHT_Z2/Ntot_Z2), color='#d95f02', linewidth=3)

plt.scatter(thMns, np.log10(NEHT_Z3/Ntot_Z3), color='#7570b3')
plt.plot(thMns, np.log10(NEHT_Z3/Ntot_Z3), color='#7570b3', linewidth=3)

plt.scatter(thMns, np.log10(NEHT_Z4/Ntot_Z4), color='#e7298a')
plt.plot(thMns, np.log10(NEHT_Z4/Ntot_Z4), color='#e7298a', linewidth=3)

plt.scatter(thMns, np.log10(NEHT_Z5/Ntot_Z5), color='gray')
plt.plot(thMns, np.log10(NEHT_Z5/Ntot_Z5), color='gray', linewidth=3, linestyle='--')




FminSv = Fmin/mJy2cgs/1000.
thMnSv = thMn/mu_as2rad 
PbaseSv = Pbase/yr2sec
plt.figtext(0.73,0.86, r"$q^{\rm{EHT}}_{\rm{min}}=%g$" %qmin_EHT, color='black', fontsize=14)
plt.figtext(0.73,0.81, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='black', fontsize=14)
plt.figtext(0.73,0.76, r"$\dot{\mathcal{M}}=%g$" %eps, color='gray', fontsize=14)
plt.figtext(0.73,0.71, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='gray', fontsize=14)
plt.figtext(0.73,0.66, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ), color='gray', fontsize=14)
plt.figtext(0.73,0.61, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='gray', fontsize=14)


# plt.figtext(0.2,0.26, r"$\eta^{\rm{bf}}=%g$" %popt[3], color='black', fontsize=14)

plt.xlabel(r"$\theta_{\rm{min}}$ [$\mu$as]")
if (logN):
	plt.ylabel(r"$\log_{10}[f_{\rm{EHT}}]$")
else:
	plt.ylabel(r"$N_{\rm{EHT}}$")


plt.xlim(thMns[0], thMns[Ng-1])
#plt.ylim(-6.0, 0.0)

plt.tight_layout()

if (TrapInt):
	if (logN):
		Savename = "log10fEHTmax_vs_thMn_NPq%g_TrapN%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc.png" %(Ntrp_P, Ntrap_z, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ)
	else:
		Savename = "fEHTmax_vs_thMn_NPq%g_TrapN%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc.png" %(Ntrp_P, Ntrap_z, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ)
else:
	if (logN):
		Savename = "log10fEHTmax_vs_thMn_reclim%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc.png" %(reclim, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ)
	else:
		Savename = "fEHTmax_vs_thMn_reclim%g_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g_eps%g_fEdd%g_amax%gpc.png" %(reclim, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng, eps, f_Edd, KQ)

Savename = Savename.replace('.', 'p')
Savename = Savename.replace('ppng', '.png')
plt.savefig(Savename)


















if (Fit):

	### plot best fit params
	matplotlib.rcParams.update({'font.size': 16})
	# for i in range(Nz):
	# 	if (i==0): 
	# 		zbfs = poptZ1.T[0]
	# 		epsbfs = poptZ1.T[1]
	# 		kQbfs = poptZ1.T[2]
	# 	elif (i==1):
	# 		zbfs = poptZ2.T[0]
	# 		epsbfs = poptZ2.T[1]
	# 		kQbfs = poptZ2.T[2]
	# 	else:
	# 		zbfs = poptZ3.T[0]
	# 		epsbfs = poptZ3.T[1]
	# 		kQbfs = poptZ3.T[2]

	zbfs = popt.T[0]
	epsbfs = popt.T[1]
	kQbfs = popt.T[2]
	plt.figure(figsize=[10,6])



	plt.subplot(221)
	plt.scatter(thMns, zbfs, color='black')
	plt.plot(thMns, zbfs, color='black', linewidth=3)

	FminSv = Fmin/mJy2cgs
	thMnSv = thMn/mu_as2rad 
	PbaseSv = Pbase/yr2sec
	plt.figtext(0.2,0.31, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=14)
	plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=14)
	plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=14)


	plt.xlabel(r"$\theta_{\rm{min}}$ [$\mu$as]")
	plt.ylabel(r"$z$")

	plt.xlim(thMns[0], thMns[Ng-1])



	plt.subplot(222)
	plt.scatter(thMns, np.log10(epsbfs), color='black')
	plt.plot(thMns, np.log10(epsbfs), color='black', linewidth=3)



	plt.xlabel(r"$\theta_{\rm{min}}$ [$\mu$as]")
	plt.ylabel(r"$\log_{10}[\dot{\mathcal{M}}$]")


	plt.xlim(thMns[0], thMns[Ng-1])




	plt.subplot(223)
	plt.scatter(thMns, np.log10(kQbfs), color='black')
	plt.plot(thMns, np.log10(kQbfs), color='black', linewidth=3)



	plt.xlabel(r"$\theta_{\rm{min}}$ [$\mu$as]")
	plt.ylabel(r"$\log_{10}[a_{\rm{max}}/\rm{pc}]$")


	plt.xlim(thMns[0], thMns[Ng-1])

	plt.tight_layout()





	plt.subplot(224)
	if (fit_MdEff):
		Mdeffbfs = popt.T[3]
		plt.scatter(thMns, np.log10(Mdeffbfs), color='black')
		plt.plot(thMns, np.log10(Mdeffbfs), color='black', linewidth=3)
	else:
		plt.scatter(thMns, np.log10(MdEff*thMns/thMns), color='black')
		plt.plot(thMns, np.log10(MdEff*thMns/thMns), color='black', linewidth=3)


	plt.xlabel(r"$\theta_{\rm{min}}$ [$\mu$as]")
	plt.ylabel(r"$\log_{10}[\eta$]")


	plt.xlim(thMns[0], thMns[Ng-1])

	plt.tight_layout()


	icnt = i+1
	Savename = "BestFitParams_Z%g_vs_thMn_qminEHT%g_qminPOP%g_Fmin%gmJy_Pbase%gyr_N%g.png" %(icnt, qmin_EHT, qmin_POP, FminSv, PbaseSv, Ng)


	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)






