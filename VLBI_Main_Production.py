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
# import VLBI_IntFuncs_V2 as IFs 
# from VLBI_IntFuncs_V2 import *
# import GWB_IntFuncs as hGWB 
# from GWB_IntFuncs import *

import GWB_IntFuncs_fEddDist as hGWB 
from GWB_IntFuncs_fEddDist import *

###PLOTTING OPTIONS
CntPlt_CumZ = True
CntPlt_DZ = False

##INT options
TrapInt = True



## method optioncs
fEdd_Dist = True


if (fEdd_Dist):
	import VLBI_IntFuncs_fEddDist as IFs 
	from VLBI_IntFuncs_fEddDist import *
else:
	import VLBI_IntFuncs_V2 as IFs 
	from VLBI_IntFuncs_V2 import *





#low dep on Mdot, got rid of fedd with dist func, qmin and Npc remain
Mmx_Pbase = True
fEdd_Npc = False
Mdot_Npc = True
qmin_Npc = True

qmin_Mdot = False
qmin_fEdd = False
fEdd_Mdot = False


if (fEdd_Dist):
	fEdd_Npc = False
	fEdd_Mdot = False
	qmin_fEdd = False







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
##Instrument params
Fmin = 1.0 * mJy2cgs
thMn = 20.0 * mu_as2rad 
Pbase = 10.0*yr2sec














#Accretion Params
KQ = 10.**(-1.0) ## sets number of pc at which RL turns on
eps = 1.0#10**(-3.75)  ## sets migration (accretion rate in CBD pushing two together)

# if (Lmx==24.0):	
# 	f_Edd = 0.0001
# 	fbin = 1.00
# else:
# 	f_Edd = 0.1  ## sets L to M connection (accretion rate onto shining BH)

fbin = 1.0
f_Edd = 10.**(-4.12)  
MdEff = 0.1

## binary pop params
qmin_EHT = 0.01   ### qminof EHT sample
qmin_POP = np.minimum(qmin_EHT, 0.01)  ### qmin of all MBHBS 

zeval = 0.2  #eval at this z
zmax = 3.0 ### integrateo out to zmax=5.0














###Cosmology
h=0.7
Om = 0.3
OL=0.7




###
Mmx = 1.*10.**10 ## This is for LLAGN - choose L, and randomly draw f_edd can give very large inferred masses - so limit by obs knowledge of BH mass
Mmax = 1.*10.**10*Msun
Mmin= 0.0*10.**5*Msun ## no really necessary, but in case wan tot limit this



###SED scaling
chi  = 0.5 #Elvis 1994 Radio Loud AGN - nuL_nu propto nu^(0.9) and #0.1  #chi L_(408MHz) = L_mm (specifc fluxes)



###DEFUNCT
DZ = 1.0
##FREE (and KEY) PARAMETERS
##Overall pop params (keep 1)
xi = 1.0







nuVbnd = c/(5.45*10**(-5))
F0opt  = 3.636*10**(-20)*nuVbnd 
maglim = 24.5
Fmin_opt = 10.**(-maglim/2.5)*F0opt/nuVbnd 

fPTA = 7. * 10.**(-9) ##see sesana, haiman, Kocsis, kelley
hPTA = 3.*10.**(-15)  # best case at above freq from PPTA



if (CntPlt_CumZ):
	Ng = 10
	#DZ = 0.2
	#zeval = 0.37  ##Max z count

	fEdds = np.linspace(-5.0, 1.0, Ng)
	epss = np.linspace(-3.0, 2.0, Ng)
	KQs = np.linspace(-2.5, 1.0, Ng)
	thMns = np.linspace(-1.0, 2.0, Ng) 
	MdEffs = np.linspace(-2., 0., Ng)
	Pbasez = np.linspace(0.0, np.log10(2.*Pbase/yr2sec), Ng)
	Mmxz = np.linspace(5.0, 10.5, Ng)
	zs = np.linspace(-2.0, 0.0, Ng)
	qmins = np.linspace(-3., 0.0, Ng)



	if (fEdd_Mdot):

				## CONTOUR PLOT TOTALs

		Ntot_grid = np.zeros([Ng,Ng])
		hGWB_grid = np.zeros([Ng,Ng])
		


		if (TrapInt):
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  np.maximum(1.e-3,-IntzZ_Trap_OptNEHT([10.**epss[i], KQ], zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, 10.**fEdds[j], xi, fbin, h, Om, OL))

					hGWB_grid[j][i] =  -IntzZ_Trap_GWB_f([10.**epss[i],  KQ], fPTA, zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, 10.**fEdds[j], xi, fbin, h, Om, OL)



		else:
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_OptNEHT([10.**epss[i], KQ], zmax, Mmx, Fmin, chi, thMn, qmin, Pbase, 10.**fEdds[j], xi, fbin, h, Om, OL))
					


		



		fig = plt.figure(figsize=[7.5,6.1])
		cnt = plt.contourf(epss, fEdds, np.log10(Ntot_grid), 200, cmap = "viridis")


		#cbar = plt.colorbar(cnt)
		N_lms = plt.contour(epss, fEdds, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		plt.clabel(N_lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		GWB_lms = plt.contour(epss, fEdds, np.log10(hGWB_grid), colors='white', levels = [-16, -15, -14.5, -14, 13])
		plt.clabel(GWB_lms, fmt = r'$10^{%g}$', colors = 'white', fontsize=14)	



		plt.axvline(x=np.log10(eps), color='chartreuse', linewidth=2, linestyle="--")
		plt.axhline(y=np.log10(f_Edd), color='chartreuse', linewidth=2, linestyle="--")
		plt.scatter(np.log10(eps), np.log10(f_Edd), color='chartreuse', marker='o', s=30)



		plt.ylabel(r'$\rm{log}_{10}[f_{\rm{Edd}}]$')
		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')



		FminSv = Fmin/mJy2cgs/1000.
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		Lmx_cgs = Lmx +7.0
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)

		#plt.figtext(0.15,0.86, r"$L^{\rm{max}}_{mm}=10^{%g}$ erg s$^{-1}$" %Lmx_cgs, color='yellow', fontsize=15)
		plt.figtext(0.18,0.87, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
		plt.figtext(0.18,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		if (FminSv*1000.0 <= 100):
			FminSv = FminSv*1000.0
			plt.figtext(0.18,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		else:
			plt.figtext(0.18,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.18,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
		plt.figtext(0.18,0.67, r"$a_{\rm{max}}=%g$ pc" %KQ, color='yellow', fontsize=15)
		plt.figtext(0.18,0.62, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)


		#plt.figtext(0.15,0.66, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)


		#plt.figtext(0.2,0.61, r"$\chi=%g$" %chi, color='yellow', fontsize=15)


		#plt.ylim(KQs[0], KQs[len(KQs)-1])
		#plt.xlim(fEdds[0], fEdds[len(fEdds)-1])

		plt.tight_layout()

		Savename = 'save'
		if (TrapInt):
			Savename = 'CumZ_Mdot_vs_fEdd_%gx%g_amax%gpc_qminEHT%g_qminPOP%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_TrapInt%g.png'%(Ng,Ng, KQ, qmin_EHT, qmin_POP, zmax, FminSv, thMnSv, PbaseSv, Ntrap_z)
		else:
			Savename = 'CumZ_Mdot_vs_fEdd_%gx%g_amax%gpc_qminEHT%g_qminPOP%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, KQ, qmin_EHT, qmin_POP, zmax, FminSv, thMnSv, PbaseSv, reclim)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)







	if (fEdd_Npc):

				## CONTOUR PLOT TOTALs

		Ntot_grid = np.zeros([Ng,Ng])
		hGWB_grid = np.zeros([Ng,Ng])

		ttots_mn = np.zeros([Ng,Ng])
		ttots_mx = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)

		RSGmxG = np.zeros([Ng,Ng])
		RSGmnG = np.zeros([Ng,Ng])

		aTmx = np.zeros(Ng)
		aTmn = np.zeros(Ng)

		if (TrapInt):
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  np.maximum(1.e-3,-IntzZ_Trap_OptNEHT([eps, 10.**KQs[j]], zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, 10.**fEdds[i], xi, fbin, h, Om, OL))
					#hGWB_grid[j][i] =  -IntzZ_Trap_GWB_f([eps,  10.**KQs[j]], fPTA, zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, MdEff, 10.**fEdds[i], xi, fbin, h, Om, OL)


				
				
				RSGmx[i] = RSGff(10.**epss[i], 1.e9*Msun, MdEff)/pc2cm
				RSGmn[i] = RSGes(10.**epss[i], 1.e9*Msun, MdEff)/pc2cm
				#aTmx[i] = aTrans(10.**10*Msun, 1.0, MdEff, 10.**epss[i])/pc2cm
				aTmn[i] = aTrans(10.**9*Msun, qsofq(qmin_EHT), MdEff, 10.**epss[i])/pc2cm

		else:
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_OptNEHT([eps, 10.**KQs[j]], zmax, Mmx, Fmin, chi, thMn, qmin, Pbase, 10.**fEdds[i], xi, fbin, h, Om, OL))
					

				RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm
				aTmx[i] = aTrans(10.**10*Msun, 1.0, MdEff, 10.**epss[i])/pc2cm
				aTmn[i] = aTrans(10.**6*Msun, qsofq(qmin_EHT), MdEff, 10.**epss[i])/pc2cm


		

		Log_thmxs = np.log10(10.**KQs*pc2cm/Dang(zeval, h, Om, OL) * 180./np.pi * 3600 * 1.e6)


		fig = plt.figure(figsize=[7.5,6.1])
		ax = fig.add_subplot(111)
		# if (Lmx==24.0):
		# 	#plt.title(r'LLAGN, $\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$' %zmax, fontsize=15)
		# 	plt.title(r'LLAGN, $\log_{10}{\left[ N_{\rm{VLBI}} \right]}$' fontsize=15)
		# else:
		# 	#plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$' %zmax, fontsize=15)
		plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$', fontsize=16)

		ax.contourf(fEdds, Log_thmxs, np.log10(Ntot_grid), 200, cmap = "viridis")
		ax2 = ax.twinx()
		ax2.contourf(fEdds, KQs, np.log10(Ntot_grid), 200, cmap = "viridis")

		#cbar = plt.colorbar(cnt)
		lms = plt.contour(fEdds, KQs, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])		
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=15)	

		# GWB_lms = plt.contour(epss, fEdds, np.log10(hGWB_grid), colors='white', levels = [-16, -15, -14.5, -14, 13])
		# plt.clabel(GWB_lms, fmt = r'$10^{%g}$', colors = 'white', fontsize=14)	


		plt.plot(fEdds, np.log10(RSGmx), color='red', linestyle=":")
		plt.plot(fEdds, np.log10(RSGmn), color='red', linestyle=":")
		plt.fill_between(fEdds, np.log10(RSGmx),np.log10(RSGmn), color='red', alpha=0.1)

		#plt.plot(epss, np.log10(aTmx), color='green', linewidth=2, linestyle=":" )
		#plt.plot(epss, np.log10(aTmn), color='yellow', linewidth=2, linestyle="--" )




		#plt.axvline(x=0.0, color='black', linewidth=1, linestyle="--")


		plt.axvline(x=np.log10(f_Edd), color='chartreuse', linewidth=2, linestyle="--")
		plt.axhline(y=np.log10(KQ), color='chartreuse', linewidth=2, linestyle="--")
		plt.scatter(np.log10(f_Edd), np.log10(KQ), color='chartreuse', marker='o', s=30)



		ax.set_xlabel(r'$\rm{log}_{10}[f_{\rm{Edd}}]$')
		#plt.ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
		ax2.set_ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
		ax.set_ylabel(r'$\rm{log}_{10}[\theta_{\rm{max}}(z=%g)/\mu\rm{as}]$' %zeval)


		FminSv = Fmin/mJy2cgs/1000.
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		Lmx_cgs = Lmx +7.0
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)

		#plt.figtext(0.15,0.86, r"$L^{\rm{max}}_{mm}=10^{%g}$ erg s$^{-1}$" %Lmx_cgs, color='yellow', fontsize=15)
		plt.figtext(0.15,0.87, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
		plt.figtext(0.15,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		if (FminSv*1000.0 <= 100):
			FminSv = FminSv*1000.0
			plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		else:
			plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.15,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
		plt.figtext(0.15,0.67, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
		plt.figtext(0.15,0.62, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)

		#plt.figtext(0.15,0.66, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)

		#plt.figtext(0.2,0.61, r"$\chi=%g$" %chi, color='yellow', fontsize=15)


		plt.ylim(KQs[0], KQs[len(KQs)-1])
		plt.xlim(fEdds[0], fEdds[len(fEdds)-1])

		plt.tight_layout()

		Savename = 'save'
		if (TrapInt):
			if (Lmx==24.0):
				Savename = 'CumZ_Npc_vs_fEdd_LLAGN_%gx%g_eps%g_qminEHT%g_qminPOP%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_TrapInt%g.png'%(Ng,Ng, eps, qmin_EHT, qmin_POP, zmax, FminSv, thMnSv, PbaseSv, Ntrap_z)
			else:
				Savename = 'CumZ_Npc_vs_fEdd_%gx%g_eps%g_qminEHT%g_qminPOP%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_TrapInt%g.png'%(Ng,Ng, eps, qmin_EHT, qmin_POP, zmax, FminSv, thMnSv, PbaseSv, Ntrap_z)
		else:
			Savename = 'CumZ_Npc_vs_fEdd_%gx%g_eps%g_qminEHT%g_qminPOP%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, eps, qmin_EHT, qmin_POP, zmax, FminSv, thMnSv, PbaseSv, reclim)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)








	if (Mmx_Pbase):
		print "Mmx vs Pbase"

		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])

		ttots_mn = np.zeros([Ng,Ng])
		ttots_mx = np.zeros([Ng,Ng])

		aseps = np.zeros([Ng,Ng])
		lgstrns = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)

		if (TrapInt):
			if (Lmx==24.0):
				#f_Edd = 10.**(-4)
				for i in range(0,Ng):
					for j in range(0,Ng):		
						Ntot_grid[j][i] =  max(1.e-3,-IntzZ_Trap_OptNEHT([eps, KQ], zmax, 10.**Mmxz[i], Fmin, chi, thMn, qmin_EHT, qmin_POP, 10.**Pbasez[j]*yr2sec, f_Edd, xi, fbin, h, Om, OL))
						aseps[j][i] = asep(10.**Pbasez[j]*yr2sec, 10.**Mmxz[i]*Msun)/pc2cm
						ttots_mn[j][i] = t_tot(aseps[j][i], 10.**Mmxz[i]*Msun, qsofq(0.01), MdEff, eps)/yr2sec
						ttots_mx[j][i] = t_tot(aseps[j][i], 10.**Mmxz[i]*Msun, qsofq(1.0), MdEff, eps)/yr2sec
						lgstrns[j][i] = np.log10(hGWB.hPTA(10.**Pbasez[j]*yr2sec, 10.**Mmxz[i]*Msun, 0.1, zeval, h, Om, OL))

			else:
				for i in range(0,Ng):
					for j in range(0,Ng):		
						Ntot_grid[j][i] =  max(1.e-3,-IntzZ_Trap_OptNEHT([eps, KQ], zmax, 10.**Mmxz[i], Fmin, chi, thMn, qmin_EHT, qmin_POP, 10.**Pbasez[j]*yr2sec, f_Edd, xi, fbin, h, Om, OL))
						aseps[j][i] = asep(10.**Pbasez[j]*yr2sec, 10.**Mmxz[i]*Msun)/pc2cm
						ttots_mn[j][i] = t_tot(aseps[j][i], 10.**Mmxz[i]*Msun, qsofq(0.01), MdEff, eps)/yr2sec
						ttots_mx[j][i] = t_tot(aseps[j][i], 10.**Mmxz[i]*Msun, qsofq(1.0), MdEff, eps)/yr2sec
						lgstrns[j][i] = np.log10(hGWB.hPTA(10.**Pbasez[j]*yr2sec, 10.**Mmxz[i]*Msun, 0.1, zeval, h, Om, OL))



		else:
			for i in range(0,Ng):
				for j in range(0,Ng):		
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_OptNEHT([eps, KQ], zmax, 10.**Mmxz[i], Fmin, chi, thMn, qmin, 10.**Pbasez[j]*yr2sec, f_Edd, xi, fbin, h, Om, OL))
					aseps[j][i] = asep(10.**Pbasez[j]*yr2sec, 10.**Mmxz[i]*Msun)/pc2cm

		PGWq0p01 = PTrans(10.**Mmxz*Msun, qsofq(0.01), MdEff, eps)/yr2sec
		PGWq1    = PTrans(10.**Mmxz*Msun, qsofq(1.0), MdEff, eps)/yr2sec
		Pminz0p1 = PminRes(10.**Mmxz*Msun, thMn, 0.1, h, Om, OL)/yr2sec
		Pminz1 = PminRes(10.**Mmxz*Msun, thMn, 0.5, h, Om, OL)/yr2sec

		



		fig = plt.figure(figsize=[7.5,6.1])
		# if (Lmx==24.0):
		# 	plt.title(r"LLAGN, $\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$" %(zmax), fontsize=20)
		# else:
		# 	plt.title(r"$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$" %(zmax), fontsize=20)

		plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$', fontsize=16)

		cnt = plt.contourf(Mmxz, Pbasez, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		#cbar = plt.colorbar(cnt)
		lms = plt.contour(Mmxz, Pbasez, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		# alms= plt.contour(Mmxz, Pbasez, np.log10(aseps), cmap = "viridis", levels = [-3.0, -2.0, -1.0, 0.0])
		# plt.clabel(alms, fmt = r'$10^{%g}$', colors = 'blue', fontsize=14)

		alms= plt.contour(Mmxz, Pbasez, aseps, colors="white", levels = [0.001, 0.003, 0.01, 0.03, 0.1])
		plt.clabel(alms, fmt = r'$a=%g$ pc', colors="white", fontsize=14)


		alms= plt.contour(Mmxz, Pbasez, lgstrns, colors="cyan", levels = [-26., -24., -22., -20., -18., -16., -15., -14., -13., -12., -11.])
		plt.clabel(alms, fmt = r'$h=10^{%g}$', colors="cyan", fontsize=14)

		#plt.plot(Mmxz, np.log10(PGWq0p01), color="yellow", linestyle=":")
		#plt.plot(Mmxz, np.log10(PGWq1), color="yellow", linestyle="--")

		#plt.plot(Mmxz, np.log10(Pminz0p1), color="cyan", linestyle=":")
		#plt.plot(Mmxz, np.log10(Pminz1), color="cyan", linestyle="--")
		#plt.plot(Mmxz, np.log10(ttots_mn/yr2sec), color="yellow")
		#plt.plot(Mmxz, np.log10(ttots_mx/yr2sec), color="cyan")
		lmsmn = plt.contour(Mmxz, Pbasez, np.log10(ttots_mn), colors="yellow", levels = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0 ])
		lmsmx = plt.contour(Mmxz, Pbasez, np.log10(ttots_mx), colors="cyan", levels = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0 ])
		plt.clabel(lmsmn, fmt = r'$t_{\rm{res}} = 10^{%g} yr$', colors="yellow", fontsize=12, linestyle='--')	
		plt.clabel(lmsmx, fmt = r'$t_{\rm{res}} = 10^{%g} yr$', colors="cyan", fontsize=12, linestyle=':')	


		plt.ylabel(r'$\rm{log}_{10}[P_{\rm{base}}/yr]$')
		plt.xlabel(r'$\rm{log}_{10}[M_{\rm{max}}/M_{\odot}]$')



		FminSv = Fmin/mJy2cgs/1000.
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		Lmx_cgs = Lmx + 7.0


		if (fEdd_Dist):
			Savename = 'Cumz_Pbase_vs_Mmx_fDistLLAGN_%gx%g_qminEHT%g_qminPOP%g_amax%g_eps%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_zmax%g_Lmx%g_Trap%g.png'%(Ng,Ng, qmin_EHT, qmin_POP, KQ, eps, FminSv, thMnSv, PbaseSv, zmax, Lmx, Ntrap_z)

			#plt.figtext(0.2,0.51, r"$\chi=%g$" %chi, color='yellow', fontsize=15)
			#plt.figtext(0.15,0.52, r"$L^{\rm{max}}_{mm}=10^{%g}$ erg s$^{-1}$" %Lmx_cgs, color='yellow', fontsize=15)
			plt.figtext(0.15,0.47, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
			plt.figtext(0.15,0.42, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.15,0.37, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.15,0.37, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.32, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
			plt.figtext(0.15,0.27, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ), color='yellow', fontsize=15)
			plt.figtext(0.15,0.22, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)
		else:
			Savename = 'Cumz_Pbase_vs_Mmx_fdelta_%gx%g_qminEHT%g_qminPOP%g_amax%g_eps%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_zmax%g_Lmx%g_Trap%g.png'%(Ng,Ng, qmin_EHT, qmin_POP, KQ, eps, FminSv, thMnSv, PbaseSv, zmax, Lmx, Ntrap_z)	#plt.figtext(0.2,0.51, r"$\chi=%g$" %chi, color='yellow', fontsize=15)
			
			#plt.figtext(0.15,0.52, r"$L^{\rm{max}}_{mm}=10^{%g}$ erg s$^{-1}$" %Lmx_cgs, color='yellow', fontsize=15)
			plt.figtext(0.15,0.47, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
			plt.figtext(0.15,0.42, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.15,0.37, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.15,0.37, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.32, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
			plt.figtext(0.15,0.27, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ), color='yellow', fontsize=15)
			plt.figtext(0.15,0.22, r"$f_{\rm{Edd}}=10^{%g}$" %np.log10(f_Edd), color='yellow', fontsize=15)
			plt.figtext(0.15,0.17, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)


		plt.ylim(Pbasez[0], Pbasez[len(KQs)-1])

		plt.tight_layout()


		# if (TrapInt):
		# 	if (Lmx==24.0):
		# 		Savename = 'Cumz_Pbase_vs_Mmx_LLAGN_%gx%g_qminEHT%g_qminPOP%g_amax%g_eps%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_zmax%g_Lmx%g_Trap%g.png'%(Ng,Ng, qmin_EHT, qmin_POP, KQ, eps, FminSv, thMnSv, PbaseSv, zmax, Lmx, Ntrap_z)
		# 	else:
		# 		Savename = 'Cumz_Pbase_vs_Mmx_%gx%g_qminEHT%g_qminPOP%g_amax%g_eps%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_zmax%g_Lmx%g_Trap%g.png'%(Ng,Ng, qmin_EHT, qmin_POP, KQ, eps, FminSv, thMnSv, PbaseSv, zmax, Lmx, Ntrap_z)

		# else:
		# 	Savename = 'Cumz_Pbase_vs_Mmx_%gx%g_qminEHT%g_qminPOP%g_amax%g_eps%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_zmax%g_Lmx%g_reclim%g.png'%(Ng,Ng, qmin_EHT, qmin_POP, KQ, eps, FminSv, thMnSv, PbaseSv, zmax, Lmx, reclim)

		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
















	if (Mdot_Npc):

		## CONTOUR PLOT TOTALs

		Ntot_grid = np.zeros([Ng,Ng])
		ttots_mn = np.zeros([Ng,Ng])
		ttots_mx = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)

		aTmx = np.zeros(Ng)
		aTmn = np.zeros(Ng)

		if (TrapInt):
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_Trap_OptNEHT([10.**epss[i], 10.**KQs[j]], zmax, Mmx, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, f_Edd, xi, fbin, h, Om, OL))
					ttots_mn[j][i] = t_tot(10.**KQs[j]*pc2cm, 10.**9.0*Msun, qsofq(0.1), MdEff, 10.**epss[i])/yr2sec
					#ttots_mx[j][i] = t_tot(10.**KQs[j]*pc2cm, 10.**9*Msun, qsofq(0.1), MdEff, 10.**epss[i])/yr2sec

				RSGmx[i] = RSGff(10.**epss[i], 1.e9*Msun, MdEff)/pc2cm
				RSGmn[i] = RSGes(10.**epss[i], 1.e9*Msun, MdEff)/pc2cm
				#RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm
				#aTmx[i] = aTrans(10.**10*Msun, 1.0, MdEff, 10.**epss[i])/pc2cm
				aTmn[i] = aTrans(10.**8*Msun, qsofq(qmin_EHT), MdEff, 10.**epss[i])/pc2cm

		else:
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_OptNEHT([10.**epss[i], 10.**KQs[j]], zmax, Mmx, Fmin, chi, thMn, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL))
					ttots[j][i] = t_tot(10.**KQs[j]*pc2cm, 10.**9*Msun, 1.0, MdEff, 10.**epss[i])/yr2sec

				RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm
				aTmx[i] = aTrans(10.**10*Msun, 1.0, MdEff, 10.**epss[i])/pc2cm
				aTmn[i] = aTrans(10.**6*Msun, qsofq(qmin_EHT), MdEff, 10.**epss[i])/pc2cm


		

		Log_thmxs = np.log10(10.**KQs*pc2cm/Dang(zeval, h, Om, OL) * 180./np.pi * 3600 * 1.e6)


		fig = plt.figure(figsize=[7.5,6.1])
		ax = fig.add_subplot(111)

		#plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$' %zmax, fontsize=15)
		plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$', fontsize=16)

		ax.contourf(epss, Log_thmxs, np.log10(Ntot_grid), 200, cmap = "viridis")
		ax2 = ax.twinx()
		ax2.contourf(epss, KQs, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		#cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, KQs, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=15)	

		plt.plot(epss, np.log10(RSGmx), color='red' )
		plt.plot(epss, np.log10(RSGmn), color='red' )
		plt.fill_between(epss,np.log10(RSGmx),np.log10(RSGmn), color='red', alpha=0.1)

		#plt.plot(epss, np.log10(aTmx), color='green', linewidth=2, linestyle=":" )
		#plt.plot(epss, np.log10(aTmn), color='yellow', linewidth=2, linestyle="--" )

		lmsmn = plt.contour(epss, KQs, np.log10(ttots_mn), colors="white", levels = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
		#lmsmx = plt.contour(epss, KQs, np.log10(ttots_mx), colors="cyan", levels = [7.0, 8.0])
		plt.clabel(lmsmn, fmt = r'$10^{%g}$', colors="white", fontsize=12, linestyle='--')	
		#plt.clabel(lmsmx, fmt = r'$10^{%g}$', colors="cyan", fontsize=12, linestyle=':')	



		#plt.axvline(x=0.0, color='black', linewidth=1, linestyle="--")


		plt.axvline(x=np.log10(eps), color='chartreuse', linewidth=2, linestyle="--")
		plt.axhline(y=np.log10(KQ), color='chartreuse', linewidth=2, linestyle="--")
		plt.scatter(np.log10(eps), np.log10(KQ), color='chartreuse', marker='o', s=30)



		ax.set_xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		#plt.ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
		ax2.set_ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
		ax.set_ylabel(r'$\rm{log}_{10}[\theta_{\rm{max}}(z=%g)/\mu\rm{as}]$' %zeval)


		FminSv = Fmin/mJy2cgs/1000.
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		Lmx_cgs = Lmx +7.0
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)

		#plt.figtext(0.15,0.86, r"$L^{\rm{max}}_{mm}=10^{%g}$ erg s$^{-1}$" %Lmx_cgs, color='yellow', fontsize=15)
		
		if (fEdd_Dist):
			plt.figtext(0.15,0.87, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
			plt.figtext(0.15,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.67, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)
		else:
			plt.figtext(0.15,0.87, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
			plt.figtext(0.15,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.67, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)
			plt.figtext(0.15,0.62, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)

		#plt.figtext(0.2,0.61, r"$\chi=%g$" %chi, color='yellow', fontsize=15)


		plt.ylim(KQs[0], KQs[len(KQs)-1])
		plt.xlim(epss[0], epss[len(epss)-1])

		plt.tight_layout()

		Savename = 'save'
		if (TrapInt):
			Savename = 'CumZ_Npc_vs_Mdot_%gx%g_qminEHT%g_qminPOP%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_TrapInt%g.png'%(Ng,Ng, qmin_EHT, qmin_POP, zmax, FminSv, thMnSv, PbaseSv, Ntrap_z)
		else:
			Savename = 'CumZ_Npc_vs_Mdot_%gx%g_qminEHT%g_qminPOP%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, qmin_EHT, qmin_POP, zmax, FminSv, thMnSv, PbaseSv, reclim)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)



if (qmin_Npc):

		## CONTOUR PLOT TOTALs
		ttots_mn = np.zeros([Ng,Ng])
		ttots_mx = np.zeros([Ng,Ng])
		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)

		aTmx = np.zeros(Ng)
		aTmn = np.zeros(Ng)

		if (TrapInt):
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_Trap_OptNEHT([eps, 10.**KQs[j]], zmax, Mmx, Fmin, chi, thMn, 10.**qmins[i], np.minimum(qmin_POP, 10.**qmins[i]), Pbase, f_Edd, xi, fbin, h, Om, OL))
					ttots_mn[j][i] =   t_tot(10.**KQs[j]*pc2cm, 10.**9.0*Msun, 10.**qmins[i], MdEff, eps)/yr2sec


				RSGmx[i] = RSGff(10.**epss[i], 1.e9*Msun, MdEff)/pc2cm
				RSGmn[i] = RSGes(10.**epss[i], 1.e9*Msun, MdEff)/pc2cm
				
				
		else:
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_OptNEHT([10.**epss[i], 10.**KQs[j]], zmax, Mmx, Fmin, chi, thMn, 10.**qmins[i], np.minimum(qmin_POP, 10.**qmins[i]), Pbase, f_Edd, xi, fbin, h, Om, OL))
					ttots_mn[j][i] =   t_tot(10.**KQs[j]*pc2cm, 10.**9.0*Msun, 10.**qmins[i], MdEff, eps)/yr2sec

				#RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm
				#aTmx[i] = aTrans(Mmax, qsofq(10**qmins), MdEff, 10.**epss[i])/pc2cm
				#aTmn[i] = aTrans(Mmin, qsofq(10**qmins), MdEff, 10.**epss[i])/pc2cm

		Log_thmxs = np.log10(10.**KQs*pc2cm/Dang(zeval, h, Om, OL) * 180./np.pi * 3600 * 1.e6)


		fig = plt.figure(figsize=[7.5,6.1])
		ax = fig.add_subplot(111)
		#plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$' %zmax, fontsize=15)
		plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$', fontsize=16)

		ax.contourf(qmins, Log_thmxs, np.log10(Ntot_grid), 200, cmap = "viridis")
		ax2 = ax.twinx()
		ax2.contourf(qmins, KQs, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		#cbar = plt.colorbar(cnt)
		lms = plt.contour(qmins, KQs, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		# plt.plot(qmins, np.log10(RSGmx), color='red' )
		# plt.plot(qmins, np.log10(RSGmn), color='red' )
		# plt.fill_between(qmins,np.log10(RSGmx),np.log10(RSGmn), color='red', alpha=0.3)

		# plt.plot(qmins, np.log10(aTmx), color='green', linewidth=2, linestyle="--" )
		# plt.plot(qmins, np.log10(aTmn), color='green', linewidth=2, linestyle="--" )

		plt.plot(epss, np.log10(RSGmx), color='red' )
		plt.plot(epss, np.log10(RSGmn), color='red' )
		plt.fill_between(epss,np.log10(RSGmx),np.log10(RSGmn), color='red', alpha=0.1)


		plt.axvline(x=0.0, color='black', linewidth=1, linestyle="--")


		plt.axvline(x=np.log10(qmin_EHT), color='chartreuse', linewidth=2, linestyle="--")
		plt.axhline(y=np.log10(KQ), color='chartreuse', linewidth=2, linestyle="--")
		plt.scatter(np.log10(qmin_EHT), np.log10(KQ), color='chartreuse', marker='o', s=30)

		lmsmn = plt.contour(qmins, KQs, np.log10(ttots_mn), colors="white", levels = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
		plt.clabel(lmsmn, fmt = r'$10^{%g}$', colors="white", fontsize=12, linestyle='--')	





		ax.set_xlabel(r'$\rm{log}_{10}[q^{\rm{Vmin}}_s]$')
		#plt.ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
		ax2.set_ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
		ax.set_ylabel(r'$\rm{log}_{10}[\theta_{\rm{max}}/\mu\rm{as}]$')


		FminSv = Fmin/mJy2cgs/1000.
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		if (fEdd_Dist):
			Savename = 'CumZ_Npc_vs_qmin_fdistLLAGN_%gx%g_Mdot%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_Lmx%g_TrapInt%g.png'%(Ng,Ng, eps, zmax, FminSv, thMnSv, PbaseSv, Lmx, Ntrap_z)
			
			plt.figtext(0.2,0.87, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
			plt.figtext(0.2,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.2,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
			plt.figtext(0.2,0.67, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)
		else:
			Savename = 'CumZ_Npc_vs_qmin_fdelta_%gx%g_Mdot%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_Lmx%g_TrapInt%g.png'%(Ng,Ng, eps, zmax, FminSv, thMnSv, PbaseSv, Lmx, Ntrap_z)

			plt.figtext(0.15,0.87, r"$q^{\rm{Vmin}}_{s}=%g$" %qmin_EHT, color='yellow', fontsize=15)
			plt.figtext(0.15,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.15,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
			plt.figtext(0.15,0.67, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)
			plt.figtext(0.15,0.62, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)


			#plt.figtext(0.2,0.61, r"$\chi=%g$" %chi, color='yellow', fontsize=15)


		plt.ylim(KQs[0], KQs[len(KQs)-1])
		plt.xlim(qmins[0], qmins[len(epss)-1])

		plt.tight_layout()


		# if (TrapInt):
		# 	Savename = 'CumZ_Npc_vs_qmin_%gx%g_Mdot%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_Lmx%g_TrapInt%g.png'%(Ng,Ng, eps, zmax, FminSv, thMnSv, PbaseSv, Lmx, Ntrap_z)
		# else:
		# 	Savename = 'CumZ_Npc_vs_qmin_%gx%g_Mdot%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, eps, zmax, FminSv, thMnSv, PbaseSv, reclim)
		
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)














if (qmin_Mdot):

		## CONTOUR PLOT TOTALs

		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)

		aTmx = np.zeros(Ng)
		aTmn = np.zeros(Ng)

		if (TrapInt):
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_Trap_OptNEHT([10.**epss[j], KQ], zmax, Mmx, Fmin, chi, thMn, 10.**qmins[i], np.minimum(qmin_POP, 10.**qmins[i]), Pbase, f_Edd, xi, fbin, h, Om, OL))
		else:
		 print "Not sure what you doing here"	

	



		fig = plt.figure(figsize=[7.5,6.1])
		ax = fig.add_subplot(111)
		#plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$' %zmax, fontsize=15)
		plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$', fontsize=16)

		ax.contourf(qmins, epss, np.log10(Ntot_grid), 200, cmap = "viridis")

		lms = plt.contour(qmins, epss, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
	
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		

		# plt.axvline(x=0.0, color='black', linewidth=1, linestyle="--")


		plt.axvline(x=np.log10(qmin_EHT), color='chartreuse', linewidth=2, linestyle="--")
		plt.axhline(y=np.log10(eps), color='chartreuse', linewidth=2, linestyle="--")
		plt.scatter(np.log10(qmin_EHT), np.log10(eps), color='chartreuse', marker='o', s=30)



		ax.set_xlabel(r'$\rm{log}_{10}[q_{\rm{min}}]$')
		ax.set_ylabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')



		FminSv = Fmin/mJy2cgs/1000.
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)

		#plt.figtext(0.2,0.87, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
		if (fEdd_Dist):
			plt.figtext(0.2,0.87, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ),  color='yellow', fontsize=15)
			plt.figtext(0.2,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.2,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
			plt.figtext(0.2,0.67, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)
		else:
			plt.figtext(0.2,0.87, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ),  color='yellow', fontsize=15)
			plt.figtext(0.2,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
			if (FminSv*1000.0 <= 100):
				FminSv = FminSv*1000.0
				plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
			else:
				plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
			plt.figtext(0.2,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
			plt.figtext(0.2,0.67, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)
			plt.figtext(0.2,0.62, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)

		#plt.figtext(0.2,0.61, r"$\chi=%g$" %chi, color='yellow', fontsize=15)


		plt.ylim(epss[0], epss[len(epss)-1])
		plt.xlim(qmins[0], qmins[len(qmins)-1])

		plt.tight_layout()


		if (TrapInt):
			Savename = 'CumZ_Mdot_vs_qmin_%gx%g_fEdd%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_TrapInt%g.png'%(Ng,Ng, f_Edd, zmax, FminSv, thMnSv, PbaseSv, Ntrap_z)
		else:
			Savename = 'CumZ_Mdot_vs_qmin_%gx%g_fEdd%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, f_Edd, zmax, FminSv, thMnSv, PbaseSv, reclim)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)















if (qmin_fEdd):

		## CONTOUR PLOT TOTALs

		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)

		aTmx = np.zeros(Ng)
		aTmn = np.zeros(Ng)

		if (TrapInt):
			for i in range(0,Ng):
				for j in range(0,Ng):
					Ntot_grid[j][i] =  max(1.e-3,-IntzZ_Trap_OptNEHT([10.**eps, KQ], zmax, Mmx, Fmin, chi, thMn, 10.**qmins[i], np.minimum(qmin_POP, 10.**qmins[i]), Pbase, 10.**fEdds[j], xi, fbin, h, Om, OL))
		else:
		 print "Not sure what you doing here"	

	



		fig = plt.figure(figsize=[7.5,6.1])
		ax = fig.add_subplot(111)
		#plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$, $z_{\rm{max}}=%g$' %zmax, fontsize=15)
		plt.title(r'$\log_{10}{\left[ N_{\rm{VLBI}} \right]}$', fontsize=16)
		ax.contourf(qmins, fEdds, np.log10(Ntot_grid), 200, cmap = "viridis")

		lms = plt.contour(qmins, fEdds, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
	
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		

		# plt.axvline(x=0.0, color='black', linewidth=1, linestyle="--")


		plt.axvline(x=np.log10(qmin_EHT), color='chartreuse', linewidth=2, linestyle="--")
		plt.axhline(y=np.log10(f_Edd), color='chartreuse', linewidth=2, linestyle="--")
		plt.scatter(np.log10(qmin_EHT), np.log10(f_Edd), color='chartreuse', marker='o', s=30)



		ax.set_xlabel(r'$\rm{log}_{10}[q_{\rm{min}}]$')
		ax.set_ylabel(r'$\rm{log}_{10}[f_{\rm{Edd}}]$')



		FminSv = Fmin/mJy2cgs/1000.
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)

		plt.figtext(0.2,0.87, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
		plt.figtext(0.2,0.82, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		if (FminSv*1000.0 <= 100):
			FminSv = FminSv*1000.0
			plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		else:
			plt.figtext(0.2,0.77, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.72, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
		#plt.figtext(0.2,0.67, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)
		plt.figtext(0.2,0.67, r"$a_{\rm{max}}=10^{%g}$ pc" %np.log10(KQ),  color='yellow', fontsize=15)
		plt.figtext(0.2,0.62, r"$f_{\rm{bin}}=%g$" %fbin, color='yellow', fontsize=15)

		#plt.figtext(0.2,0.61, r"$\chi=%g$" %chi, color='yellow', fontsize=15)


		plt.ylim(fEdds[0], fEdds[len(epss)-1])
		plt.xlim(qmins[0], qmins[len(qmins)-1])

		plt.tight_layout()


		if (TrapInt):
			Savename = 'CumZ_fEdd_vs_qmin_%gx%g_Mdot%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_TrapInt%g.png'%(Ng,Ng, eps, zmax, FminSv, thMnSv, PbaseSv, Ntrap_z)
		else:
			Savename = 'CumZ_fEdd_vs_qmin_%gx%g_Mdot%g_zmax%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, eps, zmax, FminSv, thMnSv, PbaseSv, reclim)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)


























# if (CntPlt_DZ):
# 	Ng = 10
# 	#DZ = 0.2
# 	#zeval = 0.37  ##Max z count

# 	fEdds = np.linspace(-3.5, 0., Ng)
# 	epss = np.linspace(-3.5, 2.0, Ng)
# 	KQs = np.linspace(-3., 0.0, Ng)
# 	thMns = np.linspace(-1.0, 2.0, Ng) 
# 	MdEffs = np.linspace(-2., 0., Ng)
# 	Pbasez = np.linspace(0.0, np.log10(Pbase/yr2sec), Ng)
# 	Mmxz = np.linspace(5.0, 10.5, Ng)
# 	zs = np.linspace(-2.0, 0.0, Ng)
# 	qmins = np.linspace(-2., 0.0, Ng)






	# if (Mmx_Pbase):
	# 	print "Mmx vs Pbase"

	# 	Ntot_grid = np.zeros([Ng,Ng])
	# 	Ftot_grid = np.zeros([Ng,Ng])
	# 	RSGmx = np.zeros(Ng)
	# 	RSGmn = np.zeros(Ng)

	# 	for i in range(0,Ng):
	# 		for j in range(0,Ng):		
	# 			Ntot_grid[j][i] =  max(1.e-3, NtotDZ_GWgas(zeval, 10.**Mmxz[i], DZ, Fmin, chi, thMn, qmin, eps, f_Edd, 10.**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL))
	# 			#Ntot_grid[j][i] =  -IntzZ_OptNEHT([eps, KQ], zeval, 10.**Mmxz[i], Fmin, chi, thMn, qmin, 10.**Pbasez[j]*yr2sec, f_Edd, xi, fbin, h, Om, OL)
	
	# 	aseps = asep(10.**Pbasez[j]*yr2sec, 10.**Mmxz[i]*Msun)
		

	# 	plt.figure()
	# 	plt.title(r"$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g$" %(zeval), fontsize=18)
	# 	cnt = plt.contourf(Mmxz, Pbasez, np.log10(Ntot_grid), 200, cmap = "viridis")
	# 	#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
	# 	#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

	# 	cbar = plt.colorbar(cnt)
	# 	lms = plt.contour(Mmxz, Pbasez, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
	# 	#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
	# 	#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
	# 	plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	


	# 	plt.contour(Mmxz, Pbasez, np.log10(aseps), cmap = "viridis", levels = [-3.0, -2.0, -1.0, 0.0])
	# 	plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'blue', fontsize=14)



	# 	plt.ylabel(r'$\rm{log}_{10}[P_{\rm{base}}/yr]$')
	# 	plt.xlabel(r'$\rm{log}_{10}[M_{\rm{max}}/M_{\odot}]$')



	# 	FminSv = Fmin/mJy2cgs/1000.
	# 	thMnSv = thMn/mu_as2rad 
	# 	PbaseSv = Pbase/yr2sec


	# 	#plt.figtext(0.2,0.51, r"$\chi=%g$" %chi, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.46, r"$q_{\rm{min}}=%g$" %qmin, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.41, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.36, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.31, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.26, r"$a_{\rm{max}}=%g$ pc" %KQ, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.21, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)


	# 	plt.ylim(0.0, np.log10(Pbase/yr2sec))

	# 	plt.tight_layout()

	# 	Savename = 'dz_Pbase_vs_Mmx_%gx%g_amax%g_eps%g_qmin%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, KQ, eps, qmin_EHT, FminSv, thMnSv, PbaseSv, reclim)

	# 	Savename = Savename.replace('.', 'p')
	# 	Savename = Savename.replace('ppng', '.png')
	# 	plt.savefig(Savename)

















	# if (Mdot_Npc):

	# 	## CONTOUR PLOT TOTALs

	# 	Ntot_grid = np.zeros([Ng,Ng])
	# 	Ftot_grid = np.zeros([Ng,Ng])
	# 	RSGmx = np.zeros(Ng)
	# 	RSGmn = np.zeros(Ng)

	# 	aTmx = np.zeros(Ng)
	# 	aTmn = np.zeros(Ng)

	# 	for i in range(0,Ng):
	# 		for j in range(0,Ng):
	# 			Ntot_grid[j][i] =  max(1.e-3, NtotDZ_GWgas(zeval, Mmx, DZ, Fmin, chi, thMn, qmin, 10.**epss[i], f_Edd, Pbase, 10.**KQs[j], MdEff, xi, fbin, h, Om, OL) )
	# 			#if (TrapInt):
	# 			#Ntot_grid[j][i] =  -IntzZ_Trap_OptNEHT([10.**epss[i], 10.**KQs[j]], zeval, Mmx, Fmin, chi, thMn, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
	# 			#else:
	# 			#	Ntot_grid[j][i] =  -IntzZ_OptNEHT([10.**epss[i], 10.**KQs[j]], zeval, Mmx, Fmin, chi, thMn, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)

	# 			RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
	# 			RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm
	# 			aTmx[i] = aTrans(Mmax, 1.0, MdEff, 10.**epss[i])/pc2cm
	# 			aTmn[i] = aTrans(Mmin, qmin, MdEff, 10.**epss[i])/pc2cm


	# 	Log_thmxs = np.log10(10.**KQs*pc2cm/Dang(zeval, h, Om, OL) * 180./np.pi * 3600 * 1.e6)


	# 	fig = plt.figure()
	# 	ax = fig.add_subplot(111)
	# 	plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g$' %zeval, fontsize=15)
	# 	ax.contourf(epss, Log_thmxs, np.log10(Ntot_grid), 200, cmap = "viridis")
	# 	ax2 = ax.twinx()
	# 	ax2.contourf(epss, KQs, np.log10(Ntot_grid), 200, cmap = "viridis")
	# 	#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
	# 	#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

	# 	#cbar = plt.colorbar(cnt)
	# 	lms = plt.contour(epss, KQs, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
	# 	#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
	# 	#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
	# 	plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

	# 	plt.plot(epss, np.log10(RSGmx), color='red' )
	# 	plt.plot(epss, np.log10(RSGmn), color='red' )
	# 	plt.fill_between(epss,np.log10(RSGmx),np.log10(RSGmn), color='red', alpha=0.3)

	# 	plt.plot(epss, np.log10(aTmx), color='green', linewidth=2, linestyle="--" )
	# 	plt.plot(epss, np.log10(aTmn), color='green', linewidth=2, linestyle="--" )


	# 	plt.axvline(x=0.0, color='black', linewidth=2, linestyle="--")



	# 	ax.set_xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
	# 	#plt.ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
	# 	ax2.set_ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')
	# 	ax.set_ylabel(r'$\rm{log}_{10}[\theta_{\rm{max}}/\mu\rm{as}]$')


	# 	FminSv = Fmin/mJy2cgs/1000.
	# 	thMnSv = thMn/mu_as2rad 
	# 	PbaseSv = Pbase/yr2sec
	# 	# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
	# 	# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
	# 	# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
	# 	# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


	# 	plt.figtext(0.2,0.86, r"$q_{\rm{min}}=%g$" %qmin, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.81, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.76, r"$F_{\rm{min}}=%g$ Jy" %FminSv, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.71, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
	# 	plt.figtext(0.2,0.66, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)
	# 	#plt.figtext(0.2,0.61, r"$\chi=%g$" %chi, color='yellow', fontsize=15)



	# 	plt.ylim(-3.0, 0.0)

	# 	plt.tight_layout()

	# 	#plt.show()

	# 	if (TrapInt):
	# 		plt.savefig('dz_Npc_vs_Mdot_TrapInt%g_%gx%g_zeval%g_qmin%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ntrap_z, Ng,Ng,zeval, qmin, FminSv, thMnSv, PbaseSv, reclim))
	# 	else:
	# 		plt.savefig('dz_Npc_vs_Mdot_%gx%g_zeval%g_qmin%g_Fmin%gJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng,zeval, qmin_EHT, FminSv, thMnSv, PbaseSv, reclim))



