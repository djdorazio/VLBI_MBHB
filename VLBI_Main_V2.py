import numpy as np


import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt


import scipy.integrate as intg
import math as ma
import VLBI_IntFuncs_V2 as IFs 
from VLBI_IntFuncs_V2 import *


##INT options
TrapInt = True

###PLOTTING OPTIONS
DoDiags = False

##Contour Options
CntPlt = False

Opt_DZ_CntPlt = True	
Res_Cnts = False

DZdepPlt = False


DZ_CntPlt = False

Mmx_z = False
Pbase_z = False
Mmx_Pbase = False
Mdot_Npc = True

Mdot_qmin = False
Mdot_fEdd = False
Mdot_thMn = False
Mdot_Pbase = False
Mdot_Mdeff = False

## 1-D options
Odplot = False
PltMdot = True
Pltqmin = True
PltPbase = True









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
h=0.7
Om = 0.3
OL=0.7
alpha=-5./3.
beta=8./3.
gamma = -1.
Mmx = 10.**12 ## jsut to not limit lum function - doesnt change anyhting when set at 2*10^10 number
Mmax = 2.*10.**10*Msun
Mmin= 10.**5*Msun 
thMx = 100.*mu_as2rad
PNyq = 0.01*yr2sec




###LUM FINC PARAMS
#FineShanks
#CC = ma.exp(-7.24)
#Lstr = 10.**33
#A = 3.71
#B = -5.
zlim = 1.95

###LUM FINC PARAMS  - grabage
CC = ma.exp(-6.05)#10.**(-6.05)
CCQM = ma.exp(-8.04)#10.**(-8.04)
KK = 0.06 
KKQM = 2.93
A = -0.31
AQM = -0.32
B = -1.88
BQM = -1.75
Lstr = ma.exp(25.17+7.)#10.**(32.17)
LstrQM = ma.exp(26.96+7.)#10.**(33.96)
chi  = 5.0#0.1  #chi L_14 = L_mm (specifc fluxes)
tQuasar = 10.**8 * yr2sec
xii = 1.0


zmax = 6.0#zmax_int(thMn, Mmax, h, Om, OL, Pbase)

##FREE (and KEY) PARAMETERS
fbin = 1.0
KQ = 0.01 ## sets number of pc at which RL turns on
eps = 0.001#10**(-3.75)  ## sets migration (accretion rate in CBD pushing two together)
f_Edd = 0.1  ## sets L to M connection (accretion rate onto shining BH)
MdEff = 0.1
xi = 1.0
qmin = 0.01

zeval = 0.05  #best for Pbase=10 yrs
DZ = 1.0#1.0/(4.*ma.pi)
Fmin = 10.0 * mJy2cgs
maglim = 24.5
##temp
nuVbnd = c/(5.45*10**(-5))
F0opt  = 3.636*10**(-20)*nuVbnd 
#maglim = -2.5*ma.log10(Fmin/F0opt)
Fmin_opt = 10.**(-maglim/2.5)*F0opt/nuVbnd 
thMn = 1.0 * mu_as2rad 
Pbase = 10.0*yr2sec

#Mmax = min(Mmax,  4.*ma.pi**2 * (KQ*pc2cm)**3/(G*(Pbase**2)  ))  ##tech z dep here
##Limit to Mmax for which Pbase < PRLQ (set by KQ*pc2cm)

# #ONLY if fixed upper lim on P
# PQ = 2. * ma.pi * (pc2cm)**(3./2.) /ma.sqrt(G*10.**8*Msun)


if (DZdepPlt):
	Ng = 1000
	Ftot_DZ1 = np.zeros(Ng)
	Ftot_DZ2 = np.zeros(Ng)
	Ftot_DZ3 = np.zeros(Ng)
	DZ1 = 0.1
	DZ2 = 0.2
	DZ3 = 0.3 
	zzeval = np.linspace(0.01, 2.0, Ng) 
	for i in range(0,Ng):
			Ftot_DZ1[i] =  NtotDZ_GWgas(float(zzeval[i]), DZ1, Fmin, chi, thMn, qmin, eps, f_Edd, Pbase/2., KQ, MdEff, xi, fbin, h, Om, OL)
			Ftot_DZ2[i] =  NtotDZ_GWgas(float(zzeval[i]), DZ1, Fmin, chi, thMn, qmin, eps, f_Edd, Pbase*1., KQ, MdEff, xi, fbin, h, Om, OL)			
			Ftot_DZ3[i] =  NtotDZ_GWgas(float(zzeval[i]), DZ1, Fmin, chi, thMn, qmin, eps, f_Edd, Pbase*2., KQ, MdEff, xi, fbin, h, Om, OL)

	zmax_i = np.where(Ftot_DZ2 == max(Ftot_DZ2))[0][0]
	zmx = zzeval[zmax_i]
	print "Peak z = %g" %zmx

	plt.plot(zzeval, np.log10(Ftot_DZ1), color="blue", linewidth=3)
	plt.plot(zzeval, np.log10(Ftot_DZ2), color="black", linewidth=3)
	plt.plot(zzeval, np.log10(Ftot_DZ3), color="red", linewidth=3)

	# plt.plot(zzeval, Ftot_DZ1, color="blue", linewidth=3)
	# plt.plot(zzeval, Ftot_DZ2, color="black", linewidth=3)
	# plt.plot(zzeval, Ftot_DZ3, color="red", linewidth=3)


	plt.xlabel(r'$z$')
	plt.ylabel(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$ [All Sky]')

	plt.tight_layout()

	#plt.show()
	plt.savefig('FtotDZ_vs_z_N%g_qmin%g.png'%(Ng, qmin))















if (Res_Cnts):
	Ng = 500

	#DZ = 0.2
	#zeval = 0.37
	qstst = qsofq(0.01)


	Ps = np.linspace(-1., 1., Ng)
	Ms = np.linspace(6., 10., Ng)

	Npmq = np.zeros([Ng,Ng])
	DOPmx = np.zeros([Ng,Ng])
	hPTAs = np.zeros([Ng,Ng])


	## TOTAL NUMBER FROM sMALF
	dV = dVdzdOm(zeval, h, Om, OL)
	FInt_LF = FInt_smLF(zeval, chi, Fmin, h, Om, OL)
	Ntotz   = zeval*4.*ma.pi*dV*FInt_LF/(10.**6 * pc2cm)**3

	for i in range(0,Ng):
		for j in range(0,Ng):
			#Npmq[j][i] = Ms[i]
			Npmq[j][i]  = FbinPM(10.**Ps[j]*yr2sec, 10.**Ms[i]*Msun, qstst, MdEff, eps, xi) 
			DOPmx[j][i] = DopMax(10.**Ps[j]*yr2sec, 10.**Ms[i]*Msun, qstst, 1.0)
			hPTAs[j][i] = hPTA(10.**Ps[j]*yr2sec, 10.**Ms[i]*Msun, qstst,zeval,h, Om, OL)

	
	Npmq = 	Npmq*Ntotz
	PmnDN   = PminRes(10**Ms*Msun, thMn, zeval-DZ, h, Om, OL)
	PbaseDN = PbaseObs(Pbase,zeval-DZ)*Ms/Ms

	Pmns   = PminRes(10**Ms*Msun, thMn, zeval, h, Om, OL)
	Pbases = PbaseObs(Pbase,zeval)*Ms/Ms

	PmnUP   = PminRes(10**Ms*Msun, thMn, zeval+DZ, h, Om, OL)
	PbaseUP = PbaseObs(Pbase,zeval+DZ)*Ms/Ms



	


	matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
	#matplotlib.rcParams['contour.positive_linestyle'] = 'dashed'
	plt.figure()
	plt.title(r"$\log_{10}{\left[ N_{\rm{EHT}}(P,M,q) \right]}$, $z=%g\pm%g$, $q=0.3$" %(zeval, DZ), fontsize=15)
	### tres contours
	cnt = plt.contourf(Ms, Ps, np.log10(Npmq), 200, cmap = "viridis")
	cbar = plt.colorbar(cnt)
	lms = plt.contour(Ms, Ps, np.log10(Npmq), cmap = "viridis", levels = [0.0, 1., 2., 3., 4., 5., 6.])
	#lms = plt.contour(Ms, Ps, np.log10(Npmq), cmap = "viridis", levels = [-7.0, -6.0, -5.0, -4.0,-3.0,-2.0, -1.0, 0.0])
	#plt.clabel(lms, fmt = '%2.1d', colors = 'k', fontsize=14)
	plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	


	

	### DOP contours
	Dops = plt.contour(Ms, Ps, DOPmx, colors='blue', levels = [1.1, 1.5, 2.], linestyles = 'dashed', linewidths=2)
	plt.clabel(Dops, fmt = '%2.2f', colors = 'blue', fontsize=14)	

	### hPTA contours
	hs = plt.contour(Ms, Ps, np.log10(hPTAs), colors='yellow', levels = [-17., -16.],linewidths=2, linestyles='-')
	plt.clabel(hs, fmt = '%2.2f', colors = 'yellow', fontsize=14)	


### Pmins and Pmaxes
	plt.plot(Ms, np.log10(PmnDN/yr2sec), color='black', linewidth=1, linestyle="--")
	plt.plot(Ms, np.log10(PbaseDN/yr2sec), color='black', linewidth=1, linestyle="--")

	plt.plot(Ms,  np.log10(Pmns/yr2sec), color='black', linewidth=2)
	plt.plot(Ms, np.log10(Pbases/yr2sec), color='black', linewidth=2)

	plt.plot(Ms,  np.log10(PmnUP/yr2sec), color='black',  linewidth=1, linestyle="--")
	plt.plot(Ms, np.log10(PbaseUP/yr2sec), color='black', linewidth=1, linestyle="--")

	plt.fill_between(Ms, -Ms/Ms, np.log10(PmnDN/yr2sec), color='gray', alpha=0.7)
	plt.fill_between(Ms, Ms/Ms, np.log10(PbaseDN/yr2sec), color='gray', alpha=0.7)


	plt.fill_between(Ms, np.log10(PmnDN/yr2sec), np.log10(PmnUP/yr2sec), color='gray', alpha=0.25)
	plt.fill_between(Ms, np.log10(PbaseDN/yr2sec), np.log10(PbaseUP/yr2sec), color='gray', alpha=0.25)


	plt.figtext(0.2,0.65, r"-$P_{\rm{min}}, \frac{ P_{\rm{base}} }{1+z}$", color='black', fontsize=15)
	plt.figtext(0.2,0.59, r"--$D^{3-\alpha}$", color='blue', fontsize=15)
	plt.figtext(0.2,0.55, r"- $\log_{10}{\left[h_{\rm{PTA}}\right]}$", color='yellow', fontsize=15)



	plt.ylim(-1.,1.)
	plt.xlim(6.,10.)

	plt.ylabel(r'$\rm{log}_{10}[P/\rm{yr}]$')
	plt.xlabel(r'$\rm{log}_{10}[M/M_{\odot}]$')

	plt.tight_layout()

	#plt.show()
	FminSv = Fmin/mJy2cgs
	thMnSv = thMn/mu_as2rad 
	PbaseSv = Pbase/yr2sec
	plt.savefig('NPMQ_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr.png'%(Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv))

























if (Opt_DZ_CntPlt):
	Ng = 10
	#DZ = 0.2
	#zeval = 0.37  ##Max z count

	fEdds = np.linspace(-4.0, 0., Ng)
	epss = np.linspace(-3.5, 1., Ng)
	KQs = np.linspace(-3., 0.0, Ng)
	thMns = np.linspace(-1.0, 2.0, Ng) 
	MdEffs = np.linspace(-2., 0., Ng)
	Pbasez = np.linspace(0., 2., Ng)
	zs = np.linspace(0.01, 1., Ng)
	qmins = np.linspace(-2., 0.0, Ng)
	xis   = np.linspace(-2., 2., Ng)




	if (Mdot_Npc):

		## CONTOUR PLOT TOTALs


		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)
		# RSGmx2 = np.zeros([Ng,Ng])
		# RSGmn2 = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):
				#Ftot_grid[j][i] =  FtotDZ(DZ, zeval, h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, 10.**KQs[j], MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)	
				Ntot_grid[j][i] =  max(1.e-6, NtotDZ_OPTICAL_GWgas(zeval, DZ, Fmin_opt, chi, thMn, qmin, 10.**epss[i], 10.**epss[i], Pbase, 10.**KQs[j], MdEff, xi, fbin, h, Om, OL) )
				#Ftot_grid[j][i] =  max(1.e-14, NtotDZ_GWgas(zeval, DZ, Fmin, chi, thMn, qmin, 10.**epss[i], f_Edd, Pbase, 10.**KQs[j], MdEff, xi, fbin, h, Om, OL)/NtotDZ_RLF(zeval, DZ, Fmin, chi, h, Om, OL) )

				#RSGmx2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmin, MdEff)/pc2cm
				RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm


		plt.figure()
		plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g\pm%g$' %(zeval, DZ), fontsize=15)
		cnt = plt.contourf(epss, KQs, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, KQs, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		plt.plot(epss, np.log10(RSGmx), color='red' )
		plt.plot(epss, np.log10(RSGmn), color='red' )
		plt.fill_between(epss,np.log10(RSGmx),np.log10(RSGmn), color='red', alpha=0.3)

		# plt.contour(epss, KQs, np.log10(RSGmx2), cmap = "Reds")
		# plt.contour(epss, KQs, np.log10(RSGmn2), cmap = "Reds")


		plt.axvline(x=0.0, color='black', linewidth=2, linestyle="--")



		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		plt.ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')


		FminSv = Fmin_opt/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		plt.figtext(0.2,0.86, r"$q_{\rm{min}}=%g$" %qmin, color='yellow', fontsize=15)
		plt.figtext(0.2,0.81, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.76, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.71, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)


		plt.ylim(-3.0, 0.0)

		plt.tight_layout()

		#plt.show()

		plt.savefig('OPTICAL_Npc_vs_Mdot_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr.png'%(Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv))
























if (DZ_CntPlt):
	Ng = 20
	#DZ = 0.2
	#zeval = 0.37  ##Max z count

	fEdds = np.linspace(-3.5, 0., Ng)
	epss = np.linspace(-3.5, 1.0, Ng)
	KQs = np.linspace(-3., 0.0, Ng)
	thMns = np.linspace(-1.0, 2.0, Ng) 
	MdEffs = np.linspace(-2., 0., Ng)
	Pbasez = np.linspace(0.0, np.log10(Pbase/yr2sec), Ng)
	Mmxz = np.linspace(5.0, 10.5, Ng)
	zs = np.linspace(-2.0, 0.0, Ng)
	qmins = np.linspace(-2., 0.0, Ng)
	xis   = np.linspace(-2., 2., Ng)





	if (Mmx_z):
		print "Mmx vs z"


		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])

		MRBs = np.log10(MresBase(Pbase, thMn, zs, h, Om, OL)/Msun)

		for i in range(0,Ng):
			for j in range(0,Ng):		
				Ntot_grid[j][i] =  max(1.e-3, NtotDZ_GWgas(10.**zs[i], 10.**Mmxz[j], DZ, Fmin, chi, thMn, qmin, eps, eps, Pbase, KQ, MdEff, xi, fbin, h, Om, OL))

		
		plt.figure()
		#plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g\pm%g$' %(zeval, DZ), fontsize=15)
		cnt = plt.contourf(zs, Mmxz, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		cbar = plt.colorbar(cnt)
		lms = plt.contour(zs, Mmxz, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		plt.plot(zs,MRBs, color='black', linewidth=3)


		plt.ylabel(r'$\rm{log}_{10}[M_{\rm{max}}/M_{\odot}]$')
		plt.xlabel(r'$z$')


		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		#KQsv = #10.**KQ
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		plt.figtext(0.2,0.41, r"$q_{\rm{min}}=%g$" %qmin, color='yellow', fontsize=15)
		plt.figtext(0.2,0.36, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.31, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.26, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
		plt.figtext(0.2,0.21, r"$a_{\rm{max}}=%g$ pc" %KQ, color='yellow', fontsize=15)
		plt.figtext(0.2,0.16, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)


		#plt.ylim(0.0, np.log10(Pbase/yr2sec))

		plt.tight_layout()

		#plt.show()

		plt.savefig('Mmx_vs_z_%gx%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, FminSv, thMnSv, PbaseSv, reclim))















	if (Pbase_z):
		print "Pbase vs z"


		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)
		# RSGmx2 = np.zeros([Ng,Ng])
		# RSGmn2 = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):		
				Ntot_grid[j][i] =  max(1.e-3, NtotDZ_GWgas(10**zs[i], Mmx, DZ, Fmin, chi, thMn, qmin, eps, eps, 10.**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL))

				#RSGmx2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmin, MdEff)/pc2cm
				#RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm

		Pmin_M6 = PminRes(10.**6 * Msun, thMn, 10.**zs, h, Om, OL)/yr2sec * (1.+10.**zs)
		Pmin_M7 = PminRes(10.**7 * Msun, thMn, 10.**zs, h, Om, OL)/yr2sec * (1.+10.**zs)
		Pmin_M8 = PminRes(10.**8 * Msun, thMn, 10.**zs, h, Om, OL)/yr2sec * (1.+10.**zs)
		Pmin_M9 = PminRes(10.**9 * Msun, thMn, 10.**zs, h, Om, OL)/yr2sec * (1.+10.**zs)
		Pmin_M10 = PminRes(10.**10 * Msun, thMn, 10.**zs, h, Om, OL)/yr2sec * (1.+10.**zs)

		plt.figure()
		#plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g\pm%g$' %(zeval, DZ), fontsize=15)
		cnt = plt.contourf(zs, Pbasez, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		cbar = plt.colorbar(cnt)
		lms = plt.contour(zs, Pbasez, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		plt.plot(zs, np.log10(Pmin_M6), color='red', linestyle='--')
		plt.plot(zs, np.log10(Pmin_M7), color='orange', linestyle='--')
		plt.plot(zs, np.log10(Pmin_M8), color='blue', linestyle='--' )
		plt.plot(zs, np.log10(Pmin_M9), color='green', linestyle='--' )
		plt.plot(zs, np.log10(Pmin_M10), color='black', linestyle='--' )
		#plt.fill_between(zs,np.log10(Pmin_M6),np.log10(Pmin_M10), color='red', alpha=0.3)

		# plt.contour(epss, KQs, np.log10(RSGmx2), cmap = "Reds")
		# plt.contour(epss, KQs, np.log10(RSGmn2), cmap = "Reds")


		#plt.axvline(x=0.0, color='black', linewidth=2, linestyle="--")



		plt.ylabel(r'$\rm{log}_{10}[P_{\rm{base}}/yr]$')
		plt.xlabel(r'$z$')


		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		#KQsv = #10.**KQ
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		plt.figtext(0.2,0.41, r"$q_{\rm{min}}=%g$" %qmin, color='yellow', fontsize=15)
		plt.figtext(0.2,0.36, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.31, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.26, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
		plt.figtext(0.2,0.21, r"$a_{\rm{max}}=%g$ pc" %KQ, color='yellow', fontsize=15)
		#plt.figtext(0.2,0.71, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)


		plt.ylim(0.0, np.log10(Pbase/yr2sec))

		plt.tight_layout()

		#plt.show()

		plt.savefig('Pbase_vs_z_%gx%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, FminSv, thMnSv, PbaseSv, reclim))












	if (Mmx_Pbase):
		print "Mmx vs Pbase"

		zeval = 0.05
		f_Edd = 0.1
		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)
		# RSGmx2 = np.zeros([Ng,Ng])
		# RSGmn2 = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):		
				#Ntot_grid[j][i] =  max(1.e-3, NtotDZ_GWgas(zeval, 10.**Mmxz[i], DZ, Fmin, chi, thMn, qmin, eps, 0.1, 10.**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL))
				Ntot_grid[j][i] =  -IntzZ_OptNEHT([eps, KQ], zeval, 10.**Mmxz[i], Fmin, chi, thMn, qmin, 10.**Pbasez[j]*yr2sec, f_Edd, xi, fbin, h, Om, OL)
				#RSGmx2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmin, MdEff)/pc2cm
				#RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm

		
		#plt.figure(figsize=[8,7.5])
		plt.figure()
		plt.title(r"$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g$" %(zeval), fontsize=18)
		cnt = plt.contourf(Mmxz, Pbasez, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		cbar = plt.colorbar(cnt)
		lms = plt.contour(Mmxz, Pbasez, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	



		plt.ylabel(r'$\rm{log}_{10}[P_{\rm{base}}/yr]$')
		plt.xlabel(r'$\rm{log}_{10}[M_{\rm{max}}/M_{\odot}]$')



		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		#KQsv = #10.**KQ
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		plt.figtext(0.2,0.46, r"$q_{\rm{min}}=%g$" %qmin, color='yellow', fontsize=15)
		plt.figtext(0.2,0.41, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.36, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.31, r"$\dot{\mathcal{M}}=%g$" %eps, color='yellow', fontsize=15)
		plt.figtext(0.2,0.26, r"$a_{\rm{max}}=%g$ pc" %KQ, color='yellow', fontsize=15)
		plt.figtext(0.2,0.21, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)


		plt.ylim(0.0, np.log10(Pbase/yr2sec))

		plt.tight_layout()

		#plt.show()

		plt.savefig('Pbase_vs_Mmx_%gx%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng, FminSv, thMnSv, PbaseSv, reclim))



















	if (Mdot_Npc):

		## CONTOUR PLOT TOTALs

		zeval = 0.1
		Ntot_grid = np.zeros([Ng,Ng])
		Ftot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)
		# RSGmx2 = np.zeros([Ng,Ng])
		# RSGmn2 = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):
				#Ftot_grid[j][i] =  FtotDZ(DZ, zeval, h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, 10.**KQs[j], MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)	
				#Ntot_grid[j][i] =  max(1.e-6, NtotDZ_GWgas(zeval, Mmx, DZ, Fmin, chi, thMn, qmin, 10.**epss[i], f_Edd, Pbase, 10.**KQs[j], MdEff, xi, fbin, h, Om, OL) )
				
				#Ntot_grid[j][i] =  -IntzZ_OptNEHT([10.**epss[i], 10.**KQs[j]], zeval, Mmx, Fmin, chi, thMn, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
				#if (TrapInt):
				Ntot_grid[j][i] =  -IntzZ_Trap_OptNEHT([10.**epss[i], 10.**KQs[j]], zeval, Mmx, Fmin, chi, thMn, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)
				#else:
				#	Ntot_grid[j][i] =  -IntzZ_OptNEHT([10.**epss[i], 10.**KQs[j]], zeval, Mmx, Fmin, chi, thMn, qmin, Pbase, f_Edd, xi, fbin, h, Om, OL)


				#Ftot_grid[j][i] =  max(1.e-14, NtotDZ_GWgas(zeval, DZ, Fmin, chi, thMn, qmin, 10.**epss[i], f_Edd, Pbase, 10.**KQs[j], MdEff, xi, fbin, h, Om, OL)/NtotDZ_RLF(zeval, DZ, Fmin, chi, h, Om, OL) )

				#RSGmx2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmin, MdEff)/pc2cm
				RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm


		plt.figure()
		plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g$' %zeval, fontsize=15)
		cnt = plt.contourf(epss, KQs, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, KQs, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, np.log10(Ftot_grid), cmap = "viridis", levels = [-6.0,-5.0,-4.0, -3.0, -2.0, -1.0, 0.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		plt.plot(epss, np.log10(RSGmx), color='red' )
		plt.plot(epss, np.log10(RSGmn), color='red' )
		plt.fill_between(epss,np.log10(RSGmx),np.log10(RSGmn), color='red', alpha=0.3)

		# plt.contour(epss, KQs, np.log10(RSGmx2), cmap = "Reds")
		# plt.contour(epss, KQs, np.log10(RSGmn2), cmap = "Reds")


		plt.axvline(x=0.0, color='black', linewidth=2, linestyle="--")



		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		plt.ylabel(r'$\rm{log}_{10}[a_{\rm{max}}/\rm{pc}]$')


		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		plt.figtext(0.2,0.86, r"$q_{\rm{min}}=%g$" %qmin, color='yellow', fontsize=15)
		plt.figtext(0.2,0.81, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.76, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.71, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='yellow', fontsize=15)
		plt.figtext(0.2,0.66, r"$f_{\rm{Edd}}=%g$" %f_Edd, color='yellow', fontsize=15)


		plt.ylim(-3.0, 0.0)

		plt.tight_layout()

		#plt.show()

		if (TrapInt):
			plt.savefig('Npc_vs_Mdot_TrapInt%g_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ntrap_z, Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv, reclim))
		else:
			plt.savefig('Npc_vs_Mdot_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr_reclim%g.png'%(Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv, reclim))



	if (Mdot_fEdd):

		## CONTOUR PLOT TOTALs

		zeval = 0.05
		Ntot_grid = np.zeros([Ng,Ng])
		RSGmx = np.zeros(Ng)
		RSGmn = np.zeros(Ng)
		# RSGmx2 = np.zeros([Ng,Ng])
		# RSGmn2 = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):
				#Ftot_grid[j][i] =  FtotDZ(DZ, zeval, h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, 10.**KQs[j], MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)	
				#Ftot_grid[j][i] =  max(1.e-6, NtotDZ_GWgas(zeval, Mmx, DZ, Fmin, chi, thMn, qmin, 10.**epss[i], 10**fEdds[j], Pbase, KQ, MdEff, xi, fbin, h, Om, OL) )
				Ntot_grid[j][i] =  -IntzZ_OptNEHT([10.**epss[i], KQ], zeval, Mmx, Fmin, chi, thMn, qmin, Pbase, 10.**fEdds[j], xi, fbin, h, Om, OL)


				#RSGmx2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn2[j][i] = RSG2(10.**KQs[j], 10.**epss[i], Mmin, MdEff)/pc2cm
				#RSGmx[i] = RSG(10.**epss[i], Mmax, MdEff)/pc2cm
				#RSGmn[i] = RSG(10.**epss[i], Mmin, MdEff)/pc2cm


		plt.figure()
		plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g$' %(zeval, DZ), fontsize=15)
		cnt = plt.contourf(epss, fEdds, np.log10(Ntot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, fEdds, np.log10(Ntot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	

		


		plt.axvline(x=0.0, color='black', linewidth=2, linestyle="--")



		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		plt.ylabel(r'$\rm{log}_{10}[f_{\rm{Edd}}]$')


		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		# plt.figtext(0.2,0.36, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		# plt.figtext(0.2,0.31, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		# plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		plt.figtext(0.6,0.41, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		plt.figtext(0.6,0.36, r"$\theta_{\rm{min}}=%g \mu$as" %thMnSv, color='black', fontsize=15)
		plt.figtext(0.6,0.31, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		plt.figtext(0.6,0.26, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)
		plt.figtext(0.2,0.21, r"$a_{\rm{max}}=%g$ pc" %KQ, color='yellow', fontsize=15)


		#plt.ylim(-3.0, 0.0)

		plt.tight_layout()

		#plt.show()

		plt.savefig('fEdd_vs_Mdot_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr.png'%(Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv))






	if (Mdot_thMn):

		## CONTOUR PLOT TOTALs


		Ftot_grid = np.zeros([Ng,Ng])


		for i in range(0,Ng):
			for j in range(0,Ng):
				#Ftot_grid[j][i] =  FtotDZ(DZ, zeval, h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, 10.**KQs[j], MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)	
				Ftot_grid[j][i] =  max(1.e-14, NtotDZ_GWgas(zeval, Mmx, DZ, Fmin, chi, 10.**thMns[j]* mu_as2rad , qmin, 10.**epss[i], f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL))
				


		plt.figure()
		plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g\pm%g$' %(zeval, DZ), fontsize=15)
		cnt = plt.contourf(epss, thMns, np.log10(Ftot_grid), 200, cmap = "viridis")
		#cnt = plt.contourf(epss, KQs, Ftot_grid, 200, cmap = "viridis")

		cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, thMns, np.log10(Ftot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		#lms = plt.contour(epss, KQs, Ftot_grid, cmap = "viridis", levels = [1.0, 10.0, 100.0, 1.e3, 1.e4, 1.e5])
		plt.clabel(lms, fmt = r'$10^{%g}$', colors = 'k', fontsize=14)	



		plt.axvline(x=0.0, color='black', linewidth=2, linestyle="--")



		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		plt.ylabel(r'$\rm{log}_{10}[\theta_{\rm{min}}/\mu as]$')


		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		plt.figtext(0.2,0.31, r"$q_{\rm{min}}=%g$" %qmin, color='black', fontsize=15)
		plt.figtext(0.2,0.26, r"$F_{\rm{min}}=%g$ mJy" %FminSv, color='black', fontsize=15)
		plt.figtext(0.2,0.21, r"$P_{\rm{base}}=%g$ yr" %PbaseSv, color='black', fontsize=15)


		#plt.ylim(-3.0, 1.0)

		plt.tight_layout()

		#plt.show()

		plt.savefig('thMinn_vs_Mdot_%gx%g_zeval%g_DZ%g_Fmin%gmJy_KQ%g_Pbase%gyr.png'%(Ng,Ng,zeval, DZ, FminSv, KQ, PbaseSv))








	if (Mdot_Mdeff):

		## CONTOUR PLOT TOTALs


		Ftot_grid = np.zeros([Ng,Ng])
		#Args_gr = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):
				#Ftot_grid[j][i] =  FtotDZ(DZ, zeval, h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, 10.**KQs[j], MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)	
				Ftot_grid[j][i] =  max(1.e-14, NtotDZ_GWgas(zeval, Mmx, DZ, Fmin, chi, thMn, qmin, 10.**epss[i], f_Edd, Pbase, KQ, 10**MdEffs[j], xi, fbin, h, Om, OL) )
		
		plt.figure()
		plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g\pm%g$' %(zeval, DZ), fontsize=16)
		cnt = plt.contourf(epss, MdEffs, np.log10(Ftot_grid), 200, cmap = "viridis")
		cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, MdEffs, np.log10(Ftot_grid), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		plt.clabel(lms, fmt = '%2.1d', colors = 'k', fontsize=14)	

		plt.axvline(x=0.0, color='black', linewidth=2, linestyle="--")

		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		plt.ylabel(r'$\rm{log}_{10}[\eta]$')

		plt.tight_layout()

		#plt.show()
		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		plt.savefig('MdEff_vs_Mdot_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr.png'%(Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv))


	if (Mdot_qmin):
		Ftot_grid_qMd = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):
				#Ftot_grid_qMd[j][i] =  FtotDZ(DZ, zeval, h, Om, OL, Mmax, Mmin, thMn, 10.**qmins[j], 10.**epss[i], Pbase, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)	
				Ftot_grid_qMd[j][i] =  max(1.e-14, NtotDZ_GWgas(zeval, Mmx, DZ, Fmin, chi, thMn, 10.**qmins[j], 10.**epss[i], Pbase, f_Edd, KQ, MdEff, xi, fbin, h, Om, OL) )


		plt.figure()
		plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g\pm%g$' %(zeval, DZ), fontsize=16)
		cnt = plt.contourf(epss, qmins, np.log10(Ftot_grid_qMd), 200, cmap = "viridis")
		cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, qmins, np.log10(Ftot_grid_qMd), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
		plt.clabel(lms, fmt = '%2.1d', colors = 'k', fontsize=14)


		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		plt.ylabel(r'$\rm{log}_{10}[q_{\rm{min}}]$')

		plt.tight_layout()

		#plt.show()
		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		plt.savefig('qmin_vs_Mdot_qmin_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr.png'%(Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv))



	if (Mdot_Pbase):
		Ftot_grid_P = np.zeros([Ng,Ng])

		for i in range(0,Ng):
			for j in range(0,Ng):
				#Ftot_grid_P[j][i] =  FtotDZ(DZ, zeval, h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], 10.**Pbasez[j]*yr2sec, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)	
				Ftot_grid_P[j][i] =  max(1.e-14, NtotDZ_GWgas(zeval, Mmx, DZ, Fmin, chi, thMn, qmin, 10.**epss[i], f_Edd, 10.**Pbasez[j]*yr2sec, KQ, MdEff, xi, fbin, h, Om, OL) )



		plt.figure()
		plt.title(r'$\log_{10}{\left[ N_{\rm{EHT}} \right]}$, $z=%g\pm%g$' %(zeval, DZ), fontsize=16)
		cnt = plt.contourf(epss, Pbasez, np.log10(Ftot_grid_P), 200, cmap = "viridis")
		cbar = plt.colorbar(cnt)
		lms = plt.contour(epss, Pbasez, np.log10(Ftot_grid_P), cmap = "viridis", levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
		plt.clabel(lms, fmt = '%2.1d', colors = 'k', fontsize=14)


		plt.xlabel(r'$\rm{log}_{10}[\dot{\mathcal{M}}]$')
		plt.ylabel(r'$\rm{log}_{10}[P_{\rm{base}}]$')

		plt.tight_layout()

		#plt.show()
		FminSv = Fmin/mJy2cgs
		thMnSv = thMn/mu_as2rad 
		PbaseSv = Pbase/yr2sec
		plt.savefig('Pbase_vs_Mdot_qmin_%gx%g_zeval%g_DZ%g_Fmin%gmJy_thMn%gmuas_Pbase%gyr.png'%(Ng,Ng,zeval, DZ, FminSv, thMnSv, PbaseSv))





















































if (Odplot):
	Nl = 10
	#epss = np.linspace(-6., 0., Nl)
	epss = np.linspace(-4., 0., Nl)
	if (PltMdot):
		Ftot_grid_Md1 = np.zeros(Nl)
		Ftot_grid_Md2 = np.zeros(Nl)
		Ftot_grid_Md3 = np.zeros(Nl)

		for i in range(0,Nl):
			Args_gr1 = (h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, KQ, 0.1, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_Md1[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr1)[0]	

			Args_gr2 = (h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, KQ, 0.2, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_Md2[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr2)[0]	

			Args_gr3 = (h, Om, OL, Mmax, Mmin, thMn, qmin, 10.**epss[i], Pbase, PQ, KQ, 0.3, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_Md3[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr3)[0]	





		plt.figure()

		cnt1 = plt.plot(epss, np.log10(Ftot_grid_Md1), color='black', linewidth=3, linestyle="-")
		cnt2 = plt.plot(epss, np.log10(Ftot_grid_Md2), color='blue', linewidth=3, linestyle="--")
		cnt3 = plt.plot(epss, np.log10(Ftot_grid_Md3), color='red', linewidth=3, linestyle=":")



		plt.xlabel(r'$\rm{log}_{10}\left[\dot{\mathcal{M}}\right]$')
		plt.ylabel(r'$\rm{log}_{10}\left[N_{\rm{tot}}\right]$ (All Sky)')

		plt.legend( [ cnt1[0], cnt2[0], cnt3[0] ], (r'$\eta=0.1$', r'$\eta=0.2$', r'$\eta=0.3$'), loc='upper right', fontsize=18)


		plt.tight_layout()

		#plt.show()
		plt.savefig('Ftot_vs_Mdot_N%g_zmax%g.png'%(Nl,zmax))





	if (Pltqmin):
		qmins = np.linspace(-2, 0., Nl)#np.linspace(-3., ma.log10(0.9), Nl)
		eps1 = 0.1
		eps2 = 0.01
		eps3 = 0.001

		Ftot_grid_q1 = np.zeros(Nl)
		Ftot_grid_q2 = np.zeros(Nl)
		Ftot_grid_q3 = np.zeros(Nl)

		for i in range(0,Nl):
			Args_gr1 = (h, Om, OL, Mmax, Mmin, thMn, 10.**qmins[i], eps1, Pbase, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_q1[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr1)[0]	

			Args_gr2 = (h, Om, OL, Mmax, Mmin, thMn, 10.**qmins[i], eps2, Pbase, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_q2[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr2)[0]	

			Args_gr3 = (h, Om, OL, Mmax, Mmin, thMn, 10.**qmins[i], eps3, Pbase, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_q3[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr3)[0]	





		plt.figure()

		cnt1 = plt.plot(qmins, np.log10(Ftot_grid_q1), color='black', linewidth=3, linestyle="-")
		cnt2 = plt.plot(qmins, np.log10(Ftot_grid_q2), color='blue', linewidth=3, linestyle="--")
		cnt3 = plt.plot(qmins, np.log10(Ftot_grid_q3), color='red', linewidth=3, linestyle=":")



		plt.xlabel(r'$\rm{log}_{10}\left[q_{\rm{min}}\right]$')
		plt.ylabel(r'$\rm{log}_{10}\left[N_{\rm{tot}}\right]$ (All Sky)')

		plt.legend( [ cnt1[0], cnt2[0], cnt3[0] ], (r'$\rm{log}_{10}\dot{\mathcal{M}}=-1$', r'$\rm{log}_{10}\dot{\mathcal{M}}=-2$', r'$\rm{log}_{10}\dot{\mathcal{M}}=-3$'), loc='upper right', fontsize=18)



		plt.tight_layout()

		#plt.show()
		plt.savefig('Ftot_vs_qmin_N%g_zmax%g.png'%(Nl,zmax))





	if (PltPbase):
		Pbasez = np.linspace(-1, 2., Nl)
		eps1 = 0.1
		eps2 = 0.01
		eps3 = 0.001

		Ftot_grid_P1 = np.zeros(Nl)
		Ftot_grid_P2 = np.zeros(Nl)
		Ftot_grid_P3 = np.zeros(Nl)

		for i in range(0,Nl):
			Args_gr1 = (h, Om, OL, Mmax, Mmin, thMn, qmin, eps1, 10.**Pbasez[i]*yr2sec, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_P1[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr1)[0]	

			Args_gr2 = (h, Om, OL, Mmax, Mmin, thMn, qmin, eps2, 10.**Pbasez[i]*yr2sec, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_P2[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr2)[0]	

			Args_gr3 = (h, Om, OL, Mmax, Mmin, thMn, qmin, eps3, 10.**Pbasez[i]*yr2sec, PQ, KQ, MdEff, xi, fbin, Fmin, CC, CCQM, KK, KKQM, A, AQM, B, BQM, Lstr, LstrQM, chi, zlim)
			Ftot_grid_P3[i] = 4.*ma.pi * intg.quad(Fbin_Integrand_GWgas, 0., zmax,  args=Args_gr3)[0]	





		plt.figure()

		cnt1 = plt.plot(Pbasez, np.log10(Ftot_grid_P1), color='black', linewidth=3, linestyle="-")
		cnt2 = plt.plot(Pbasez, np.log10(Ftot_grid_P2), color='blue', linewidth=3, linestyle="--")
		cnt3 = plt.plot(Pbasez, np.log10(Ftot_grid_P3), color='red', linewidth=3, linestyle=":")



		plt.xlabel(r'$\rm{log}_{10}\left[P_{\rm{base}}/\rm{yr}\right]$')
		plt.ylabel(r'$\rm{log}_{10}\left[N_{\rm{tot}}\right]$ (All Sky)')

		plt.legend( [ cnt1[0], cnt2[0], cnt3[0] ], (r'$\rm{log}_{10}\dot{\mathcal{M}}=-1$', r'$\rm{log}_{10}\dot{\mathcal{M}}=-2$', r'$\rm{log}_{10}\dot{\mathcal{M}}=-3$'), loc='upper right', fontsize=18)



		plt.tight_layout()

		#plt.show()
		plt.savefig('Ftot_vs_Pbase_N%g_zmax%g.png'%(Nl,zmax))













if (DoDiags):


	eps1 = 0.01
	eps2 = 0.001
	eps3 = 0.0001
	#eps4 = 1.e-6

	


	zz = np.linspace(0.001, 2.0, 40)
	#zz = np.arange(0.01, zmax, zmax/20.)
	DZ = 1./(4.*ma.pi)

	zlog = np.log10(zz)

	
	fb1 = np.zeros(len(zz))
	fb2 = np.zeros(len(zz))
	fb3 = np.zeros(len(zz))
	

	## TOTAL NUMBER FROM sMALF
	dV = dVdzdOm(zz, h, Om, OL)
	FInt_LF = FInt_smLF(zz, chi, Fmin, h, Om, OL)
	Ntotz   = zz*4.*ma.pi*dV*FInt_LF/(10.**6 * pc2cm)**3
	

	Ftotz1 = np.zeros(len(zz))
	Ftotz2 = np.zeros(len(zz))
	Ftotz3 = np.zeros(len(zz))
	#Ftotz4 = np.zeros(len(zz))
	for i in range(len(zz)):
		Ftotz1[i] = NtotDZ_GWgas(zz[i], Mmx, DZ, Fmin, chi, 1.0*thMn, qmin, eps1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
		Ftotz2[i] = NtotDZ_GWgas(zz[i], Mmx, DZ, Fmin, chi, 1.0*thMn, qmin, eps1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
		Ftotz3[i] = NtotDZ_GWgas(zz[i], Mmx, DZ, Fmin, chi, 1.0*thMn, qmin, eps1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
		fb1[i] = FtotDZ_GWgas(zz[i], Mmx, DZ, Fmin, chi, 1.0*thMn, qmin, eps1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
		fb2[i] = FtotDZ_GWgas(zz[i], Mmx, DZ, Fmin, chi, 1.0*thMn, qmin, eps1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
		fb3[i] = FtotDZ_GWgas(zz[i], Mmx, DZ, Fmin, chi, 1.0*thMn, qmin, eps1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	


	

	# Pmin = np.zeros(len(zz))
	# Pmax = np.zeros(len(zz))
	# for i in range(len(zz)):
	# 	Pmin[i] = max(Pmn[i], PNyq)
	# 	Pmax[i] = min(Pmx[i], Pbase)

	# 	if (Pmn[i] >= Pbase):
	# 		print "Pmax less than Pmin (taken care of)"


	plt.figure(figsize=[12,6])

	# plt.subplot(221)
	# plt.plot(zz,dV/(10**6 * pc2cm)**3)


	# plt.subplot(221)
	# plt.plot(zz,LF)

	plt.subplot(221)
	plt.plot(zz,np.log10(FInt_LF), color='black', linewidth=3)
	plt.ylabel(r'$\rm{log}_{10}[\rm{N} \rm{Mpc}^{-3}]$')
	plt.xlabel(r'z')

	plt.subplot(222)
	plt.plot(zz,np.log10(fb1), color='black', linewidth=3)
	plt.ylabel(r'$\rm{log}_{10}[f_{\rm{bin}}]$')
	plt.xlabel(r'z')



	plt.subplot(223)
	plt.plot(zz, Ntotz, color='black', linewidth=3)
	#plt.ylabel(r'$\rm{log}_{10}[N_{\rm{mm}}]$')
	plt.ylabel(r'$N_{\rm{mm}}$')
	plt.xlabel(r'z')

	plt.subplot(224)
	plt.plot(zz, np.log10(Ftotz1), color='black', linewidth=3)
	plt.ylabel(r'$\rm{log}_{10}[f_{\rm{EHT}}]$')
	plt.xlabel(r'z')

	# plt.subplot(221)
	# plt.plot(zz, Pmn/yr2sec)
	# plt.plot(zz, Pbase/(1.+zz)/yr2sec)

	plt.tight_layout()

	plt.savefig("NEHT.png")
	#plt.show()


	#### vary eps
	plt.figure(figsize=[12,6])



	plt.subplot(211)
	plt.plot(zz,np.log10(fb1), color='black', linewidth=3)
	plt.plot(zz,np.log10(fb2), color='blue', linewidth=3)
	plt.plot(zz,np.log10(fb3), color='red', linewidth=3)
	#plt.plot(zz,np.log10(fb4), color='green', linewidth=3)
	plt.ylabel(r'$\rm{log}_{10}[f_{\rm{bin}}]$')
	plt.xlabel(r'z')



	plt.subplot(212)
	plt.plot(zz, np.log10(Ftotz1), color='black', linewidth=3)
	plt.plot(zz, np.log10(Ftotz2), color='blue', linewidth=3)
	plt.plot(zz, np.log10(Ftotz3), color='red', linewidth=3)
	#plt.plot(zz, np.log10(Ftotz4), color='green', linewidth=3)
	plt.ylabel(r'$\rm{log}_{10}[f_{\rm{EHT}}]$')
	plt.xlabel(r'z')




	# plt.subplot(211)
	# plt.plot(zz,(fb1), color='black', linewidth=3)
	# plt.plot(zz,(fb2), color='blue', linewidth=3)
	# plt.plot(zz,(fb3), color='red', linewidth=3)
	# plt.plot(zz,(fb4), color='green', linewidth=3)
	# plt.ylabel(r'$\rm{log}_{10}[f_{\rm{bin}}]$')
	# plt.xlabel(r'z')



	# plt.subplot(212)
	# plt.plot(zz, (Ftotz1), color='black', linewidth=3)
	# plt.plot(zz, (Ftotz2), color='blue', linewidth=3)
	# plt.plot(zz, (Ftotz3), color='red', linewidth=3)
	# plt.plot(zz, (Ftotz4), color='green', linewidth=3)
	# plt.ylabel(r'$\rm{log}_{10}[f_{\rm{EHT}}]$')
	# plt.xlabel(r'z')



	plt.tight_layout()

	plt.savefig("Vary_eps.png")






	

	# plt.figure()


	# plt.plot(zz, (Ftotz1), color='black', linewidth=3)
	
	# #plt.ylabel(r'$\rm{log}_{10}[f_{\rm{EHT}}]$')
	# #plt.xlabel(r'z')




	# plt.tight_layout()

	# plt.savefig("daigs_1.png")






