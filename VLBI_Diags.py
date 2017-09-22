import numpy as np


import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt


import scipy.integrate as intg
import math as ma
import VLBI_IntFuncs_V2 as IFs 
from VLBI_IntFuncs_V2 import *


##INT options
TrapInt = True

###PLOTTING OPTIONS











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

###Cosmology
h=0.7
Om = 0.3
OL=0.7


Mmx = 40000.*10.**10 ## jsut to not limit lum function - doesnt change anyhting when set at 2*10^10 number
Mmax = 40000.*10.**10*Msun
Mmin= 0.**5*Msun 



###SED scaling
chi  = 0.5 #Elvis 1994 Radio Loud AGN - nulnu propto nu^(0.9) and #0.1  #chi L_(408MHz) = L_mm (specifc fluxes)



###DEFUNCT
DZ = 1.0

##FREE (and KEY) PARAMETERS
##Overall pop params (keep 1)
xi = 1.0
fbin = 1.0



#Accretion Params
KQ = 10.**(-1.0) ## sets number of pc at which RL turns on
eps = 1.0#10**(-3.75)  ## sets migration (accretion rate in CBD pushing two together)
f_Edd = 1.0  ## sets L to M connection (accretion rate onto shining BH)
MdEff = 0.1

## binary pop params
qmin_EHT = 0.1   ### qminof EHT sample
qmin_POP = 0.01  ### qmin of all MBHBS
zeval = 0.5  #eval at this z
zmax = 5.0 ### integrateo out to zmax=5.0

##Insturment params
Fmin = 10.0 * mJy2cgs
thMn = 1.0 * mu_as2rad 
Pbase = 10.0*yr2sec


nuVbnd = c/(5.45*10**(-5))
F0opt  = 3.636*10**(-20)*nuVbnd 
maglim = 24.5
Fmin_opt = 10.**(-maglim/2.5)*F0opt/nuVbnd 







##Lum func * dN/dVdz

Lmms = np.linspace(Lmin(0.01, h, Om, OL, Fmin)-1.0, 30., 200.)
LLs = [24., 25., 26., 27., 28., 30., 32.]
zzs = np.linspace(-2.0, np.log10(8.0), 200.)
ZZs = [0.1, 0.5, 1.0, 2.0, 3.0]

# LumDV1 = NtotDZ_Integrand(Lmms, 0.1, chi, h, Om, OL)
# LumDV2 = NtotDZ_Integrand(Lmms, 1.0, chi, h, Om, OL)
# LumDV3 = NtotDZ_Integrand(Lmms, 2.0, chi, h, Om, OL)
# LumDV4 = NtotDZ_Integrand(Lmms, 3.0, chi, h, Om, OL)
# LumDV5 = NtotDZ_Integrand(Lmms, 4.0, chi, h, Om, OL)

LumDV1 = smLF(Lmms, ZZs[0], chi)
LumDV2 = smLF(Lmms, ZZs[1], chi)
LumDV3 = smLF(Lmms, ZZs[2], chi)
LumDV4 = smLF(Lmms, ZZs[3], chi)
LumDV5 = smLF(Lmms, ZZs[4], chi)


zLumDV1 = smLF(LLs[0], 10.**zzs, chi)
zLumDV2 = smLF(LLs[1], 10.**zzs, chi)
zLumDV3 = smLF(LLs[2], 10.**zzs, chi)
zLumDV4 = smLF(LLs[3], 10.**zzs, chi)
zLumDV5 = smLF(LLs[4], 10.**zzs, chi)





LumVol1 = 4.*np.pi*ZZs[0]*dVdzdOm(ZZs[0], h, Om, OL) * smLF(Lmms, ZZs[0], chi)/(10.**6 * pc2cm)**3
LumVol2 = 4.*np.pi*ZZs[1]*dVdzdOm(ZZs[1], h, Om, OL) * smLF(Lmms, ZZs[1], chi)/(10.**6 * pc2cm)**3
LumVol3 = 4.*np.pi*ZZs[2]*dVdzdOm(ZZs[2], h, Om, OL) * smLF(Lmms, ZZs[2], chi)/(10.**6 * pc2cm)**3
LumVol4 = 4.*np.pi*ZZs[3]*dVdzdOm(ZZs[3], h, Om, OL) * smLF(Lmms, ZZs[3], chi)/(10.**6 * pc2cm)**3
LumVol5 = 4.*np.pi*ZZs[4]*dVdzdOm(ZZs[4], h, Om, OL) * smLF(Lmms, ZZs[4], chi)/(10.**6 * pc2cm)**3

zLumVol1 = 4.*np.pi*10.**zzs*dVdzdOm(10.**zzs, h, Om, OL) *smLF(LLs[0], 10.**zzs, chi)/(10.**6 * pc2cm)**3
zLumVol2 = 4.*np.pi*10.**zzs*dVdzdOm(10.**zzs, h, Om, OL) *smLF(LLs[1], 10.**zzs, chi)/(10.**6 * pc2cm)**3
zLumVol3 = 4.*np.pi*10.**zzs*dVdzdOm(10.**zzs, h, Om, OL) *smLF(LLs[2], 10.**zzs, chi)/(10.**6 * pc2cm)**3
zLumVol4 = 4.*np.pi*10.**zzs*dVdzdOm(10.**zzs, h, Om, OL) *smLF(LLs[3], 10.**zzs, chi)/(10.**6 * pc2cm)**3
zLumVol5 = 4.*np.pi*10.**zzs*dVdzdOm(10.**zzs, h, Om, OL) *smLF(LLs[4], 10.**zzs, chi)/(10.**6 * pc2cm)**3
zLumVol6 = 4.*np.pi*10.**zzs*dVdzdOm(10.**zzs, h, Om, OL) *smLF(LLs[5], 10.**zzs, chi)/(10.**6 * pc2cm)**3
zLumVol7 = 4.*np.pi*10.**zzs*dVdzdOm(10.**zzs, h, Om, OL) *smLF(LLs[6], 10.**zzs, chi)/(10.**6 * pc2cm)**3




# ### F(M,z)
Ms = np.linspace(5.0, np.log10(Mmx), 200.)
Mbs = Lmm2Mbn(Lmms, Mmx, 1.0)

eps_CBD4 = 10.0
eps_CBD3 = 1.0
eps_CBD2 = 0.1
eps_CBD1 = 0.01
FofM4 = np.empty(len(Ms))
FofM3 = np.empty(len(Ms))
FofM2 = np.empty(len(Ms))
FofM1 = np.empty(len(Ms))

FofL4 = np.empty(len(Lmms))
FofL3 = np.empty(len(Lmms))
FofL2 = np.empty(len(Lmms))
FofL1 = np.empty(len(Lmms))

Foz7 = np.empty(len(zzs))
Foz6 = np.empty(len(zzs))
Foz5 = np.empty(len(zzs))
Foz4 = np.empty(len(zzs))
Foz3 = np.empty(len(zzs))
Foz2 = np.empty(len(zzs))
Foz1 = np.empty(len(zzs))


IntGL4 = np.empty(len(Lmms))
IntGL3 = np.empty(len(Lmms))
IntGL2 = np.empty(len(Lmms))
IntGL1 = np.empty(len(Lmms))

ML1 = Lmm2Mbn(LLs[0], Mmx, f_Edd)
ML2 = Lmm2Mbn(LLs[1], Mmx, f_Edd)
ML3 = Lmm2Mbn(LLs[2], Mmx, f_Edd)
ML4 = Lmm2Mbn(LLs[3], Mmx, f_Edd)
ML5 = Lmm2Mbn(LLs[4], Mmx, f_Edd)
ML6 = Lmm2Mbn(LLs[5], Mmx, f_Edd)
ML7 = Lmm2Mbn(LLs[6], Mmx, f_Edd)

for i in range(len(zzs)):
	Foz2[i] = fbin_GWgas(10.**zzs[i], ML2*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD4, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	Foz3[i] = fbin_GWgas(10.**zzs[i], ML3*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD3, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	Foz4[i] = fbin_GWgas(10.**zzs[i], ML4*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD2, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	Foz5[i] = fbin_GWgas(10.**zzs[i], ML5*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD1, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	Foz6[i] = fbin_GWgas(10.**zzs[i], ML6*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD1, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	Foz7[i] = fbin_GWgas(10.**zzs[i], ML7*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD1, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)


for i in range(len(Ms)):
	FofM4[i] = fbin_GWgas(zeval, 10.**Ms[i]*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD4, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	FofM3[i] = fbin_GWgas(zeval, 10.**Ms[i]*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD3, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	FofM2[i] = fbin_GWgas(zeval, 10.**Ms[i]*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD2, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	FofM1[i] = fbin_GWgas(zeval, 10.**Ms[i]*Msun, thMn, qmin_EHT, qmin_POP, eps_CBD1, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	
for i in range(len(Lmms)):
	FofL4[i] = FbinofLmm(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD4, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	FofL3[i] = FbinofLmm(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD3, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	FofL2[i] = FbinofLmm(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD2, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	FofL1[i] = FbinofLmm(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
#			   FbinofLmm(Lmm, z, Mmx, chi, thMn, qmin, eps, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL):

	IntGL4[i] = Fbin_Integrand_GWgas(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD4, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	IntGL3[i] = Fbin_Integrand_GWgas(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD3, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	IntGL2[i] = Fbin_Integrand_GWgas(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD2, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)
	IntGL1[i] = Fbin_Integrand_GWgas(Lmms[i], zeval, Mmx, chi, thMn, qmin_EHT, qmin_POP, eps_CBD1, f_Edd, Pbase, KQ, MdEff, xi, fbin, h, Om, OL)






# #NEHT(z)
Nofz4 = np.empty(len(zzs))
Nofz3 = np.empty(len(zzs))
Nofz2 = np.empty(len(zzs))
Nofz1 = np.empty(len(zzs))

for i in range(len(zzs)):
	Nofz1[i] = OptNEHT_Trap_dL(10.**zzs[i], Mmx, eps_CBD1, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL)
	Nofz2[i] = OptNEHT_Trap_dL(10.**zzs[i], Mmx, eps_CBD2, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL)
	Nofz3[i] = OptNEHT_Trap_dL(10.**zzs[i], Mmx, eps_CBD3, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL)
	Nofz4[i] = OptNEHT_Trap_dL(10.**zzs[i], Mmx, eps_CBD4, f_Edd, KQ, Fmin, chi, thMn, qmin_EHT, qmin_POP, Pbase, xi, fbin, h, Om, OL)


#LumDV4 = NtotDZ_Integrand(Lmms, 3.0, chi, h, Om, OL)

### INtegrated over Lmm
# NtotDL_Trap_RLF(z, Fmin, chi, h, Om, OL):



Mb1 = Lmm2Mbn(Lmms, Mmx, 1.e-5)
Mb2 = Lmm2Mbn(Lmms, Mmx, 1.e-3)
Mb3 = Lmm2Mbn(Lmms, Mmx, 0.1)
Mb4 = Lmm2Mbn(Lmms, Mmx, 1.0)
Mb5 = Lmm2Mbn(Lmms, Mmx, 10.0)

plt.figure()


p1 = plt.plot(Lmms , np.log10(Mb1), color='#1b9e77', linewidth=3)
p2 = plt.plot(Lmms , np.log10(Mb2), color='#d95f02', linewidth=3)
p3 = plt.plot(Lmms , np.log10(Mb3), color='#7570b3', linewidth=3, linestyle='--')
p4 = plt.plot(Lmms , np.log10(Mb4), color='#e7298a', linewidth=3)
p5 = plt.plot(Lmms , np.log10(Mb5), color='gray', linewidth=3)
#plt.axhline(0.0, color='black', linestyle=':')
plt.xlabel(r'$\rm{log}_{10}[L_{bol}]$')
plt.ylabel(r'$\rm{log}_{10}[M_{bin}/M_{\odot}]$')

plt.figlegend([p1[0], p2[0], p3[0], p4[0], p5[0]],(r'$f_{Edd} = 1.e-5$', r'$f_{Edd} = 1.e-3$', r'$f_{Edd} = 0.1$', r'$f_{Edd} = 1.0$',r'$f_{Edd} = 10.0$'),(0.69, 0.18), fontsize=12)

plt.tight_layout()

plt.savefig("Mbin_vs_Lmm.png")









#plt.figure(figsize=[12,4])
plt.figure()

### Plot LumFncX Vol element
plt.subplot(121)
pl1 = plt.plot(Lmms, np.log10(LumVol1), color='#1b9e77', linewidth=3)
pl2 = plt.plot(Lmms, np.log10(LumVol2), color='#d95f02', linewidth=3)
pl3 = plt.plot(Lmms, np.log10(LumVol3), color='#7570b3', linewidth=3)
pl4 = plt.plot(Lmms, np.log10(LumVol4), color='#e7298a', linewidth=3)
pl5 = plt.plot(Lmms, np.log10(LumVol5), color='gray', linewidth=3)
plt.ylabel(r'$\rm{log}_{10}[\frac{d^2N}{dLdV} \rm{Mpc}^{-3}]$')
plt.xlabel(r'$\rm{log}_{10}[L_{mm}]$')

plt.subplot(122)
pl1 = plt.plot(np.log10(Mbs), np.log10(LumVol1), color='#1b9e77', linewidth=3)
pl2 = plt.plot(np.log10(Mbs), np.log10(LumVol2), color='#d95f02', linewidth=3)
pl3 = plt.plot(np.log10(Mbs), np.log10(LumVol3), color='#7570b3', linewidth=3)
pl4 = plt.plot(np.log10(Mbs), np.log10(LumVol4), color='#e7298a', linewidth=3)
pl5 = plt.plot(np.log10(Mbs), np.log10(LumVol5), color='gray', linewidth=3)


plt.axhline(0.0, color='black', linestyle=':')

plt.figlegend([pl1[0], pl2[0], pl3[0], pl4[0], pl5[0]],(r'$z = 0.1$', r'$z = 0.5$', r'$z = 1.0$',r'$z = 2.0$',r'$z = 3.0$'),(0.81, 0.71), fontsize=12)

plt.ylabel(r'$\rm{log}_{10}[\frac{d^2N}{dLdV} \frac{d^2V}{dzd\Omega}]$')
plt.xlabel(r'$\rm{log}_{10}[M/M_{\odot}]$')



plt.tight_layout()


plt.savefig("Lmm_LumFuncXVolElement_vs_L_vs_M.png")















#Lm func times col element
plt.figure()

pl1 = plt.plot(np.log10(Mbs), np.log10(LumDV1), color='#1b9e77', linewidth=3)
pl2 = plt.plot(np.log10(Mbs), np.log10(LumDV2), color='#d95f02', linewidth=3)
pl3 = plt.plot(np.log10(Mbs), np.log10(LumDV3), color='#7570b3', linewidth=3)
pl4 = plt.plot(np.log10(Mbs), np.log10(LumDV4), color='#e7298a', linewidth=3)
pl5 = plt.plot(np.log10(Mbs), np.log10(LumDV5), color='gray', linewidth=3)

plt.axhline(0.0, color='black', linestyle=':')

plt.figlegend([pl1[0], pl2[0], pl3[0], pl4[0], pl5[0]],(r'$z = 0.1$', r'$z = 0.5$', r'$z = 1.0$',r'$z = 2.0$',r'$z = 3.0$'),(0.81, 0.71), fontsize=12)

plt.ylabel(r'$\rm{log}_{10}[\frac{d^2N}{dLdV} \rm{Mpc}^{-3}]$')
# plt.xlabel(r'$\rm{log}_{10}[L_{mm}]$')
plt.xlabel(r'$\rm{log}_{10}[M/M_{\odot}]$')



plt.tight_layout()


plt.savefig("Lmm_LumFunc_vs_L.png")














#plt.figure(figsize=[12,4])
plt.figure()

### Plot LumFnc time Vole element
#plt.subplot(121)
pl1 = plt.plot(np.log10(Lmm2Mbn(Lmms, Mmx, 1.0)), np.log10(LumDV1), color='#1b9e77', linewidth=3)
pl2 = plt.plot(np.log10(Lmm2Mbn(Lmms, Mmx, 1.0)), np.log10(LumDV2), color='#d95f02', linewidth=3)
pl3 = plt.plot(np.log10(Lmm2Mbn(Lmms, Mmx, 1.0)), np.log10(LumDV3), color='#7570b3', linewidth=3)
pl4 = plt.plot(np.log10(Lmm2Mbn(Lmms, Mmx, 1.0)), np.log10(LumDV4), color='#e7298a', linewidth=3)
#pl5 = plt.plot(np.log10(Lmm2Mbn(Lmms, Mmx, 0.1)), np.log10(LumDV5), color='gray', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')

plt.figlegend([pl1[0], pl2[0], pl3[0], pl4[0], pl5[0]],(r'$z = 0.1$', r'$z = 0.5$', r'$z = 1.0$',r'$z = 2.0$',r'$z = 3.0$'),(0.81, 0.71), fontsize=12)

plt.ylabel(r'$\rm{log}_{10}[\frac{d^2N}{dLdV} \rm{Mpc}^{-3}]$')
plt.xlabel(r'$\rm{log}_{10}[M]$')



plt.tight_layout()


plt.savefig("LumFunc_vs_Mbin.png")









plt.figure()

#plt.subplot(122)
ppl1 = plt.plot(zzs, np.log10(zLumDV1), color='#1b9e77', linewidth=3)
ppl2 = plt.plot(zzs, np.log10(zLumDV2), color='#d95f02', linewidth=3)
ppl3 = plt.plot(zzs, np.log10(zLumDV3), color='#7570b3', linewidth=3)
ppl4 = plt.plot(zzs, np.log10(zLumDV4), color='#e7298a', linewidth=3)
ppl5 = plt.plot(zzs, np.log10(zLumDV5), color='gray', linewidth=3)

plt.axhline(0.0, color='black', linestyle=':')

plt.figlegend([ppl1[0], ppl2[0], ppl3[0], ppl4[0], ppl5[0]],(r'$\log_{10}{L_{mm}} = 24$', r'$\log_{10}{L_{mm}} = 25$', r'$\log_{10}{L_{mm}} = 26$', r'$\log_{10}{L_{mm}} = 26$', r'$\log_{10}{L_{mm}} = 28$'),(0.74, 0.695), fontsize=12)

plt.ylim(-12, -3.)
plt.ylabel(r'$\rm{log}_{10}[\frac{d^2N}{dLdV} \rm{Mpc}^{-3}]$')
plt.xlabel(r'$\rm{log}_{10}[z]$')



plt.tight_layout()


plt.savefig("Lmm_LumFunc_vs_z.png")












plt.figure(figsize=[10,8])

plt.subplot(221)
#ppl1 = plt.plot(10.**zzs, np.log10(zFoz1), color='#1b9e77', linewidth=3)
ppl2 = plt.plot(10.**zzs, np.log10(Foz2), color='#d95f02', linewidth=3)
ppl3 = plt.plot(10.**zzs, np.log10(Foz3), color='#7570b3', linewidth=3)
ppl4 = plt.plot(10.**zzs, np.log10(Foz4), color='#e7298a', linewidth=3)
ppl5 = plt.plot(10.**zzs, np.log10(Foz5), color='#66a61e', linewidth=3)
ppl6 = plt.plot(10.**zzs, np.log10(Foz6), color='#e6ab02', linewidth=3)
ppl7 = plt.plot(10.**zzs, np.log10(Foz7), color='#a6761d', linewidth=3)
plt.ylabel(r'$\rm{log}_{10}[\mathcal{F}]$')
plt.xlabel(r'$z$')
plt.ylim(-5, 0.)
# plt.xlim(0.0,3.0)




plt.subplot(222)
#ppl1 = plt.plot(10.**zzs, np.log10(zLumVol1), color='#1b9e77', linewidth=3)
ppl2 = plt.plot(10.**zzs, np.log10(zLumVol2), color='#d95f02', linewidth=3)
ppl3 = plt.plot(10.**zzs, np.log10(zLumVol3), color='#7570b3', linewidth=3)
ppl4 = plt.plot(10.**zzs, np.log10(zLumVol4), color='#e7298a', linewidth=3)
ppl5 = plt.plot(10.**zzs, np.log10(zLumVol5), color='#66a61e', linewidth=3)
ppl6 = plt.plot(10.**zzs, np.log10(zLumVol6), color='#e6ab02', linewidth=3)
ppl7 = plt.plot(10.**zzs, np.log10(zLumVol7), color='#a6761d', linewidth=3)
plt.ylabel(r'$\rm{log}_{10}[4 \pi z \frac{d^2N}{dLdV}\frac{d^2V}{dzd\Omega}]$')
plt.xlabel(r'$z$')
plt.ylim(-2, 8.)
# plt.xlim(0.0,3.0)
plt.axhline(0.0, color='black', linestyle=':')



plt.subplot(223)
#ppl1 = plt.plot(10.**zzs, np.log10(zFoz1), color='#1b9e77', linewidth=3)
ppl2 = plt.plot(10.**zzs, np.log10(Foz2*zLumVol2), color='#d95f02', linewidth=3)
ppl3 = plt.plot(10.**zzs, np.log10(Foz3*zLumVol3), color='#7570b3', linewidth=3)
ppl4 = plt.plot(10.**zzs, np.log10(Foz4*zLumVol4), color='#e7298a', linewidth=3)
ppl5 = plt.plot(10.**zzs, np.log10(Foz5*zLumVol5), color='#66a61e', linewidth=3)
ppl6 = plt.plot(10.**zzs, np.log10(Foz6*zLumVol6), color='#e6ab02', linewidth=3)
ppl7 = plt.plot(10.**zzs, np.log10(Foz7*zLumVol7), color='#a6761d', linewidth=3)

plt.ylim(-5, 5.)
# plt.xlim(0.0,3.0)
plt.ylabel(r'$\rm{log}_{10}[4 \pi z \frac{d^2N}{dLdV}\frac{d^2V}{dzd\Omega} \times \mathcal{F}]$')
plt.xlabel(r'$z$')


plt.figlegend([ ppl2[0], ppl3[0], ppl4[0], ppl5[0], ppl6[0], ppl7[0]],(r'$\log_{10}{L_{mm}} = %g$' %LLs[1], r'$\log_{10}{L_{mm}} =  %g$' %LLs[2], r'$\log_{10}{L_{mm}} =  %g$' %LLs[3], r'$\log_{10}{L_{mm}} =  %g$' %LLs[4], r'$\log_{10}{L_{mm}} =  %g$' %LLs[5], r'$\log_{10}{L_{mm}} =  %g$' %LLs[6] ),loc="lower right", fontsize=12)

plt.tight_layout()


plt.savefig("Lmm_LumFuncXVolElement_vs_z.png")










#4 panel


plt.figure(figsize=[10,8])

# plt.subplot(221)
# pl1 = plt.plot(Ms, np.log10(FofM1), color='#1b9e77', linewidth=3)
# pl2 = plt.plot(Ms, np.log10(FofM2), color='#d95f02', linewidth=3)
# pl3 = plt.plot(Ms, np.log10(FofM3), color='#7570b3', linewidth=3)
# pl4 = plt.plot(Ms, np.log10(FofM4), color='#e7298a', linewidth=3)
# plt.axhline(0.0, color='black', linestyle=':')
# #plt.ylim(-12, -3.)
# plt.figlegend([pl1[0], pl2[0], pl3[0], pl4[0]],(r'$\dot{\mathcal{M}} = 0.01$', r'$\dot{\mathcal{M}} = 0.1$', r'$\dot{\mathcal{M}} = 1.0$',r'$\dot{\mathcal{M}} = 10.0$'),'upper right')
# plt.ylabel(r'$\rm{log}_{10}[\mathcal{F}(M, %g)]$' %zeval)
# plt.xlabel(r'$\rm{log}_{10}[M_{bin}/M_{\odot}]$')

plt.subplot(221)
pl1 = plt.plot(Ms, FofM1, color='#1b9e77', linewidth=3)
pl2 = plt.plot(Ms, FofM2, color='#d95f02', linewidth=3)
pl3 = plt.plot(Ms, FofM3, color='#7570b3', linewidth=3)
pl4 = plt.plot(Ms, FofM4, color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
#plt.ylim(-12, -3.)
plt.figlegend([pl1[0], pl2[0], pl3[0], pl4[0]],(r'$\dot{\mathcal{M}} = 0.01$', r'$\dot{\mathcal{M}} = 0.1$', r'$\dot{\mathcal{M}} = 1.0$',r'$\dot{\mathcal{M}} = 10.0$'),(0.12,0.82), fontsize=11)
plt.ylabel(r'$\rm{log}_{10}[\mathcal{F}(M, %g)]$' %zeval)
plt.xlabel(r'$\rm{log}_{10}[M_{bin}/M_{\odot}]$')



plt.subplot(222)
plt.plot(Lmms, np.log10(FofL1), color='#1b9e77', linewidth=3)
plt.plot(Lmms, np.log10(FofL2), color='#d95f02', linewidth=3)
plt.plot(Lmms, np.log10(FofL3), color='#7570b3', linewidth=3)
plt.plot(Lmms, np.log10(FofL4), color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
plt.ylim(-6.0, 0.0)
plt.ylabel(r'$\rm{log}_{10}[\mathcal{F}(L, %g)]$' %zeval)
plt.xlabel(r'$\rm{log}_{10}[L_{mm}]$')


# plt.subplot(222)
# plt.plot(Lmms, FofL1, color='#1b9e77', linewidth=3)
# plt.plot(Lmms, FofL2, color='#d95f02', linewidth=3)
# plt.plot(Lmms, FofL3, color='#7570b3', linewidth=3)
# plt.plot(Lmms, FofL4, color='#e7298a', linewidth=3)
# plt.axhline(0.0, color='black', linestyle=':')
# #plt.ylim(-12, -3.)
# plt.ylabel(r'$\rm{log}_{10}[\mathcal{F}(L, %g)]$' %zeval)
# plt.xlabel(r'$\rm{log}_{10}[L_{mm}]$')





plt.subplot(223)
plt.plot(Lmms, np.log10(IntGL1), color='#1b9e77', linewidth=3)
plt.plot(Lmms, np.log10(IntGL2), color='#d95f02', linewidth=3)
plt.plot(Lmms, np.log10(IntGL3), color='#7570b3', linewidth=3)
plt.plot(Lmms, np.log10(IntGL4), color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
#plt.ylim(-12, -3.)
plt.ylabel(r'$\rm{log}_{10}[ \frac{d^2V}{dz d\Omega} \frac{d^2N}{dLdV} \mathcal{F}(L, %g)]$' %zeval)
plt.xlabel(r'$\rm{log}_{10}[L_{mm}]$')



plt.subplot(224)
plt.plot(zzs, np.log10(Nofz1), color='#1b9e77', linewidth=3)
plt.plot(zzs, np.log10(Nofz2), color='#d95f02', linewidth=3)
plt.plot(zzs, np.log10(Nofz3), color='#7570b3', linewidth=3)
plt.plot(zzs, np.log10(Nofz4), color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
plt.ylim(-3,6)
plt.ylabel(r'$\rm{log}_{10}[\rm{N_{\rm{EHT}}}(z)]$')
plt.xlabel(r'$\rm{log}_{10}[z]$')




plt.tight_layout()

thMnsv = thMn/mu_as2rad
Pbasesv = Pbase/yr2sec

plt.savefig("NEHT_FprobDiags_vsLmm_fEdd%g_Mdeff%g_amax%gpc_thMn%gmuas_qminEHT%g_qminPOP%g_Pbase%gyr.png" %(f_Edd, MdEff, KQ, thMnsv, qmin_EHT, qmin_POP, Pbasesv))









###SAME AS ABOVE BUT AGAINST MASS(Lmm)

plt.figure(figsize=[10,8])

plt.subplot(221)
pl1 = plt.plot(Ms, FofM1, color='#1b9e77', linewidth=3)
pl2 = plt.plot(Ms, FofM2, color='#d95f02', linewidth=3)
pl3 = plt.plot(Ms, FofM3, color='#7570b3', linewidth=3)
pl4 = plt.plot(Ms, FofM4, color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
#plt.ylim(-12, -3.)
plt.figlegend([pl1[0], pl2[0], pl3[0], pl4[0]],(r'$\dot{\mathcal{M}} = 0.01$', r'$\dot{\mathcal{M}} = 0.1$', r'$\dot{\mathcal{M}} = 1.0$',r'$\dot{\mathcal{M}} = 10.0$'),(0.12,0.82), fontsize=11)
plt.ylabel(r'$\rm{log}_{10}[\mathcal{F}(M, %g)]$' %zeval)
plt.xlabel(r'$\rm{log}_{10}[M_{bin}/M_{\odot}]$')



plt.subplot(222)
plt.plot(np.log10(Mbs), np.log10(FofL1), color='#1b9e77', linewidth=3)
plt.plot(np.log10(Mbs), np.log10(FofL2), color='#d95f02', linewidth=3)
plt.plot(np.log10(Mbs), np.log10(FofL3), color='#7570b3', linewidth=3)
plt.plot(np.log10(Mbs), np.log10(FofL4), color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
plt.ylim(-6.0, 0.0)
plt.ylabel(r'$\rm{log}_{10}[\mathcal{F}(L, %g)]$' %zeval)
plt.xlabel(r'$\rm{log}_{10}[M_{bin}/M_{\odot}]$')






plt.subplot(223)
plt.plot(np.log10(Mbs), np.log10(IntGL1), color='#1b9e77', linewidth=3)
plt.plot(np.log10(Mbs), np.log10(IntGL2), color='#d95f02', linewidth=3)
plt.plot(np.log10(Mbs), np.log10(IntGL3), color='#7570b3', linewidth=3)
plt.plot(np.log10(Mbs), np.log10(IntGL4), color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
#plt.ylim(-12, -3.)
plt.ylabel(r'$\rm{log}_{10}[ \frac{d^2V}{dz d\Omega} \frac{d^2N}{dLdV} \mathcal{F}(L, %g)]$' %zeval)
plt.xlabel(r'$\rm{log}_{10}[M_{bin}/M_{\odot}]$')



plt.subplot(224)
plt.plot(zzs, np.log10(Nofz1), color='#1b9e77', linewidth=3)
plt.plot(zzs, np.log10(Nofz2), color='#d95f02', linewidth=3)
plt.plot(zzs, np.log10(Nofz3), color='#7570b3', linewidth=3)
plt.plot(zzs, np.log10(Nofz4), color='#e7298a', linewidth=3)
plt.axhline(0.0, color='black', linestyle=':')
plt.ylim(-3,6)
plt.ylabel(r'$\rm{log}_{10}[\rm{N_{\rm{EHT}}}(z)]$')
plt.xlabel(r'$\rm{log}_{10}[z]$')




plt.tight_layout()

thMnsv = thMn/mu_as2rad
Pbasesv = Pbase/yr2sec

plt.savefig("NEHT_FprobDiags_vsMass_fEdd%g_Mdeff%g_amax%gpc_thMn%gmuas_qminEHT%g_qminPOP%g_Pbase%gyr.png" %(f_Edd, MdEff, KQ, thMnsv, qmin_EHT, qmin_POP, Pbasesv))


