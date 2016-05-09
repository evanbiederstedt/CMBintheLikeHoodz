
# coding: utf-8

# In[1]:
#
#
# hundred_samples = np.linspace(0.05, 0.5, num=100)
# 
# Planck found \Omega_CDM
# GAVO simulated map set at \Omega_CDM = 0.122
# CAMB default below at omch2=0.122
#


# In[2]:

#
# First output 200 CAMB scalar outputs
#
# 0.005 to 0.05
#


# In[3]:


from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower


# In[4]:
"""
#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(2000, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
for name in powers: 
    print(name)


# In[5]:

# plot the total lensed CMB power spectra versus unlensed, and fractional difference
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)
# Python CL arrays are all zero based (starting at L=0), Note L=0,1 entries will be zero by default.
# The differenent CL are always in the order TT, EE, BB, TE (with BB=0 for unlensed scalar results).
ls = np.arange(totCL.shape[0])
print(ls)
#print(totCL[:30]) # print first 30 totCL

fig, ax = plt.subplots(2,2, figsize = (12,12))
ax[0,0].plot(ls,totCL[:,0], color='k')
ax[0,0].plot(ls,unlensedCL[:,0], color='r')
ax[0,0].set_title('TT')
ax[0,1].plot(ls[2:], 1-unlensedCL[2:,0]/totCL[2:,0]);
ax[0,1].set_title(r'$\Delta TT$')
ax[1,0].plot(ls,totCL[:,1], color='k')
ax[1,0].plot(ls,unlensedCL[:,1], color='r')
ax[1,0].set_title(r'$EE$')
ax[1,1].plot(ls,totCL[:,3], color='k')
ax[1,1].plot(ls,unlensedCL[:,3], color='r')
ax[1,1].set_title(r'$TE$');
for ax in ax.reshape(-1): ax.set_xlim([2,2500])
"""

# In[6]:

twohundred_samples = np.linspace(0.005, 0.05, num=200)
#print(twohundred_samples)


#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers =results.get_cmb_power_spectra(pars)
for name in powers:
    print(name)


"""
    array([ 0.005     ,  0.00522613,  0.00545226,  0.00567839,  0.00590452,
    0.00613065,  0.00635678,  0.00658291,  0.00680905,  0.00703518,
    0.00726131,  0.00748744,  0.00771357,  0.0079397 ,  0.00816583,
    0.00839196,  0.00861809,  0.00884422,  0.00907035,  0.00929648,
    0.00952261,  0.00974874,  0.00997487,  0.01020101,  0.01042714,
    0.01065327,  0.0108794 ,  0.01110553,  0.01133166,  0.01155779,
    0.01178392,  0.01201005,  0.01223618,  0.01246231,  0.01268844,
    0.01291457,  0.0131407 ,  0.01336683,  0.01359296,  0.0138191 ,
    0.01404523,  0.01427136,  0.01449749,  0.01472362,  0.01494975,
    0.01517588,  0.01540201,  0.01562814,  0.01585427,  0.0160804 ,
    0.01630653,  0.01653266,  0.01675879,  0.01698492,  0.01721106,
    0.01743719,  0.01766332,  0.01788945,  0.01811558,  0.01834171,
    0.01856784,  0.01879397,  0.0190201 ,  0.01924623,  0.01947236,
    0.01969849,  0.01992462,  0.02015075,  0.02037688,  0.02060302,
    0.02082915,  0.02105528,  0.02128141,  0.02150754,  0.02173367,
    0.0219598 ,  0.02218593,  0.02241206,  0.02263819,  0.02286432,
    0.02309045,  0.02331658,  0.02354271,  0.02376884,  0.02399497,
    0.02422111,  0.02444724,  0.02467337,  0.0248995 ,  0.02512563,
    0.02535176,  0.02557789,  0.02580402,  0.02603015,  0.02625628,
    0.02648241,  0.02670854,  0.02693467,  0.0271608 ,  0.02738693,
    0.02761307,  0.0278392 ,  0.02806533,  0.02829146,  0.02851759,
    0.02874372,  0.02896985,  0.02919598,  0.02942211,  0.02964824,
    0.02987437,  0.0301005 ,  0.03032663,  0.03055276,  0.03077889,
    0.03100503,  0.03123116,  0.03145729,  0.03168342,  0.03190955,
    0.03213568,  0.03236181,  0.03258794,  0.03281407,  0.0330402 ,
    0.03326633,  0.03349246,  0.03371859,  0.03394472,  0.03417085,
    0.03439698,  0.03462312,  0.03484925,  0.03507538,  0.03530151,
    0.03552764,  0.03575377,  0.0359799 ,  0.03620603,  0.03643216,
    0.03665829,  0.03688442,  0.03711055,  0.03733668,  0.03756281,
    0.03778894,  0.03801508,  0.03824121,  0.03846734,  0.03869347,
    0.0389196 ,  0.03914573,  0.03937186,  0.03959799,  0.03982412,
    0.04005025,  0.04027638,  0.04050251,  0.04072864,  0.04095477,
    0.0411809 ,  0.04140704,  0.04163317,  0.0418593 ,  0.04208543,
    0.04231156,  0.04253769,  0.04276382,  0.04298995,  0.04321608,
    0.04344221,  0.04366834,  0.04389447,  0.0441206 ,  0.04434673,
    0.04457286,  0.04479899,  0.04502513,  0.04525126,  0.04547739,
    0.04570352,  0.04592965,  0.04615578,  0.04638191,  0.04660804,
    0.04683417,  0.0470603 ,  0.04728643,  0.04751256,  0.04773869,
    0.04796482,  0.04819095,  0.04841709,  0.04864322,  0.04886935,
    0.04909548,  0.04932161,  0.04954774,  0.04977387,  0.05      ])
    
"""




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls0 = unlencl[:,0][2:1101]
print(len(cls0))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls1 = unlencl[:,0][2:1101]
print(len(cls1))











pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls2 = unlencl[:,0][2:1101]
print(len(cls2))










pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls3 = unlencl[:,0][2:1101]
print(len(cls3))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls4 = unlencl[:,0][2:1101]
print(len(cls4))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls5 = unlencl[:,0][2:1101]
print(len(cls5))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls6 = unlencl[:,0][2:1101]
print(len(cls6))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls7 = unlencl[:,0][2:1101]
print(len(cls7))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls8 = unlencl[:,0][2:1101]
print(len(cls8))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls9 = unlencl[:,0][2:1101]
print(len(cls9))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls10 = unlencl[:,0][2:1101]
print(len(cls10))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls11 = unlencl[:,0][2:1101]
print(len(cls11))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls12 = unlencl[:,0][2:1101]
print(len(cls12))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls13 = unlencl[:,0][2:1101]
print(len(cls13))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls14 = unlencl[:,0][2:1101]
print(len(cls14))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls15 = unlencl[:,0][2:1101]
print(len(cls15))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls16 = unlencl[:,0][2:1101]
print(len(cls16))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls17 = unlencl[:,0][2:1101]
print(len(cls17))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls18 = unlencl[:,0][2:1101]
print(len(cls18))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls19 = unlencl[:,0][2:1101]
print(len(cls19))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls20 = unlencl[:,0][2:1101]
print(len(cls20))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls21 = unlencl[:,0][2:1101]
print(len(cls21))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls22 = unlencl[:,0][2:1101]
print(len(cls22))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls23 = unlencl[:,0][2:1101]
print(len(cls23))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls24 = unlencl[:,0][2:1101]
print(len(cls24))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls25 = unlencl[:,0][2:1101]
print(len(cls25))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls26 = unlencl[:,0][2:1101]
print(len(cls26))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls27 = unlencl[:,0][2:1101]
print(len(cls27))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls28 = unlencl[:,0][2:1101]
print(len(cls28))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls29 = unlencl[:,0][2:1101]
print(len(cls29))











pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls30 = unlencl[:,0][2:1101]
print(len(cls30))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls31 = unlencl[:,0][2:1101]
print(len(cls31))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls32 = unlencl[:,0][2:1101]
print(len(cls32))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls33 = unlencl[:,0][2:1101]
print(len(cls33))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls34 = unlencl[:,0][2:1101]
print(len(cls34))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls35 = unlencl[:,0][2:1101]
print(len(cls35))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls36 = unlencl[:,0][2:1101]
print(len(cls36))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls37 = unlencl[:,0][2:1101]
print(len(cls37))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls38 = unlencl[:,0][2:1101]
print(len(cls38))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls39 = unlencl[:,0][2:1101]
print(len(cls39))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls40 = unlencl[:,0][2:1101]
print(len(cls40))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls41 = unlencl[:,0][2:1101]
print(len(cls41))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls42 = unlencl[:,0][2:1101]
print(len(cls42))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls43 = unlencl[:,0][2:1101]
print(len(cls43))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls44 = unlencl[:,0][2:1101]
print(len(cls44))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls45 = unlencl[:,0][2:1101]
print(len(cls45))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls46 = unlencl[:,0][2:1101]
print(len(cls46))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls47 = unlencl[:,0][2:1101]
print(len(cls47))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls48 = unlencl[:,0][2:1101]
print(len(cls48))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls49 = unlencl[:,0][2:1101]
print(len(cls49))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls50 = unlencl[:,0][2:1101]
print(len(cls50))


pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls51 = unlencl[:,0][2:1101]
print(len(cls51))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls52 = unlencl[:,0][2:1101]
print(len(cls52))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls53 = unlencl[:,0][2:1101]
print(len(cls53))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls54 = unlencl[:,0][2:1101]
print(len(cls54))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls55 = unlencl[:,0][2:1101]
print(len(cls55))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls56 = unlencl[:,0][2:1101]
print(len(cls56))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls57 = unlencl[:,0][2:1101]
print(len(cls57))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls58 = unlencl[:,0][2:1101]
print(len(cls58))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls59 = unlencl[:,0][2:1101]
print(len(cls59))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls60 = unlencl[:,0][2:1101]
print(len(cls60))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls61 = unlencl[:,0][2:1101]
print(len(cls61))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls62 = unlencl[:,0][2:1101]
print(len(cls62))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls63 = unlencl[:,0][2:1101]
print(len(cls63))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls64 = unlencl[:,0][2:1101]
print(len(cls64))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls65 = unlencl[:,0][2:1101]
print(len(cls65))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls66 = unlencl[:,0][2:1101]
print(len(cls66))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls67 = unlencl[:,0][2:1101]
print(len(cls67))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls68 = unlencl[:,0][2:1101]
print(len(cls68))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls69 = unlencl[:,0][2:1101]
print(len(cls69))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls70 = unlencl[:,0][2:1101]
print(len(cls70))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls71 = unlencl[:,0][2:1101]
print(len(cls71))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls72 = unlencl[:,0][2:1101]
print(len(cls72))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls73 = unlencl[:,0][2:1101]
print(len(cls73))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls74 = unlencl[:,0][2:1101]
print(len(cls74))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls75 = unlencl[:,0][2:1101]
print(len(cls75))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls76 = unlencl[:,0][2:1101]
print(len(cls76))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls77 = unlencl[:,0][2:1101]
print(len(cls77))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls78 = unlencl[:,0][2:1101]
print(len(cls78))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls79 = unlencl[:,0][2:1101]
print(len(cls79))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls80 = unlencl[:,0][2:1101]
print(len(cls80))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls81 = unlencl[:,0][2:1101]
print(len(cls81))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls82 = unlencl[:,0][2:1101]
print(len(cls82))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls83 = unlencl[:,0][2:1101]
print(len(cls83))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls84 = unlencl[:,0][2:1101]
print(len(cls84))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls85 = unlencl[:,0][2:1101]
print(len(cls85))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls86 = unlencl[:,0][2:1101]
print(len(cls86))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls87 = unlencl[:,0][2:1101]
print(len(cls87))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls88 = unlencl[:,0][2:1101]
print(len(cls88))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls89 = unlencl[:,0][2:1101]
print(len(cls89))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls90 = unlencl[:,0][2:1101]
print(len(cls90))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls91 = unlencl[:,0][2:1101]
print(len(cls91))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls92 = unlencl[:,0][2:1101]
print(len(cls92))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls93 = unlencl[:,0][2:1101]
print(len(cls93))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls94 = unlencl[:,0][2:1101]
print(len(cls94))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls95 = unlencl[:,0][2:1101]
print(len(cls95))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls96 = unlencl[:,0][2:1101]
print(len(cls96))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls97 = unlencl[:,0][2:1101]
print(len(cls97))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls98 = unlencl[:,0][2:1101]
print(len(cls98))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls99 = unlencl[:,0][2:1101]
print(len(cls99))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls100 = unlencl[:,0][2:1101]
print(len(cls100))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls101 = unlencl[:,0][2:1101]
print(len(cls101))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls102 = unlencl[:,0][2:1101]
print(len(cls102))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls103 = unlencl[:,0][2:1101]
print(len(cls103))



pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls104 = unlencl[:,0][2:1101]
print(len(cls104))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls105 = unlencl[:,0][2:1101]
print(len(cls105))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls106 = unlencl[:,0][2:1101]
print(len(cls106))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls107 = unlencl[:,0][2:1101]
print(len(cls107))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls108 = unlencl[:,0][2:1101]
print(len(cls108))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls109 = unlencl[:,0][2:1101]
print(len(cls109))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls110 = unlencl[:,0][2:1101]
print(len(cls110))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls111 = unlencl[:,0][2:1101]
print(len(cls111))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls112 = unlencl[:,0][2:1101]
print(len(cls112))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls113 = unlencl[:,0][2:1101]
print(len(cls113))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls114 = unlencl[:,0][2:1101]
print(len(cls114))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls115 = unlencl[:,0][2:1101]
print(len(cls115))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls116 = unlencl[:,0][2:1101]
print(len(cls116))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls117 = unlencl[:,0][2:1101]
print(len(cls117))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls118 = unlencl[:,0][2:1101]
print(len(cls118))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls119 = unlencl[:,0][2:1101]
print(len(cls119))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls120 = unlencl[:,0][2:1101]
print(len(cls120))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls121 = unlencl[:,0][2:1101]
print(len(cls121))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls122 = unlencl[:,0][2:1101]
print(len(cls122))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls123 = unlencl[:,0][2:1101]
print(len(cls123))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls124 = unlencl[:,0][2:1101]
print(len(cls124))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls125 = unlencl[:,0][2:1101]
print(len(cls125))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls126 = unlencl[:,0][2:1101]
print(len(cls126))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls127 = unlencl[:,0][2:1101]
print(len(cls127))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls128 = unlencl[:,0][2:1101]
print(len(cls128))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls129 = unlencl[:,0][2:1101]
print(len(cls129))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls130 = unlencl[:,0][2:1101]
print(len(cls130))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls131 = unlencl[:,0][2:1101]
print(len(cls131))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls132 = unlencl[:,0][2:1101]
print(len(cls132))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls133 = unlencl[:,0][2:1101]
print(len(cls133))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls134 = unlencl[:,0][2:1101]
print(len(cls134))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls135 = unlencl[:,0][2:1101]
print(len(cls135))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls136 = unlencl[:,0][2:1101]
print(len(cls136))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls137 = unlencl[:,0][2:1101]
print(len(cls137))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls138 = unlencl[:,0][2:1101]
print(len(cls138))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls139 = unlencl[:,0][2:1101]
print(len(cls139))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls140 = unlencl[:,0][2:1101]
print(len(cls140))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls141 = unlencl[:,0][2:1101]
print(len(cls141))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls142 = unlencl[:,0][2:1101]
print(len(cls142))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls143 = unlencl[:,0][2:1101]
print(len(cls143))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls144 = unlencl[:,0][2:1101]
print(len(cls144))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls145 = unlencl[:,0][2:1101]
print(len(cls145))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls146 = unlencl[:,0][2:1101]
print(len(cls146))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls147 = unlencl[:,0][2:1101]
print(len(cls147))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls148 = unlencl[:,0][2:1101]
print(len(cls148))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls149 = unlencl[:,0][2:1101]
print(len(cls149))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls150 = unlencl[:,0][2:1101]
print(len(cls150))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls151 = unlencl[:,0][2:1101]
print(len(cls151))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls152 = unlencl[:,0][2:1101]
print(len(cls152))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls153 = unlencl[:,0][2:1101]
print(len(cls153))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls154 = unlencl[:,0][2:1101]
print(len(cls154))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls155 = unlencl[:,0][2:1101]
print(len(cls155))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls156 = unlencl[:,0][2:1101]
print(len(cls156))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls157 = unlencl[:,0][2:1101]
print(len(cls157))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls158 = unlencl[:,0][2:1101]
print(len(cls158))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls159 = unlencl[:,0][2:1101]
print(len(cls159))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls160 = unlencl[:,0][2:1101]
print(len(cls160))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls161 = unlencl[:,0][2:1101]
print(len(cls161))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls162 = unlencl[:,0][2:1101]
print(len(cls162))








pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls163 = unlencl[:,0][2:1101]
print(len(cls163))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls164 = unlencl[:,0][2:1101]
print(len(cls164))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls165 = unlencl[:,0][2:1101]
print(len(cls165))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls166 = unlencl[:,0][2:1101]
print(len(cls166))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls167 = unlencl[:,0][2:1101]
print(len(cls167))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls168 = unlencl[:,0][2:1101]
print(len(cls168))







pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls169 = unlencl[:,0][2:1101]
print(len(cls169))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls170 = unlencl[:,0][2:1101]
print(len(cls170))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls171 = unlencl[:,0][2:1101]
print(len(cls171))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls172 = unlencl[:,0][2:1101]
print(len(cls172))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls173 = unlencl[:,0][2:1101]
print(len(cls173))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls174 = unlencl[:,0][2:1101]
print(len(cls174))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls175 = unlencl[:,0][2:1101]
print(len(cls175))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls176 = unlencl[:,0][2:1101]
print(len(cls176))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls177 = unlencl[:,0][2:1101]
print(len(cls177))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls178 = unlencl[:,0][2:1101]
print(len(cls178))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls179 = unlencl[:,0][2:1101]
print(len(cls179))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls180 = unlencl[:,0][2:1101]
print(len(cls180))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls181 = unlencl[:,0][2:1101]
print(len(cls181))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls182 = unlencl[:,0][2:1101]
print(len(cls182))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls183 = unlencl[:,0][2:1101]
print(len(cls183))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls184 = unlencl[:,0][2:1101]
print(len(cls184))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls185 = unlencl[:,0][2:1101]
print(len(cls185))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls186 = unlencl[:,0][2:1101]
print(len(cls186))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls187 = unlencl[:,0][2:1101]
print(len(cls187))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls188 = unlencl[:,0][2:1101]
print(len(cls188))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls189 = unlencl[:,0][2:1101]
print(len(cls189))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls190 = unlencl[:,0][2:1101]
print(len(cls190))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls191 = unlencl[:,0][2:1101]
print(len(cls191))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls192 = unlencl[:,0][2:1101]
print(len(cls192))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls193 = unlencl[:,0][2:1101]
print(len(cls193))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls194 = unlencl[:,0][2:1101]
print(len(cls194))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls195 = unlencl[:,0][2:1101]
print(len(cls195))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls196 = unlencl[:,0][2:1101]
print(len(cls196))




pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls197 = unlencl[:,0][2:1101]
print(len(cls197))






pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls198 = unlencl[:,0][2:1101]
print(len(cls198))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls199 = unlencl[:,0][2:1101]
print(len(cls199))





pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(514, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']
ls = np.arange(unlencl.shape[0])
print(ls)
print(len(ls))
#
# plot of spectrum

cls200 = unlencl[:,0][2:1101]
print(len(cls200))



"""
0.005
0.00522613065327
0.00545226130653
0.0056783919598
0.00590452261307
0.00613065326633
0.0063567839196
0.00658291457286
0.00680904522613
0.0070351758794
0.00726130653266
0.00748743718593
0.0077135678392
0.00793969849246
0.00816582914573
0.00839195979899
0.00861809045226
0.00884422110553
0.00907035175879
0.00929648241206
0.00952261306533
0.00974874371859
0.00997487437186
0.0102010050251
0.0104271356784
0.0106532663317
0.0108793969849
0.0111055276382
0.0113316582915
0.0115577889447
0.011783919598
0.0120100502513
0.0122361809045
0.0124623115578
0.0126884422111
0.0129145728643
0.0131407035176
0.0133668341709
0.0135929648241
0.0138190954774
0.0140452261307
0.0142713567839
0.0144974874372
0.0147236180905
0.0149497487437
0.015175879397
0.0154020100503
0.0156281407035
0.0158542713568
0.0160804020101
0.0163065326633
0.0165326633166
0.0167587939698
0.0169849246231
0.0172110552764
0.0174371859296
0.0176633165829
0.0178894472362
0.0181155778894
0.0183417085427
0.018567839196
0.0187939698492
0.0190201005025
0.0192462311558
0.019472361809
0.0196984924623
0.0199246231156
0.0201507537688
0.0203768844221
0.0206030150754
0.0208291457286
0.0210552763819
0.0212814070352
0.0215075376884
0.0217336683417
0.021959798995
0.0221859296482
0.0224120603015
0.0226381909548
0.022864321608
0.0230904522613
0.0233165829146
0.0235427135678
0.0237688442211
0.0239949748744
0.0242211055276
0.0244472361809
0.0246733668342
0.0248994974874
0.0251256281407
0.025351758794
0.0255778894472
0.0258040201005
0.0260301507538
0.026256281407
0.0264824120603
0.0267085427136
0.0269346733668
0.0271608040201
0.0273869346734
0.0276130653266
0.0278391959799
0.0280653266332
0.0282914572864
0.0285175879397
0.028743718593
0.0289698492462
0.0291959798995
0.0294221105528
0.029648241206
0.0298743718593
0.0301005025126
0.0303266331658
0.0305527638191
0.0307788944724
0.0310050251256
0.0312311557789
0.0314572864322
0.0316834170854
0.0319095477387
0.032135678392
0.0323618090452
0.0325879396985
0.0328140703518
0.033040201005
0.0332663316583
0.0334924623116
0.0337185929648
0.0339447236181
0.0341708542714
0.0343969849246
0.0346231155779
0.0348492462312
0.0350753768844
0.0353015075377
0.035527638191
0.0357537688442
0.0359798994975
0.0362060301508
0.036432160804
0.0366582914573
0.0368844221106
0.0371105527638
0.0373366834171
0.0375628140704
0.0377889447236
0.0380150753769
0.0382412060302
0.0384673366834
0.0386934673367
0.0389195979899
0.0391457286432
0.0393718592965
0.0395979899497
0.039824120603
0.0400502512563
0.0402763819095
0.0405025125628
0.0407286432161
0.0409547738693
0.0411809045226
0.0414070351759
0.0416331658291
0.0418592964824
0.0420854271357
0.0423115577889
0.0425376884422
0.0427638190955
0.0429899497487
0.043216080402
0.0434422110553
0.0436683417085
0.0438944723618
0.0441206030151
0.0443467336683
0.0445728643216
0.0447989949749
0.0450251256281
0.0452512562814
0.0454773869347
0.0457035175879
0.0459296482412
0.0461557788945
0.0463819095477
0.046608040201
0.0468341708543
0.0470603015075
0.0472864321608
0.0475125628141
0.0477386934673
0.0479648241206
0.0481909547739
0.0484170854271
0.0486432160804
0.0488693467337
0.0490954773869
0.0493216080402
0.0495477386935
0.0497738693467
0.05
"""




# In[50]:

cl_array = np.array([cls0, cls1, cls2, cls3, cls4, cls5, cls6, cls7, cls8, cls9, cls10,
                     cls11, cls12, cls13, cls14, cls15, cls16, cls17, cls18, cls19, cls20, 
                     cls21, cls22, cls23, cls24, cls25, cls26, cls27, cls28, cls29, cls30, 
                     cls31, cls32, cls33, cls34, cls35, cls36, cls37, cls38, cls39, cls40,
                     cls41, cls42, cls43, cls44, cls45, cls46, cls47, cls48, cls49, cls50,
                     cls51, cls52, cls53, cls54, cls55, cls56, cls57, cls58, cls59, cls60,
                     cls61, cls62, cls63, cls64, cls65, cls66, cls67, cls68, cls69, cls70,
                     cls71, cls72, cls73, cls74, cls75, cls76, cls77, cls78, cls79, cls80,
                     cls81, cls82, cls83, cls84, cls85, cls86, cls87, cls88, cls89, cls90,
                     cls91, cls92, cls93, cls94, cls95, cls96, cls97, cls98, cls99, cls100,
                     cls101, cls102, cls103, cls104, cls105, cls106, cls107, cls108, cls109, cls110,
                     cls111, cls112, cls113, cls114, cls115, cls116, cls117, cls118, cls119, cls120,
                     cls121, cls122, cls123, cls124, cls125, cls126, cls127, cls128, cls129, cls130,
                     cls131, cls132, cls133, cls134, cls135, cls136, cls137, cls138, cls139, cls140,
                     cls141, cls142, cls143, cls144, cls145, cls146, cls147, cls148, cls149, cls150,
                     cls151, cls152, cls153, cls154, cls155, cls156, cls157, cls158, cls159, cls160,
                     cls161, cls162, cls163, cls164, cls165, cls166, cls167, cls168, cls169, cls170,
                     cls171, cls172, cls173, cls174, cls175, cls176, cls177, cls178, cls179, cls180,
                     cls181, cls182, cls183, cls184, cls185, cls186, cls187, cls188, cls189, cls190,
                     cls191, cls192, cls193, cls194, cls195, cls196, cls197, cls198, cls199, cls200])




# In[51]:

print(cl_array.shape)


# In[52]:

f = "CAMB_cl_varyBaryon_lmax1100varyFeb2016.npy"

np.save(f, cl_array)

