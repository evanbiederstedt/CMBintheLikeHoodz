
# coding: utf-8

# In[1]:

# 
# Planck found \Omega_CDM = 0.12029
#


# In[2]:

#
# First output 40 CAMB scalar outputs
#
# 0.075 to 0.1655
#


# In[3]:


from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower


# In[4]:

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


# In[6]:

forty_samples = np.linspace(0.075, 0.1655, 40)
print(forty_samples)


# In[7]:

#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.075, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
#pars.set_for_lmax(2000, lens_potential_accuracy=0)
#calculate results for these parameters
results = camb.get_results(pars)
#get dictionary of CAMB power spectra
powers = results.get_cmb_power_spectra(pars)
# plot the total lensed CMB power spectra versus unlensed, and fractional difference
totCL=powers['total']
unlensedCL=powers['unlensed_scalar']
print(totCL.shape)
# Python CL arrays are all zero based (starting at L=0), Note L=0,1 entries will be zero by default.
# The differenent CL are always in the order TT, EE, BB, TE (with BB=0 for unlensed scalar results).
ls = np.arange(totCL.shape[0])
print(ls)
print(unlensedCL.shape)


# In[8]:

print(unlensedCL.shape[0])


# In[9]:

#
# 0.075
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.075, mnu=0.06, omk=0, tau=0.06)
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

cls0 = unlencl[:,0][2:1601]
print(len(cls0))


# In[10]:

#
# 0.07732051
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.07732051, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls1 = unlencl[:,0][2:1601]


# In[11]:

#
# 0.07964103
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.07964103, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls2 = unlencl[:,0][2:1601]


# In[12]:

#
# 0.08196154
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.08196154, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls3 = unlencl[:,0][2:1601]


# In[13]:

#
# 0.08428205
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.08428205, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls4 = unlencl[:,0][2:1601]


# In[14]:

#
# 0.08660256
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.08660256, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls5 = unlencl[:,0][2:1601]


# In[15]:

#
# 0.08892308
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.08892308, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls6 = unlencl[:,0][2:1601]


# In[16]:

#
# 0.09124359
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.09124359, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls7 = unlencl[:,0][2:1601]


# In[17]:

#
# 0.0935641
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.0935641, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls8 = unlencl[:,0][2:1601]


# In[18]:

#
# 0.09588462
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.09588462, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls9 = unlencl[:,0][2:1601]


# In[19]:

#
# 0.09820513
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.09820513, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls10 = unlencl[:,0][2:1601]


# In[20]:

#
# 0.10052564
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.10052564, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls11 = unlencl[:,0][2:1601]


# In[21]:

#
# 0.10284615
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.10284615, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls12 = unlencl[:,0][2:1601]


# In[22]:

#
# 0.10516667
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.10516667, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls13 = unlencl[:,0][2:1601]


# In[23]:

#
# 0.10748718
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.10748718, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls14 = unlencl[:,0][2:1601]


# In[24]:

#
# 0.10980769
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.10980769, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls15 = unlencl[:,0][2:1601]


# In[25]:

#
# 0.11212821
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.11212821, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls16 = unlencl[:,0][2:1601]


# In[26]:

#
# 0.11444872
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.11444872, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls17 = unlencl[:,0][2:1601]


# In[27]:

#
# 0.11676923
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.11676923, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls18 = unlencl[:,0][2:1601]


# In[28]:

#
# 0.11908974
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.11908974, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls19 = unlencl[:,0][2:1601]


# In[29]:

#
# 0.12141026
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12141026, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls20 = unlencl[:,0][2:1601]


# In[30]:

#
# 0.12373077
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12373077, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls21 = unlencl[:,0][2:1601]


# In[31]:

#
# 0.12605128
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12605128, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls22 = unlencl[:,0][2:1601]


# In[32]:

#
# 0.12837179
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.12837179, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls23 = unlencl[:,0][2:1601]


# In[33]:

#
# 0.13069231
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.13069231, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls24 = unlencl[:,0][2:1601]


# In[34]:

#
# 0.13301282
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.13301282, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls25 = unlencl[:,0][2:1601]


# In[35]:

#
# 0.13533333
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.13533333, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls26 = unlencl[:,0][2:1601]


# In[36]:

#
# 0.13765385
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.13765385, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls27 = unlencl[:,0][2:1601]


# In[37]:

#
# 0.1399743
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.1399743, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls28 = unlencl[:,0][2:1601]


# In[38]:

#
# 0.14229487
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.14229487, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls29 = unlencl[:,0][2:1601]


# In[39]:

#
# 0.14461538
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.14461538, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls30 = unlencl[:,0][2:1601]


# In[40]:

#
# 0.1469359
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.1469359, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls31 = unlencl[:,0][2:1601]


# In[41]:

#
# 0.14925641
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.14925641, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls32 = unlencl[:,0][2:1601]


# In[42]:

#
# 0.15157692
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.15157692, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls33 = unlencl[:,0][2:1601]


# In[43]:

#
# 0.15389744
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.15389744, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls34 = unlencl[:,0][2:1601]


# In[44]:

#
# 0.15621795
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.15621795, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls35 = unlencl[:,0][2:1601]


# In[45]:

#
# 0.15853846
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.15853846, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls36 = unlencl[:,0][2:1601]


# In[46]:

#
# 0.16085897
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.16085897, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls37 = unlencl[:,0][2:1601]


# In[47]:

#
# 0.16317949
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.16317949, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls38 = unlencl[:,0][2:1601]


# In[48]:

#
# 0.1655
#
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.1655, mnu=0.06, omk=0, tau=0.06)
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
# plot = plt.plot(ls, cl_075, color='r')
cls39 = unlencl[:,0][2:1601]


# In[ ]:




# In[50]:

cl_array = np.array([cls0, cls1, cls2, cls3, cls4, cls5, cls6, cls7, cls8, cls9, cls10,
                     cls11, cls12, cls13, cls14, cls15, cls16, cls17, cls18, cls19, cls20, 
                     cls21, cls22, cls23, cls24, cls25, cls26, cls27, cls28, cls29, cls30, 
                     cls31, cls32, cls33, cls34, cls35, cls36, cls37, cls38, cls39])


# In[51]:

print(cl_array.shape)


# In[52]:

f = "CAMB_cl_varyCDMlmax1600.npy"

np.save(f, cl_array)

