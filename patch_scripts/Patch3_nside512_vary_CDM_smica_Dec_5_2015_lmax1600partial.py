
# coding: utf-8
# December 5, 2015

# In[6]:

from __future__ import (division, print_function, absolute_import)


# In[7]:

import math
import matplotlib.pyplot as plt 
import numpy as np
import healpy as hp
import pyfits as pf
import astropy as ap
import os
from scipy.special import eval_legendre  ##special scipy function


# In[8]:

# http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.io.readsav.html
# http://www.astrobetter.com/blog/2009/11/24/read-idl-save-files-into-python/


# In[4]:

import scipy
import scipy.io



# In[6]:

import os
os.getcwd()
os.chdir("/Users/evanbiederstedt/downloads")





# In[11]:

# In[8]:

patch_file = scipy.io.readsav('listpix_patch3.sav')


# In[ ]:




# In[9]:

type(patch_file)


# In[10]:

arr3 = patch_file['listpix_patch3']
#print(arr3)


# In[11]:

type(arr3)


# In[12]:

print(len(arr3)) # pixels total 12476


# In[13]:

smica_map = "COM_CompMap_CMB-smica_2048_R1.20.fits"


# In[ ]:




# In[14]:

nside=512
npix = 12*(nside**2) #total number of pixels, npix
LMAX = ((2*nside)) #maximum l of the power spectrum C_l
heal_npix = hp.nside2npix(nside) # Healpix calculated npix

print("The total number of pixels is " + str(npix))
print("The maximum ell of the power spectrum C_l set to lmax = 2*nside " +str(LMAX))
print("Healpix tells me total number of pixels npix is equal to " + str(heal_npix))


# In[15]:

mapread_smica = hp.read_map(smica_map, field=0)
#hp.mollview(mapread_camb512)
#hp.mollview(mapread_smica)
print("CMB map, Noise map")
smica_noise = hp.read_map(smica_map, field=1)
#hp.mollview(smica_noise)


# In[16]:

print(mapread_smica[:20])
print(smica_noise[:20])


# In[17]:

smica512 = hp.pixelfunc.ud_grade(mapread_smica, 512)
noise512 = hp.pixelfunc.ud_grade(smica_noise, 512)
print(smica512[:20])
print(noise512[:20])


# In[18]:

print(len(smica512))
print(len(noise512))


# In[ ]:




# In[19]:

# rename array for convenience
tempval = smica512

# Data:
#     tempval      # the array of pixel values, (3145728,)


# In[20]:

print(len(tempval))
print(tempval.shape)
tempval[:10]


# In[21]:

#
# We only wish to use the pixels defined in our patch
# These pixel indices are listed in arr3 such that total number pixels total 12476
#
# arr3: this defines pixel indices within patch
# 
# To access pixel indices within array of CMB pixels, just use tempval[arr3]
#
patch=smica512[arr3]
noisepatch = noise512[arr3]


# In[22]:

print(len(patch))
print(len(noisepatch))


# In[23]:

print(patch[:30])
print(noisepatch[:30])


# In[ ]:




# In[12]:

# For lmax = 1600, we must create an array of ell values, i.e. [0 1 2 3....1599 1600]
ell = np.arange(1601)
#print(ell)
# 
# Subtract the monopole and dipole, l=0, l=1
ellval = ell[2:]
#print(ellval)


# In[ ]:




# In[13]:

PlM_50 = "cl_varyCDMlmax1600ptPlMat50.npy"
PlM_100 = "cl_varyCDMlmax1600ptPlMat100.npy"
PlM_150 = "cl_varyCDMlmax1600ptPlMat150.npy"


# In[14]:

data1 = np.load(PlM_50)
data2 = np.load(PlM_100)
data3 = np.load(PlM_150)


# In[15]:

print(data1.shape)
print(data2.shape)
print(data3.shape)


# In[16]:

type(data1)


# In[ ]:

ff = "CAMB_cl_varyCDMlmax1600.npy"

cell_array = np.load(ff)




# In[ ]:




# In[ ]:




# In[17]:


PlMat_total = np.concatenate((data1, data2, data3))


# In[18]:

PlMat_total.shape


# In[ ]:




# In[19]:

PlMat = PlMat_total


# In[20]:

PlMat[2]


# In[ ]:

# Step 3: (2*l +1)/4pi from l=2 to l=lmax
#          [5/4pi 7/4pi 9/4pi 11/4pi .... 65/4pi ]
norm = ((2*ellval + 1))/(4*math.pi)
print(len(ellval))
print(norm.shape)
print(norm[2])


# In[ ]:


# Step 4: multiply 
#         [5/4pi*P_2(M) + 7/4pi*P_3(M) +...... + 65/4pi*P_32(M)]
#
# multiply PlMat by (2*l+1)/4pi, i.e. norm
norm_matrix = norm[:, None, None] * PlMat
# [5/4pi * P_2(M)  7/4pi * P_3(M) ....   65/4pi * P_32(M)]


# In[ ]:

print(norm_matrix.shape)


# In[ ]:

print(PlMat.shape)


# In[ ]:

# Step 5: multiply by theoretical CAMB values, [C_2 C_3    C_31 C_32]
#         [5/4pi**C_2* P_2(M) + 7/4pi*C_3* P_3(M) +...... + 65/4pi*C_32* P_32(M)]



# In[ ]:

# define pixel-value arrays
mT = np.matrix(patch)     # mT.shape = (1, 3072)
m = np.matrix(patch).T    # m.shape = (3072, 1)
Npix2pi = (len(patch))*2*math.pi  # LF constant
print(mT.shape)
print(m.shape)
print(Npix2pi)




# In[ ]:

tempp = patch
noise = noisepatch


def LogLF(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms + logdetC[1] + Npix2pi

def modelfit(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms

def logdet(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return logdetC[1] 


def squaredLogLF(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = (noise**2) * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms + logdetC[1] + Npix2pi

def squared_modelfit(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = (noise**2) * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms

def squared_logdet(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = (noise**2) * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return logdetC[1] 



def noiselessLogLF(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij #+ Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms + logdetC[1] + Npix2pi

def noiselessmodelfit(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij #+ Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms

def noiselesslogdet(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij #+ Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return logdetC[1] 


# In[ ]:

def noiseonlyLogLF(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms + logdetC[1] + Npix2pi

def noiseonlymodelfit(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms

def noiseonlylogdet(cell):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return logdetC[1] 


# In[ ]:

forty_samples = np.linspace(0.015, 0.035, num=40)


# In[ ]:

logLF_40 = [LogLF(cell_array[i]) for i in range(40)]


# In[ ]:

modelfit_terms = [modelfit(cell_array[i])  for i in range(40)]


# In[ ]:

logdet_terms = [logdet(cell_array[i]) for i in range(40)]


# In[ ]:

sqlogLF_40 = [squaredLogLF(cell_array[i]) for i in range(40)]


# In[ ]:

sqmodelfit_terms = [squared_modelfit(cell_array[i]) for i in range(40)]


# In[ ]:

sqlogdet_terms = [squared_logdet(cell_array[i]) for i in range(40)]


# In[ ]:

noise_logLF = [noiselessLogLF(cell_array[i]) for i in range(40)]


# In[ ]:

noise_modelfits = [noiselessmodelfit(cell_array[i]) for i in range(40)]


# In[ ]:

noise_logdet = [noiselesslogdet(cell_array[i]) for i in range(40)]


# In[ ]:

onlynoise_logLF = [noiseonlyLogLF(cell_array[i]) for i in range(40)]


# In[ ]:

onlynoise_modelfits = [noiseonlymodelfit(cell_array[i]) for i in range(40)]


# In[ ]:

onlynoise_logdet = [noiseonlylogdet(cell_array[i]) for i in range(40)]


# In[ ]:




# In[ ]:

modelfit_terms


# In[ ]:

plt.plot(forty_samples, logLF_40)
plt.title("-2loglF ouput, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo1.png")

# In[ ]:

plt.plot(forty_samples, modelfit_terms)
plt.title("only model fit terms, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo2.png")


# In[ ]:

plt.plot(forty_samples, logdet_terms)
plt.title("only logdetC, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo3.png")


# In[ ]:

plt.plot(forty_samples, sqlogLF_40)
plt.title("squared noise, -2logLF, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo4.png")


# In[ ]:

plt.plot(forty_samples, sqmodelfit_terms)
plt.title("squared noise, model fit terms, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo5.png")

# In[ ]:

plt.plot(forty_samples, sqlogdet_terms)
plt.title("squared noise, logdet C terms, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo6.png")

# In[ ]:

plt.plot(forty_samples, onlynoise_logLF)
plt.title("Sij=0, -2logLF, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo7.png")

# In[ ]:

plt.plot(forty_samples, onlynoise_modelfits)
plt.title("Sij=0, model fit terms, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo8.png")

# In[ ]:

plt.plot(forty_samples, onlynoise_logdet)
plt.title("Sij=0, logdet C, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo9.png")

# In[ ]:




# In[ ]:

plt.plot(forty_samples, noise_logLF)
plt.title("Nij=0, -2logLF, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo10.png")


# In[ ]:

plt.plot(forty_samples, noise_modelfits)
plt.title("Nij=0, model fit terms, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo11.png")


# In[ ]:

plt.plot(forty_samples, noise_logdet)
plt.title("Nij=0, logdet C, SMICA Planck map")
plt.ylabel("-2logLF")
plt.xlabel("Omega_CDM")
plt.axvline(x = 0.12029, color = 'r')
plt.savefig("foo12.png")


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



