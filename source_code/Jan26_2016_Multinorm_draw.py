
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib inline')


# In[2]:

from __future__ import (division, print_function, absolute_import)


import math
import matplotlib.pyplot as plt 
import numpy as np
import healpy as hp
import pyfits as pf
import astropy as ap
import scipy
import scipy.io
from scipy.special import eval_legendre  ##special scipy function
import os

os.getcwd()
os.chdir("/Users/evanbiederstedt/downloads")

patch_file = scipy.io.readsav('listpix_patch3.sav')  # import "patch indices" from IDL
arr3 = patch_file['listpix_patch3']

# print(len(arr3)) # pixels total 768


# In[3]:

gavo_map = "GAVO_nside512_omegab_004.fits"
# this is GAVO simulated data, http://gavo.mpa-garching.mpg.de/planck/


# In[4]:

# set Healpix parameters: nside, lmax

nside=512
npix = 12*(nside**2) #total number of pixels, npix
LMAX = ((2*nside)) #maximum l of the power spectrum C_l
heal_npix = hp.nside2npix(nside) # Healpix calculated npix

print("The total number of pixels is " + str(npix))
print("The maximum ell of the power spectrum C_l set to lmax = 2*nside " +str(LMAX))
print("Healpix tells me total number of pixels npix is equal to " + str(heal_npix))


# In[5]:

gavo_signal = hp.read_map(gavo_map, field=0)  # field=0 is temperature
print("CMB map, Noise map")
gavo_noise = hp.read_map(gavo_map, field=1)    # field=1 is noise


# In[6]:

# simply show maps
hp.mollview(gavo_signal)
hp.mollview(gavo_noise)


# In[7]:

#########################
#
# We only wish to use the pixels defined in our patch
# These pixel indices are listed in arr3 such that total number pixels total 12476
#
# arr3: this defines pixel indices within patch
# 
# To access pixel indices within array of CMB pixels, just use tempval[arr3]
#
#########################

patch =  gavo_signal[arr3]
noisepatch = gavo_noise[arr3]

print(len(patch))
print(len(noisepatch))


# In[8]:

#######################
#
# The log-likelihood
#
# -2lnL \propto m^T C^-1 m + ln det C + N ln (2pi)
#
# First term, m^T C^-1 m is the "model fit term"
# Second term, lndetC is the "complexity penalty"
# Third term, N ln 2pi, a constant
#
# m = tempval
# C = Sij
#
#######################



# Here, m = patch


#####################
#
# Dotproductmatrix
#
#####################


# Next, create the matrix, n_i /cdot n_j
# solely using Healpy routines, i.e. taking the dot product of the vectors
# The result is "dotproductmatrix"

## healpy.pixelfunc.pix2vec(nside, ipix, nest=False)
## 
## will give three arrays
## arrays of all x values, all y values, all z values
## RING scheme default
# len()=3
# type()=tuple


vecval = hp.pix2vec(nside, arr3) #Nside = 512, type()=tuple


len(vecval)


vecvalx = vecval[0] #len() = 12476
vecvaly = vecval[1]
vecvalz = vecval[2]



# First arrange arrays vertically
# numpy.vstack = Stack arrays in sequence vertically (row wise), input sequence of arrays
totalvecval = np.vstack((vecvalx, vecvaly, vecvalz)) #type()=numpy.ndarray


print(totalvecval)
print("*******")

trans = totalvecval.T #transpose

print(trans)
print("*******")

dotproductmatrix = trans.dot(totalvecval) #take the dot product
print(dotproductmatrix.shape) # = (npix, npix) = (12476, 12476)
print("*******")

print(dotproductmatrix[:10][:10])
print("*******")


# In[9]:

# you can check all negative eigenvalues are VERY small, order e-15
#
# eigenvals = np.linalg.eigvals(dotproductmatrix) 
#
# for i in eigenvals: 
#    if i < 0:
#        print(i)
#


# In[10]:

#
# The following procedure is for the angular power spectrum, C^th_ell
# However, we are using some cosmological parameter, /alpha
#
#
# =========================================================
# =========================================================
#
# \Sum_l (2*l + 1)/4pi C^th_l P_l (dotproductmatrix)
# sum from l=2 to l=lmax
#
# arrays l = [2 3 4 .... lmax]
#        C_l = [C_2 C_3 .... C_lmax]
#
# The correct way to do the summation (for e.g. lmax=32):
# 
# Step 1: calculate the matrix
#            M = dotproductmatrix
#
# Step 2: evaluate the function P_l(x) for each entry of the matrix
#         OUTPUT: [P_2(M) P_3(M) P_4(M) .... P_lmax(M) ]
#
# Step 3: (2*l +1)/4pi from l=2 to l=lmax
#          [5/4pi 7/4pi 9/4pi 11/4pi .... 65/4pi ]
#
# Step 4: multiply 
#         [5/4pi*P_2(M) + 7/4pi*P_3(M) +...... + 65/4pi*P_32(M)]
#
#
# Step 5: multiply by theoretical CAMB values, [C_2 C_3    C_31 C_32]
#         [5/4pi**C_2* P_2(M) + 7/4pi*C_3* P_3(M) +...... + 65/4pi*C_32* P_32(M)]
#
# Step 6: This is an array of S_ij for each theory C_l, l=2 to l=32
#         
#
#
# =========================================================
# =========================================================


# In[11]:

# For lmax = 1100, we must create an array of ell values, i.e. [0 1 2 3....1599 1600]
ell = np.arange(1101)
#print(ell)
# 
# Subtract the monopole and dipole, l=0, l=1
ellval = ell[2:]
#print(ellval)


# In[12]:

#
# Vary Baryon, patch 3
#

PlM_50 = "cl_varyBaryonlmax1100patch3PlMat50.npy"
PlM_100 = "cl_varyBaryonlmax1100patch3PlMat100.npy"
PlM_150 = "cl_varyBaryonlmax1100patch3PlMat150.npy"



data1 = np.load(PlM_50)
data2 = np.load(PlM_100)
data3 = np.load(PlM_150)




print(data1.shape)
print(data2.shape)
print(data3.shape)



# see script CAMB_vary_OmegaBaryon_Dec_9_lmax1100.py

ff = "CAMB_cl_varyBaryonlmax1100vary.npy"

cell_array_loaded = np.load(ff)

cell_array = cell_array_loaded*(1e10)

print(cell_array.shape)


# In[13]:

PlMat_total = np.concatenate((data1, data2, data3))   # this is P_2(M), P_3(M), ..., P_lmax (M)

PlMat_total.shape


PlMat = PlMat_total


# Step 3: (2*l +1)/4pi from l=2 to l=lmax
#          [5/4pi 7/4pi 9/4pi 11/4pi .... 65/4pi ]
norm = ((2*ellval + 1))/(4*math.pi)
print(len(ellval))
print("******")
print(norm.shape)
print("*****")
print(PlMat.shape)


# In[14]:

#
# Check sizes of eigenvalues
# 
# eigenvals_PlMat = np.linalg.eigvals(PlMat)
#
# eigenvals_PlMat.shape
#
# for row in eigenvals_PlMat.real:
#     for i in row:
#         if (i<0 and i<-1e-10):
#             print(i)
#           
# -1.00032390644e-10
# -1.00260768658e-10
# -1.00488591982e-10
# -1.00715150533e-10
#
# The rest are smaller than this, on the order -1e-15
#


# In[15]:

# Multiply to get:
#         [5/4pi*P_2(M) + 7/4pi*P_3(M) +...... + 65/4pi*P_32(M)]
#
# multiply PlMat by (2*l+1)/4pi, i.e. norm
norm_matrix = norm[:, None, None] * PlMat
# [5/4pi * P_2(M)  7/4pi * P_3(M) ....   65/4pi * P_32(M)]


# In[16]:

#
# Check sizes of eigenvalues
#
# eigenvals_norm_matrix = np.linalg.eigvals(norm_matrix) 
#
# for row in eigenvals_PlMat.real:
#     for i in row:
#         if (i<0 and i<-1e-10):
#             print(i)
#           
# -1.02747972205e-08
# -1.03136456948e-08
# -1.03523982409e-08
# -1.03911221835e-08
#
# The rest are smaller than this, on the order -1e-15
#   


# In[17]:

# define pixel-value arrays
mT = np.matrix(patch)     # mT.shape = (1, 768)
m = np.matrix(patch).T    # m.shape = (768, 1)
Npix2pi = (len(patch))*2*math.pi  # LF constant
print(mT.shape)
print(m.shape)
print(Npix2pi)


# In[18]:

#
# vary Omega_Baryon
# 
# Planck found \Omega_B = 0.02234
# GAVO simulated map set at \Omega_B = 0.04
# CAMB default below at ombh2=0.022
#
forty_samples = np.linspace(0.005, 0.10, num=40)
print(forty_samples)


# In[19]:

#
# Step A: Set up one covariance matrix C=S+N at the "standard" set of cosmological parameters
#  pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
#


# In[ ]:




# In[ ]:




# In[20]:

import camb
from camb import model, initialpower
#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
pars.set_for_lmax(1100, lens_potential_accuracy=0);

#calculate results for these parameters
results = camb.get_results(pars)

#get dictionary of CAMB power spectra
powers =results.get_cmb_power_spectra(pars)
unlencl = powers['unlensed_scalar']

cls0 = unlencl[:,0][2:1101]
print(len(cls0))


# In[21]:

# 
# Multiply by 1e10 to avoid underflow
#

tempp = patch*(1e10)
noise = noisepatch*(1e10)
cell = cls0*(1e10)


# In[22]:

# Construct the matrix using definition above, def LogLF()
#
# norm_matrix is (2*l+1)/4pi * P_ell(Mat)
# elementwise LF is (2*l+1)/4pi * C^th_ell * P_ell(Mat)
#
# # 5/4pi * C^th_2 * P_2(Mat), 7/4pi * C^th_3 * P_3(Mat), 9/4pi * C^th_4 * P_4(Mat), ..., 
CellPellM = cell[:, None, None] * norm_matrix


# In[23]:

print(CellPellM.shape)


# In[24]:

#
# eigenvals_CellPellM = np.linalg.eigvals(CellPellM) 
# 
# for row in eigenvals_CellPellM.real:
#     for i in row:
#         if (i<0 and i<-1e-8):
#             print(i)
#             
# -1.44546933082e-08
# -1.12224623041e-08
# -1.52809195714e-08
#


# In[25]:

# Do the summation
Sij =  np.sum(CellPellM, axis=0) # now one matrix


# In[26]:

#
# eigenvals_Sij = np.linalg.eigvals(Sij) 
#
# for i in eigenvals_Sij.real:
#     if (i<0 and i<-1e-6):
#         print(i)
#        
# -5.21725495676e-06
# -5.35688881031e-06
# -5.33350330933e-06
# -4.91069663263e-06
#


# In[27]:

print(Sij.shape)


# In[28]:

id_matrix = np.identity(len(noise))
Nij = noise**2 * id_matrix
Cij = Sij + Nij
print(Cij.shape)


# In[29]:

# 
# With squared noise, noise**2, the matrix is positive semidefinite. Otherwise, it fails DRAMATICALLY
# 
# eigenvals_Cij = np.linalg.eigvals(Cij) 
# 
# for i in eigenvals_Cij.real:
#     if (i<0):
#         print(i)
#
#
# Otherwise, with noise not squared, it fails DRAMATICALLY:
#
# A sample of eigenvalues:
# -1366.09930314
# -34.3949953695
# -186.967111398
# -396.917101696
# -261.826549045
# -379.420389554
# -234.189298314
# -531.082474565
# -1573.04857094
# -78.087152717
# -79.6528937648
#
#


# In[30]:

print("Parameters set for covariance matrix Cij are H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06")
print("*********")
print(Cij)


# In[31]:

Cij.shape


# In[32]:

len(noise)


# In[33]:

#
#
########################
# Code test for Hogg: 
########################
#
# Step (1): generate 1000 patches at a set of cosmological parameters
#
zero_mean = np.zeros(len(noise))
draws = np.random.multivariate_normal(zero_mean, Cij, 1000)
#
# each 'draws'  is a "data point" , i.e. 
# a 768-vector which is an intensity field that the CAMB model plus noise could have generated
# 


# In[34]:

print(zero_mean.shape)
print(draws.shape)
print(draws[0].shape)


# In[35]:

#
# If we generate data from our model, then the likelihood MUST, on average, peak at the correct parameters.  
# 
# So generate (from the C=S+N matrix) 1000 patches at one set of cosmological parameters, 
# and compute the logLF as a function of one of those parameters, 
# and show that we have a peak at the right parameters.
#
#


# In[36]:

#
# LogLF function of cell, i.e. theoretical Boltzmann values C_ell and tempp, i.e. pixel values of the patch
#
def LogLF_twoargs(cell, tempp):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms + logdetC[1] + Npix2pi


# In[37]:

def modelfit_twoargs(cell, tempp):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
     #logdetC = np.linalg.slogdet(Cij)
    return model_fit_terms


# In[38]:

def logdet_twoargs(cell, tempp):
    # norm_matrix is (2*l+1)/4pi * P_ell(Mat)
    CellPellM = cell[:, None, None] * norm_matrix # elementwise (2*l+1)/4pi * C^th_ell * P_ell(Mat)
    Sij =  np.sum(CellPellM, axis=0) # now one matrix
    id_matrix = np.identity(len(tempp))
    Nij = noise * id_matrix
    Cij = Sij + Nij
    model_fit_terms = np.array([np.dot(tempp.T , (np.linalg.solve(Cij, tempp)) )])
    logdetC = np.linalg.slogdet(Cij)
    return logdetC[1] 


# In[39]:

#
# set 'cell' with 'cell_array', i.e. Boltzmann code values, i.e. values from "CAMB_cl_varyBaryonlmax1100vary.npy", 
# from Python3.4 script CAMB_vary_OmegaBaryon_Dec_9_lmax1100.py
# 
# These are Boltzmann code results of varying Omega_baryon from 0.005 to 0.05
# forty_samples = np.linspace(0.005, 0.05, num=40)
# i.e. 
# pars.set_cosmology(H0=67.5, ombh2=0.005, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
# 0.005
# pars.set_cosmology(H0=67.5, ombh2=0.00615385, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
# 0.00615385
# pars.set_cosmology(H0=67.5, ombh2=0.00730769, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
# 0.00730769
# etc. etc. 
#


# In[40]:

#
# 'tempp' is 'draws' i.e. 
# draws = np.random.multivariate_normal(zero_mean, Cij, 1000)
# where draws.shape = (1000, 768), i.e. each 'draw' is vector (768,)
#


# In[41]:

# draw_points = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in range(1000)])


# In[42]:

LogLF_draw0 = np.asarray([LogLF_twoargs(cell_array[i], draws[0]) for i in range(40)])


# In[43]:

LogLF_draw012 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in range(3)])


# In[44]:

#print(draw_points.shape)
print("*******")
print(LogLF_draw0.shape)
print("*******")
print(LogLF_draw012.shape)


# In[45]:

cell_array.shape


# In[46]:

cell_array[0].shape


# In[ ]:




# In[47]:

draws.shape


# In[48]:

draws[0].shape


# In[49]:

LogLF_draw01234 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in range(5)])


# In[50]:

LogLF_draw01234.shape


# In[ ]:




# In[51]:

LogLF_draw01234[0].shape


# In[52]:

LogLF_draw01234[0][0].shape


# In[53]:

draw_points_100 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(0,100)])


# In[54]:

draw_points_100.shape


# In[ ]:

f1 = "draw_points_100.npy"
np.save(f1, draw_points_100)


# In[ ]:

draw_points_200 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(100,200)])


# In[ ]:

f2 = "draw_points_200.npy"
np.save(f2, draw_points_200)


# In[ ]:

draw_points_300 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(200,300)])


# In[ ]:

f3 = "draw_points_300.npy"
np.save(f3, draw_points_300)


# In[ ]:

draw_points_400 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(300,400)])


# In[ ]:

f4 = "draw_points_400.npy"
np.save(f4, draw_points_400)


# In[ ]:

draw_points_500 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(400,500)])


# In[ ]:

f5 = "draw_points_500.npy"
np.save(f5, draw_points_500)


# In[ ]:

draw_points_600 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(500,600)])


# In[ ]:

f6 = "draw_points_600.npy"
np.save(f6, draw_points_600)


# In[ ]:

draw_points_700 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(600,700)])


# In[ ]:

f7 = "draw_points_700.npy"
np.save(f7, draw_points_700)


# In[ ]:

draw_points_800 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(700,800)])


# In[ ]:

f8 = "draw_points_800.npy"
np.save(f8, draw_points_800)


# In[ ]:

draw_points_900 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(800,900)])


# In[ ]:

f9 = "draw_points_900.npy"
np.save(f9, draw_points_900)


# In[ ]:

draw_points_1000 = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in np.arange(900,1000)])


# In[ ]:

f10 = "draw_points_1000.npy"
np.save(f10, draw_points_1000)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



