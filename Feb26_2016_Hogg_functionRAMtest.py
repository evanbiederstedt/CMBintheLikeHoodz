
# coding: utf-8


import math
import matplotlib.pyplot as plt 
import numpy as np
import healpy as hp
import astropy as ap
import scipy
import scipy.io
from scipy.special import eval_legendre  ##special scipy function
import os

os.getcwd()
os.chdir("/Users/evanbiederstedt/downloads")

# For lmax = 1100, we must create an array of ell values, i.e. [0 1 2 3....1599 1600]
ell = np.arange(1101)
#print(ell)
# 
# Subtract the monopole and dipole, l=0, l=1
ellval = ell[2:]
#print(ellval)

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



# see script CAMB_vary_OmegaB_lmax1100_Feb2016.py

ff = "CAMB_cl_varyBaryon_lmax1100varyFeb2016.npy"

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




# Multiply to get:
#         [5/4pi*P_2(M) + 7/4pi*P_3(M) +...... + 65/4pi*P_32(M)]
#
# multiply PlMat by (2*l+1)/4pi, i.e. norm
norm_matrix = norm[:, None, None] * PlMat
# [5/4pi * P_2(M)  7/4pi * P_3(M) ....   65/4pi * P_32(M)]



"""

# define pixel-value arrays
mT = np.matrix(patch)     # mT.shape = (1, 768)
m = np.matrix(patch).T    # m.shape = (768, 1)
Npix2pi = (len(patch))*2*math.pi  # LF constant
print(mT.shape)
print(m.shape)
print(Npix2pi)
"""


#
# vary Omega_Baryon
# 
# Planck found \Omega_B = 0.02234
# GAVO simulated map set at \Omega_B = 0.04
# CAMB default below at ombh2=0.022
#
twohundred_samples = np.linspace(0.005, 0.05, num=200)
print(twohundred_samples)

#
# Step A: Set up one covariance matrix C=S+N at the "standard" set of cosmological parameters
#  pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
#


# 
# Multiply by 1e10 to avoid underflow
#

"""
tempp = patch*(1e10)
noise = noisepatch*(1e10)
cell = cls0*(1e10)
"""


fname_draw = "Feb12_2016_Multinormdraw.npy"
fname_noise = "Feb12_2016_noise.npy"

draws = np.load(fname_draw)
noise = np.load(fname_noise)



id_matrix = np.identity(len(noise))




"""
    
    id_matrix = np.identity(len(noise))
    Nij = noise**2 * id_matrix
    Cij = Sij + Nij
    print(Cij.shape)
    
"""

Npix2pi = (len(noise))*2*math.pi  # LF constant



#
# If we generate data from our model, then the likelihood MUST, on average, peak at the correct parameters.  
# 
# So generate (from the C=S+N matrix) 1000 patches at one set of cosmological parameters, 
# and compute the logLF as a function of one of those parameters, 
# and show that we have a peak at the right parameters.
#
#



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


#
# 'tempp' is 'draws' i.e. 
# draws = np.random.multivariate_normal(zero_mean, Cij, 1000)
# where draws.shape = (1000, 768), i.e. each 'draw' is vector (768,)
#




# draw_points = np.asarray([[LogLF_twoargs(cell_array[i], draws[j]) for i in range(40)] for j in range(1000)])




LogLF_draw0 = np.asarray([LogLF_twoargs(cell_array[0], draws[0]))
#print(draw_points.shape)
print("*******")
print(LogLF_draw0.shape)


f1 = "draw_points_100_RAMtest.npy"
np.save(f1, LogLF_draw0)
                          
del(LogLF_draw0)


