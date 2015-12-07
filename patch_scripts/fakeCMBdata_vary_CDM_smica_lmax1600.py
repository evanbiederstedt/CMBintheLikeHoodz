
# coding: utf-8
# December 2-3, 2015
#
# copied from  cp /Users/evanbiederstedt/Downloads/Patch3_nside512_vary_CDM_smica_lmax1600.py
#
#

# In[1]:

from __future__ import (division, print_function, absolute_import)


# In[2]:

import math
import matplotlib.pyplot as plt 
import numpy as np
import healpy as hp
import pyfits as pf
import astropy as ap
import os
from scipy.special import eval_legendre  ##special scipy function


# In[3]:

# http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.io.readsav.html
# http://www.astrobetter.com/blog/2009/11/24/read-idl-save-files-into-python/


# In[4]:

import scipy


# In[5]:

#
# scipy.io.readsav
#
# scipy.io.readsav(file_name, idict=None, python_dict=False, uncompressed_file_name=None, verbose=False)[source]
#
# Read an IDL .sav file
#
#


# In[6]:

import os
os.getcwd()
os.chdir("/Users/evanbiederstedt/downloads")


import scipy.io


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

print(len(arr3)) # pixels total 768


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

########
########
# Here is where we create fake CMB data
#
#

mu, sigma = 0, 25 # mean and standard deviation


patch =  np.random.normal(mu, sigma, 768)
noisepatch = noise512[arr3]


# In[22]:

print(len(patch))
print(len(noisepatch))


# In[23]:

print(patch[:30])
print(noisepatch[:30])


# In[24]:

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


# In[25]:

m = patch


# In[26]:

# Next, create the matrix, n_i /cdot n_j
# solely using Healpy routines, i.e. taking the dot product of the vectors
# The result is "dotproductmatrix"


# In[27]:

npix


# In[28]:

nside


# In[29]:

## healpy.pixelfunc.pix2vec(nside, ipix, nest=False)
## 
## will give three arrays
## arrays of all x values, all y values, all z values
## RING scheme default
# len()=3
# type()=tuple


# In[30]:

vecval = hp.pix2vec(nside, arr3) #Nside = 512, type()=tuple


# In[31]:

len(vecval)


# In[32]:

vecvalx = vecval[0] #len() = 12476
vecvaly = vecval[1]
vecvalz = vecval[2]


# In[33]:

# First arrange arrays vertically
# numpy.vstack = Stack arrays in sequence vertically (row wise), input sequence of arrays
totalvecval = np.vstack((vecvalx, vecvaly, vecvalz)) #type()=numpy.ndarray


# In[34]:

trans = totalvecval.T #transpose


# In[35]:

dotproductmatrix = trans.dot(totalvecval) #take the dot product
print(dotproductmatrix.shape) # = (npix, npix) = (12476, 12476)
# type(dotproductmatrix) = np.ndarray


# In[36]:

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
# The correct way to do the summation:
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


# In[37]:

#print(dotproductmatrix)


# In[ ]:




# In[38]:

# For lmax = 1600, we must create an array of ell values, i.e. [0 1 2 3....1599 1600]
ell = np.arange(1601)
#print(ell)
# 
# Subtract the monopole and dipole, l=0, l=1
ellval = ell[2:]
#print(ellval)


# In[39]:

# The correct way to do the summation:
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


# In[40]:

dotproductmatrix.shape


# In[41]:

# Step 1: calculate the matrix
M = dotproductmatrix


# In[42]:

# Step 2: evaluate the function P_l(x) for each entry of the matrix
#         OUTPUT: [P_2(M) P_3(M) P_4(M) .... P_lmax(M) ]


# In[43]:

# CODE BOTTLENECK!
# 
# Evaluate Legendre from l=2 to l=lmax for each matrix entry
# [P_2(M) P_3(M) P_4(M) .... P_lmax(M) ]
# 
# WITHOUT BROADCASTING, one would do something like 
# PlMat = []
# for i in ellval:
#    PlMat.append( eval_legendre(i, dotproductmatrix) )
#
#
# With broadcasting, we use
# PlMat = eval_legendre(ellval[:, None, None], dotproductmatrix)
# PlMat = [P_2(M) P_3(M) P_4(M) .... P_lmax(M) ]
# PlMat is an array, len()=31 of 31 3072 by 3072 matrices
# PlMat.shape = (31, 3072, 3072)


# In[44]:

#This doesn't run for lmax=512
#So, split 'ellval' into ten arrays and then sum afterwards

splitell = np.array_split(ellval, 150)
splitell[0]


# In[45]:

PlMat1 = eval_legendre(splitell[0][:, None, None], dotproductmatrix)


# In[46]:

PlMat2 = eval_legendre(splitell[1][:, None, None], dotproductmatrix)


# In[47]:

PlMat3 = eval_legendre(splitell[2][:, None, None], dotproductmatrix)


# In[48]:

PlMat4 = eval_legendre(splitell[3][:, None, None], dotproductmatrix)


# In[49]:

PlMat5 = eval_legendre(splitell[4][:, None, None], dotproductmatrix)


# In[50]:

PlMat6 = eval_legendre(splitell[5][:, None, None], dotproductmatrix)


# In[51]:

PlMat7 = eval_legendre(splitell[6][:, None, None], dotproductmatrix)


# In[52]:

PlMat8 = eval_legendre(splitell[7][:, None, None], dotproductmatrix)


# In[53]:

PlMat9 = eval_legendre(splitell[8][:, None, None], dotproductmatrix)


# In[54]:

PlMat10 = eval_legendre(splitell[9][:, None, None], dotproductmatrix)


# In[55]:

PlMat11 = eval_legendre(splitell[10][:, None, None], dotproductmatrix)


# In[56]:

PlMat12 = eval_legendre(splitell[11][:, None, None], dotproductmatrix)


# In[57]:

PlMat13 = eval_legendre(splitell[12][:, None, None], dotproductmatrix)


# In[58]:

PlMat14 = eval_legendre(splitell[13][:, None, None], dotproductmatrix)


# In[59]:

PlMat15 = eval_legendre(splitell[14][:, None, None], dotproductmatrix)


# In[60]:

PlMat16 = eval_legendre(splitell[15][:, None, None], dotproductmatrix)


# In[61]:

PlMat17 = eval_legendre(splitell[16][:, None, None], dotproductmatrix)


# In[62]:

PlMat18 = eval_legendre(splitell[17][:, None, None], dotproductmatrix)


# In[63]:

PlMat19 = eval_legendre(splitell[18][:, None, None], dotproductmatrix)


# In[64]:

PlMat20 = eval_legendre(splitell[19][:, None, None], dotproductmatrix)


# In[65]:

PlMat21 = eval_legendre(splitell[20][:, None, None], dotproductmatrix)


# In[66]:

PlMat22 = eval_legendre(splitell[21][:, None, None], dotproductmatrix)


# In[67]:

PlMat23 = eval_legendre(splitell[22][:, None, None], dotproductmatrix)


# In[68]:

PlMat24 = eval_legendre(splitell[23][:, None, None], dotproductmatrix)


# In[69]:

PlMat25 = eval_legendre(splitell[24][:, None, None], dotproductmatrix)


# In[70]:

PlMat26 = eval_legendre(splitell[25][:, None, None], dotproductmatrix)


# In[71]:

PlMat27 = eval_legendre(splitell[26][:, None, None], dotproductmatrix)


# In[72]:

PlMat28 = eval_legendre(splitell[27][:, None, None], dotproductmatrix)


# In[73]:

PlMat29 = eval_legendre(splitell[28][:, None, None], dotproductmatrix)


# In[74]:

PlMat30 = eval_legendre(splitell[29][:, None, None], dotproductmatrix)


# In[75]:

PlMat31 = eval_legendre(splitell[30][:, None, None], dotproductmatrix)


# In[76]:

PlMat32 = eval_legendre(splitell[31][:, None, None], dotproductmatrix)


# In[77]:

PlMat33 = eval_legendre(splitell[32][:, None, None], dotproductmatrix)


# In[78]:

PlMat34 = eval_legendre(splitell[33][:, None, None], dotproductmatrix)


# In[79]:

PlMat35 = eval_legendre(splitell[34][:, None, None], dotproductmatrix)


# In[80]:

PlMat36 = eval_legendre(splitell[35][:, None, None], dotproductmatrix)


# In[81]:

PlMat37 = eval_legendre(splitell[36][:, None, None], dotproductmatrix)


# In[82]:

PlMat38 = eval_legendre(splitell[37][:, None, None], dotproductmatrix)


# In[83]:

PlMat39 = eval_legendre(splitell[38][:, None, None], dotproductmatrix)


# In[84]:

PlMat40 = eval_legendre(splitell[39][:, None, None], dotproductmatrix)


# In[85]:

PlMat41 = eval_legendre(splitell[40][:, None, None], dotproductmatrix)


# In[86]:

PlMat42 = eval_legendre(splitell[41][:, None, None], dotproductmatrix)


# In[87]:

PlMat43 = eval_legendre(splitell[42][:, None, None], dotproductmatrix)


# In[88]:

PlMat44 = eval_legendre(splitell[43][:, None, None], dotproductmatrix)


# In[89]:

PlMat45 = eval_legendre(splitell[44][:, None, None], dotproductmatrix)


# In[90]:

PlMat46 = eval_legendre(splitell[45][:, None, None], dotproductmatrix)


# In[91]:

PlMat47 = eval_legendre(splitell[46][:, None, None], dotproductmatrix)


# In[92]:

PlMat48 = eval_legendre(splitell[47][:, None, None], dotproductmatrix)


# In[93]:

PlMat49 = eval_legendre(splitell[48][:, None, None], dotproductmatrix)


# In[94]:

PlMat50 = eval_legendre(splitell[49][:, None, None], dotproductmatrix)


# In[95]:

PlMat51 = eval_legendre(splitell[50][:, None, None], dotproductmatrix)


# In[96]:

PlMat52 = eval_legendre(splitell[51][:, None, None], dotproductmatrix)


# In[97]:

PlMat53 = eval_legendre(splitell[52][:, None, None], dotproductmatrix)


# In[98]:

PlMat54 = eval_legendre(splitell[53][:, None, None], dotproductmatrix)


# In[99]:

PlMat55 = eval_legendre(splitell[54][:, None, None], dotproductmatrix)


# In[100]:

PlMat56 = eval_legendre(splitell[55][:, None, None], dotproductmatrix)


# In[101]:

PlMat57 = eval_legendre(splitell[56][:, None, None], dotproductmatrix)


# In[102]:

PlMat58 = eval_legendre(splitell[57][:, None, None], dotproductmatrix)


# In[103]:

PlMat59 = eval_legendre(splitell[58][:, None, None], dotproductmatrix)


# In[104]:

PlMat60 = eval_legendre(splitell[59][:, None, None], dotproductmatrix)


# In[105]:

PlMat61 = eval_legendre(splitell[60][:, None, None], dotproductmatrix)


# In[106]:

PlMat62 = eval_legendre(splitell[61][:, None, None], dotproductmatrix)


# In[107]:

PlMat63 = eval_legendre(splitell[62][:, None, None], dotproductmatrix)


# In[108]:

PlMat64 = eval_legendre(splitell[63][:, None, None], dotproductmatrix)


# In[109]:

PlMat65 = eval_legendre(splitell[64][:, None, None], dotproductmatrix)


# In[110]:

PlMat66 = eval_legendre(splitell[65][:, None, None], dotproductmatrix)


# In[111]:

PlMat67 = eval_legendre(splitell[66][:, None, None], dotproductmatrix)


# In[112]:

PlMat68 = eval_legendre(splitell[67][:, None, None], dotproductmatrix)


# In[113]:

PlMat69 = eval_legendre(splitell[68][:, None, None], dotproductmatrix)


# In[114]:

PlMat70 = eval_legendre(splitell[69][:, None, None], dotproductmatrix)


# In[115]:

PlMat71 = eval_legendre(splitell[70][:, None, None], dotproductmatrix)


# In[116]:

PlMat72 = eval_legendre(splitell[71][:, None, None], dotproductmatrix)


# In[117]:

PlMat73 = eval_legendre(splitell[72][:, None, None], dotproductmatrix)


# In[118]:

PlMat74 = eval_legendre(splitell[73][:, None, None], dotproductmatrix)


# In[119]:

PlMat75 = eval_legendre(splitell[74][:, None, None], dotproductmatrix)


# In[120]:

PlMat76 = eval_legendre(splitell[75][:, None, None], dotproductmatrix)


# In[121]:

PlMat77 = eval_legendre(splitell[76][:, None, None], dotproductmatrix)


# In[122]:

PlMat78 = eval_legendre(splitell[77][:, None, None], dotproductmatrix)


# In[123]:

PlMat79 = eval_legendre(splitell[78][:, None, None], dotproductmatrix)


# In[124]:

PlMat80 = eval_legendre(splitell[79][:, None, None], dotproductmatrix)


# In[125]:

PlMat81 = eval_legendre(splitell[80][:, None, None], dotproductmatrix)


# In[126]:

PlMat82 = eval_legendre(splitell[81][:, None, None], dotproductmatrix)


# In[127]:

PlMat83 = eval_legendre(splitell[82][:, None, None], dotproductmatrix)


# In[128]:

PlMat84 = eval_legendre(splitell[83][:, None, None], dotproductmatrix)


# In[129]:

PlMat85 = eval_legendre(splitell[84][:, None, None], dotproductmatrix)


# In[130]:

PlMat86 = eval_legendre(splitell[85][:, None, None], dotproductmatrix)


# In[131]:

PlMat87 = eval_legendre(splitell[86][:, None, None], dotproductmatrix)


# In[132]:

PlMat88 = eval_legendre(splitell[87][:, None, None], dotproductmatrix)


# In[133]:

PlMat89 = eval_legendre(splitell[88][:, None, None], dotproductmatrix)


# In[134]:

PlMat90 = eval_legendre(splitell[89][:, None, None], dotproductmatrix)


# In[135]:

PlMat91 = eval_legendre(splitell[90][:, None, None], dotproductmatrix)


# In[136]:

PlMat92 = eval_legendre(splitell[91][:, None, None], dotproductmatrix)


# In[137]:

PlMat93 = eval_legendre(splitell[92][:, None, None], dotproductmatrix)


# In[138]:

PlMat94 = eval_legendre(splitell[93][:, None, None], dotproductmatrix)


# In[139]:

PlMat95 = eval_legendre(splitell[94][:, None, None], dotproductmatrix)


# In[140]:

PlMat96 = eval_legendre(splitell[95][:, None, None], dotproductmatrix)


# In[141]:

PlMat97 = eval_legendre(splitell[96][:, None, None], dotproductmatrix)


# In[142]:

PlMat98 = eval_legendre(splitell[97][:, None, None], dotproductmatrix)


# In[143]:

PlMat99 = eval_legendre(splitell[98][:, None, None], dotproductmatrix)


# In[144]:

PlMat100 = eval_legendre(splitell[99][:, None, None], dotproductmatrix)


# In[145]:

PlMat101 = eval_legendre(splitell[100][:, None, None], dotproductmatrix)


# In[146]:

PlMat102 = eval_legendre(splitell[101][:, None, None], dotproductmatrix)


# In[147]:

PlMat103 = eval_legendre(splitell[102][:, None, None], dotproductmatrix)


# In[148]:

PlMat104 = eval_legendre(splitell[103][:, None, None], dotproductmatrix)


# In[149]:

PlMat105 = eval_legendre(splitell[104][:, None, None], dotproductmatrix)


# In[150]:

PlMat106 = eval_legendre(splitell[105][:, None, None], dotproductmatrix)


# In[151]:

PlMat107 = eval_legendre(splitell[106][:, None, None], dotproductmatrix)


# In[152]:

PlMat108 = eval_legendre(splitell[107][:, None, None], dotproductmatrix)


# In[153]:

PlMat109 = eval_legendre(splitell[108][:, None, None], dotproductmatrix)


# In[154]:

PlMat110 = eval_legendre(splitell[109][:, None, None], dotproductmatrix)


# In[155]:

PlMat111 = eval_legendre(splitell[110][:, None, None], dotproductmatrix)


# In[156]:

PlMat112 = eval_legendre(splitell[111][:, None, None], dotproductmatrix)


# In[157]:

PlMat113 = eval_legendre(splitell[112][:, None, None], dotproductmatrix)


# In[158]:

PlMat114 = eval_legendre(splitell[113][:, None, None], dotproductmatrix)


# In[159]:

PlMat115 = eval_legendre(splitell[114][:, None, None], dotproductmatrix)


# In[160]:

PlMat116 = eval_legendre(splitell[115][:, None, None], dotproductmatrix)


# In[161]:

PlMat117 = eval_legendre(splitell[116][:, None, None], dotproductmatrix)


# In[162]:

PlMat118 = eval_legendre(splitell[117][:, None, None], dotproductmatrix)


# In[163]:

PlMat119 = eval_legendre(splitell[118][:, None, None], dotproductmatrix)


# In[164]:

PlMat120 = eval_legendre(splitell[119][:, None, None], dotproductmatrix)


# In[165]:

PlMat121 = eval_legendre(splitell[120][:, None, None], dotproductmatrix)


# In[166]:

PlMat122 = eval_legendre(splitell[121][:, None, None], dotproductmatrix)


# In[167]:

PlMat123 = eval_legendre(splitell[122][:, None, None], dotproductmatrix)


# In[168]:

PlMat124 = eval_legendre(splitell[123][:, None, None], dotproductmatrix)


# In[169]:

PlMat125 = eval_legendre(splitell[124][:, None, None], dotproductmatrix)


# In[170]:

PlMat126 = eval_legendre(splitell[125][:, None, None], dotproductmatrix)


# In[171]:

PlMat127 = eval_legendre(splitell[126][:, None, None], dotproductmatrix)


# In[172]:

PlMat128 = eval_legendre(splitell[127][:, None, None], dotproductmatrix)


# In[173]:

PlMat129 = eval_legendre(splitell[128][:, None, None], dotproductmatrix)


# In[174]:

PlMat130 = eval_legendre(splitell[129][:, None, None], dotproductmatrix)


# In[175]:

PlMat131 = eval_legendre(splitell[130][:, None, None], dotproductmatrix)


# In[176]:

PlMat132 = eval_legendre(splitell[131][:, None, None], dotproductmatrix)


# In[177]:

PlMat133 = eval_legendre(splitell[132][:, None, None], dotproductmatrix)


# In[178]:

PlMat134 = eval_legendre(splitell[133][:, None, None], dotproductmatrix)


# In[179]:

PlMat135 = eval_legendre(splitell[134][:, None, None], dotproductmatrix)


# In[180]:

PlMat136 = eval_legendre(splitell[135][:, None, None], dotproductmatrix)


# In[181]:

PlMat137 = eval_legendre(splitell[136][:, None, None], dotproductmatrix)


# In[182]:

PlMat138 = eval_legendre(splitell[137][:, None, None], dotproductmatrix)


# In[183]:

PlMat139 = eval_legendre(splitell[138][:, None, None], dotproductmatrix)


# In[184]:

PlMat140 = eval_legendre(splitell[139][:, None, None], dotproductmatrix)


# In[185]:

PlMat141 = eval_legendre(splitell[140][:, None, None], dotproductmatrix)


# In[186]:

PlMat142 = eval_legendre(splitell[141][:, None, None], dotproductmatrix)


# In[187]:

PlMat143 = eval_legendre(splitell[142][:, None, None], dotproductmatrix)


# In[188]:

PlMat144 = eval_legendre(splitell[143][:, None, None], dotproductmatrix)


# In[189]:

PlMat145 = eval_legendre(splitell[144][:, None, None], dotproductmatrix)


# In[190]:

PlMat146 = eval_legendre(splitell[145][:, None, None], dotproductmatrix)


# In[191]:

PlMat147 = eval_legendre(splitell[146][:, None, None], dotproductmatrix)


# In[192]:

PlMat148 = eval_legendre(splitell[147][:, None, None], dotproductmatrix)


# In[193]:

PlMat149 = eval_legendre(splitell[148][:, None, None], dotproductmatrix)


# In[194]:

PlMat150 = eval_legendre(splitell[149][:, None, None], dotproductmatrix)


# In[195]:

splitell[49]


# In[196]:


PlMat50 = np.concatenate((PlMat1, PlMat2, PlMat3, PlMat4, PlMat5, PlMat6, PlMat7,
                                 PlMat8, PlMat9, PlMat10, PlMat11, PlMat12, PlMat13, PlMat14, PlMat15, 
                                PlMat16, PlMat17, PlMat18, PlMat19, PlMat20, PlMat21, PlMat22, PlMat23, 
                                PlMat24, PlMat25, PlMat26, PlMat27, PlMat28, PlMat29, PlMat30, PlMat31, 
                                PlMat32, PlMat33, PlMat34, PlMat35, PlMat36, PlMat37, PlMat38, PlMat39, 
                                 PlMat40, PlMat41, PlMat42, PlMat43, PlMat44, PlMat45, PlMat46, PlMat47,
                                 PlMat48, PlMat49, PlMat50))


# In[197]:

print(PlMat50.shape)





PlMat100 = np.concatenate((PlMat51, PlMat52, PlMat53, PlMat54, PlMat55,
                           PlMat56, PlMat57, PlMat58, PlMat59, PlMat60, PlMat61, PlMat62, PlMat63,
                           PlMat64, PlMat65, PlMat66, PlMat67, PlMat68, PlMat69, PlMat70, PlMat71,
                           PlMat72, PlMat73, PlMat74, PlMat75, PlMat76, PlMat77, PlMat78, PlMat79,
                           PlMat80, PlMat81, PlMat82, PlMat83, PlMat84, PlMat85, PlMat86, PlMat87,
                           PlMat88, PlMat89, PlMat90, PlMat91, PlMat92, PlMat93, PlMat94, PlMat95,
                           PlMat96, PlMat97, PlMat98, PlMat99, PlMat100))

print(PlMat100.shape)



PlMat150 = np.concatenate((PlMat101, PlMat102, PlMat103, PlMat104, PlMat105, PlMat106, PlMat107, PlMat108,
                           PlMat109, PlMat110, PlMat111, PlMat112, PlMat113, PlMat114, PlMat115, PlMat116,
                           PlMat117, PlMat118, PlMat119,
                           PlMat120, PlMat121, PlMat122, PlMat123, PlMat124, PlMat125, PlMat126, PlMat127,
                           PlMat128, PlMat129, PlMat130, PlMat131, PlMat132, PlMat133, PlMat134, PlMat135,
                           PlMat136, PlMat137, PlMat138, PlMat139, PlMat140, PlMat141, PlMat142, PlMat143,
                           PlMat144, PlMat145, PlMat146, PlMat147, PlMat148, PlMat149, PlMat150))

print(PlMat150.shape)



#
# save this result and upload in another IPython notebook
#



PlM_50 = "cl_varyCDMlmax1600ptPlMat50fakesigma25.npy"
PlM_100 = "cl_varyCDMlmax1600ptPlMat100fakesigma25.npy"
PlM_150 = "cl_varyCDMlmax1600ptPlMat150fakesigma25.npy"


np.save(PlM_50,  PlMat50)
np.save(PlM_100, PlMat100)
np.save(PlM_150, PlMat150)



print(PlMat50.shape)
print(PlMat100.shape)
print(PlMat150.shape)







