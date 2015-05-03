
# coding: utf-8

# In[136]:

get_ipython().magic(u'matplotlib inline')

# We use CAMB generated scalar C_l values, N_side=4
# and CAMB simulated maps associated with scalar C_l values above
# N_side=4 gives N_pix=192


# In[137]:

import math
import matplotlib.pyplot as plt 
import numpy as np
import healpy as hp
import pyfits as pf
import astropy as ap
import os

np.set_printoptions(threshold='nan')  # Default is threshold=1000
# Print all is threshold = 'nan'
## Use this to print all values, disables corner printing


# In[138]:

cd ~/downloads


# In[139]:

file = "camb_cls_nside4.fits" # CAMB C_l scalars
temp = "camb_nside4.fits" # CAMB simulated maps


# In[140]:

# open a FITS file 
# type()=pyfits.hdu.hdulist.HDUList
ff = pf.open(file)


# In[141]:

### Recall there are four columns: temp, E pol, B pol, grad-temp cross terms
## first two values are zero, i.e. monopole, dipole
cls = ff[1].data  # actually second HDU, first is empty


# In[142]:

# XXX.field() references columns by 0-index
# field(0) is temperature values
# all Cl scalar temp values put into ndarray
# type()=numpy.ndarray
# len(cltemp2) = 385 (lmax=2Npix)
cltemp = cls.field(0) 


# In[143]:

# define ell values
# array from 0 to lmax, the size of map
ll = np.arange(385)


# In[144]:

# P_0 is the monopole, P_1 is the dipole
# remove 0, 1

ll_nodipoles = np.delete(ll, [0,1]) #numpy.ndarray, [2,3,..,384,385], len()=385 


# In[145]:

# First calculate the covariance matrix by the definition, i.e. 
# 
# C_ij = < \delta T_i \delta T_j > = 1/4pi \sum^{N_pix} p = 
# ( T_i(p) - mean(T_i) ) * ( T_j(p) - mean(T_j) )
#
#
# Larson, Weiland, Hinshaw, Bennett, 
# http://arxiv.org/pdf/1409.7718.pdf
#


# In[146]:

tempmap = hp.mrdfits("camb_nside4.fits") #type() = list


# In[147]:

tempdata = tempmap[0] #len()=192, type()=numpy.ndarray


# In[148]:

# First calculate the mean
#
# numpy.mean
#
# numpy.mean(a, axis=None, dtype=None, out=None, keepdims=False)
# Compute the arithmetic mean along the specified axis.


# Use either mean1 = np.mean(tempdata1) or mean1 = np.average(tempdata1)


# In[149]:

mean = np.average(tempdata) 


# In[150]:

Tpi = (tempdata - mean) #type()=numpy.ndarray


# In[151]:

Tpimatrix = np.matrix(Tpi) #Tpimatrix.shape = (1, 192)


# In[152]:

transpose = Tpimatrix.T #shape = (192, 1)


# In[153]:

print "This is the temperature pixel-pixel covariance matrix, equation (1)"
print "C_{ij} = <\Delta T_i \Delta T_j >"
print "Matrix shape 192 by 192"
cov = (1/(4*math.pi)) * ( transpose * Tpimatrix ) #covariance matrix, shape (192, 192)

print cov

print "**************"
print "**************"


# In[154]:

### Begin calculating S_ij piece by piece, in order to do the summation correctly 
#
# S_ij = sum(2ell+1) C_l P_l(dotproductmatrix)


# In[155]:

ell = np.arange(192)


# In[156]:

## will give three arrays
## arrays of all x values, all y values, all z values
## RING scheme default
# len()=3
# type()=tuple
vecval = hp.pix2vec(4, ell) #Nside = 4, type()=tuple


# In[157]:

vecvalx = vecval[0] #shape (192,)
vecvaly = vecval[1]
vecvalz = vecval[2]



# In[158]:

## First arrange arrays vertically
## numpy.vstack = Stack arrays in sequence vertically (row wise), input sequence of arrays

totalvecval = np.vstack((vecvalx, vecvaly, vecvalz)) #type()=numpy.ndarray


# In[159]:

trans = totalvecval.T #transpose


# In[160]:

dotproductmatrix = trans.dot(totalvecval) #take the dot product


# In[161]:

## We set C_l scalar TT values to lmax=2Npix
## Fix this by simply deleting high C_l values
# print np.arange(193,384)

## our indices to remove, arr

arr = np.arange(193,384) #numpy.ndarray such that [193,,...383, 384]


# In[162]:

# We screwed up above. Write lmax = Npix. Remove excess C_l values. 

# [x for i,x in enumerate(a) if i not in ind2remove]

# enumerate() is a build-in python function

# clscalar = np.delete(cltemp, [192:384])


newcls = [ x for i,x in enumerate(cltemp) if i not in arr]


# In[163]:

from scipy.special import eval_legendre  ##special scipy function


# In[164]:

## Begin calculating S_ij piece by piece, in order to do the summation correctly. 
#
# S_ij = sum(2ell+1) C_l P_l(dotproductmatrix)

# NOT QUICK!

summatrix = sum( [eval_legendre(i, dotproductmatrix) for i in ell])




# In[165]:

# matrix_total = 
# (1/(4*math.pi)) * sum((2 * ll + 1) * cltemp ) * eval_legendre(ll, matrix_dotprod)
#
# Begin with adding theoretical scalar C_l values
#
add_clvalues = sum([ i * summatrix for i in newcls ])


# In[166]:

# matrix_total = 
# (1/(4*math.pi)) * np.sum((2*ll + 1) * cltemp ) * eval_legendre(ll, matrix_dotprod)

wholematrix = sum([((2 * i) + 1) * add_clvalues for i in newcls ])



# In[167]:


covmatrix = (1/(4 * math.pi)) * wholematrix #covariance matrix for Nside=4



# In[168]:

print "covmatrix matrix S_ij"
print "covmatrix.shape = (192, 192), covmatrix2.size = 36864"

print covmatrix

print "**************"
print "**************"


# In[169]:

### We now wish to calculate SUM(2ell+1)C_l
### We do this two ways: we first calculate using scalar C_l values
### We then take map data, use healpy.sphtfunc.anafast
### healpy.sphtfunc.anafast(map1, map2=None, nspec=None, lmax=None, 
### mmax=None, iter=3, alm=False, pol=True, use_weights=False, datapath=None)
### returns cl or a list of cl’s (TT, EE, BB, TE, EB, TB for polarized input map) 
### Here, only TT cl
ell[0:4]


# In[170]:

# len(newcls) = 194
# len(ell) = 192
# ell[0:4] = array([0, 1, 2, 3])
# newcls[0:4] = [0.0, 0.0, 1.2639207e-09, 5.8827504e-10]
# monopole, dipole---both zero values
## Therefore, reshape by removing monopole term
# ell[0:4] = array([0, 1, 2, 3])

newnewcls = np.delete(newcls, [0,1])
# len(newnewcls) = 192
newnewcls[0:4]
# newnewcls[0:4] = array([  1.26392075e-09,   5.88275040e-10,   3.28673144e-10,
#                          2.06575299e-10], dtype=float32)


# In[171]:

clsvalmatrix = np.matrix(newnewcls) # clsvalmatrix.shape = (1, 192)
ellmatrix = np.matrix(ell) # ellmatrix.shape = (1,192)
# ellmatrix.T.shape = (192, 1)


# In[172]:

## Begin with C_l scalar values and simply sum
# spherharmcov = sum([((2 * i) + 1) * newcls for i in ell ])
#
### for i in ellmatrix.T:
###    spherharmcov = ((2 * i) + 1) * clsvalmatrix


# In[173]:

for i in clsvalmatrix: #clsvalmatrix.shape = (1, 192)
    harmcov = ((2*ellmatrix.T + 1)) * i #ellmatrix.T.shape = (192, 1)


# In[174]:

print "Covariance matrix from generate scalar C_l values"
print "Summing (2l+1)C_l, beginning with l=0"
print "Unclear to me whether should use l=0,1,2 or l=1,2,3"

print harmcov

print "**************"
print "**************"


# In[175]:

## Now using CAMB temperature map
# Generate C_l values from map with anafast
# Repeat above procedure
### UNDER CONSTRUCTION


# In[176]:

ell_193 = np.arange(193)


# In[177]:

newells = np.delete(ell_193, [0])
# print newells = [ 1 2 ..... 191 192]
newellmatrix = np.matrix(newells)


# In[178]:

for i in clsvalmatrix: #clsvalmatrix.shape = (1, 192)
    harmcovLone = ((2*newellmatrix.T + 1)) * i #ellmatrix.T.shape = (192, 1)


# In[179]:

print "Covariance matrix from generate scalar C_l values"
print "Summing (2l+1)C_l, beginning with l=1"

print harmcovLone

print "**************"
print "**************"


# In[180]:

# We now compute likelihood functions
#
# For equations, refer to Nside4_covariance_likelihood.tex
# "Nside = 4: Covariance and likelihood"
#


# In[181]:

# map generated from scalar C_l 
map = hp.synfast(newcls, 4)
# len(map) = 192
### healpy.sphtfunc.synfast(cls, nside, lmax=None, mmax=None, alm=False, pol=True, 
###                         pixwin=False, fwhm=0.0, sigma=None, new=False, verbose=True)
# INPUT cls, nside
# OUTPUT sky maps


# In[182]:

powerspectrum = hp.anafast(map) # create powerspectrum \hat{C}_l
# len(powerspectrum) = 12


# In[183]:

print "HEALPix anafast \hat(C)_l estimator, len()=12"

print powerspectrum

print "Difficult to see how to compute spherical harmonic likelihood with CAMB generate maps"
print "Estimator from HEALPix is len()=12, theoretical C_l is 192"
print "Test with simulated Planck data GOVA"

print "**************"
print "**************"


# In[184]:

# Compute temperature space (real space) likelihood
#
# where T is a temperature map array (3072,), S is our above matrix = entirething, 
# and N is the number of pixels in a vector
#
# First Year Wilkinson Microwave Anisotropy Probe (WMAP) Observations: 
# Parameter Estimation Methodology L. Verde, et al., 2003, ApJS, 148, 195
#
# http://lambda.gsfc.nasa.gov/product/map/dr1/
# pub_papers/firstyear/methodology/wmap_param_method.pdf
#
# print "−2 ln L \propto T*S^{-1}*T + ln detS + N ln 2\pi"


# In[185]:

tempvalues = np.matrix(tempdata) #create matrix of temperature values from CAMB map
# tempvalues.shape = (1, 192)


# In[186]:

temptranpose = tempvalues.T
# temptranpose.shape = (192, 1)


# In[187]:

## Invert covariance matrix, S_ij
## numpy.linalg.inv(a)
## Compute the (multiplicative) inverse of a matrix.
## Given a square matrix a, return the matrix ainv satisfying 
## dot(a, ainv) = dot(ainv, a) = eye(a.shape[0]).
##

inverse = np.linalg.inv(covmatrix)
inversematrix = np.matrix(inverse) # S^{-1}
# inverse.shape = (192, 192)


# In[188]:

# WARNING: DO NOT USE NP.LIGALG.DET(A)! 
# Operatioin will not succeed, gives inf or zero
#
#
#
# numpy.linalg.slogdet
# numpy.linalg.slogdet(a)[source]
# Compute the sign and (natural) logarithm of the determinant of an array.
# INPUT array
# OUTPUT (sign, ln det a)
# sign = A number representing the sign of the determinant. 
# For a real matrix, this is 1, 0, or -1. 
# For a complex matrix, this is a complex number with absolute value 1 (i.e., it is on the unit circle), or else 0.
# ln det a = The natural log of the absolute value of the determinant.
#
lnDetS = np.linalg.slogdet(inversematrix) #type()=tuple
#print lnDetS = (1.0, 2203.8440835363299)
#type(lnDetS[1]) = numpy.float64


# In[189]:

type(temptranpose)


# In[190]:

realtemplikehd = ( tempvalues * inversematrix * temptranpose ) 
realtemplikehdpt2 = ( tempvalues * inversematrix * temptranpose ) + lnDetS[1]
realtemplikehdpt3 = ( tempvalues * inversematrix * temptranpose ) + ( lnDetS[1] ) + ((192)*2*math.pi)


# In[191]:

print "Real space (temperature space) likelihood function"
print "Formalism:  −2 ln L \propto T*S^{-1}*T + ln detS + N ln 2\pi"
print "First number: T*S^(-1)*T "

print realtemplikehd

print "Second number: T*S^{-1}*T + ln detS "

print realtemplikehdpt2

print "Third number:  −2 ln L \propto T*S^{-1}*T + ln detS + N ln 2\pi"

print realtemplikehdpt3

print "**************"
print "**************"
print "END"


# In[ ]:




# In[ ]:



