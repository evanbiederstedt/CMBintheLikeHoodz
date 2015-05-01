
# coding: utf-8

# In[41]:

get_ipython().magic(u'matplotlib inline')

# We use CAMB generated scalar C_l values, set at N_side=2, N_side=4, N_side=8. 
# That is total pixels Npix=48, Npix=192, Npix=768
# This amounts to lmax=96, lmax=384, lmax=1536


# In[42]:

import math
import matplotlib.pyplot as plt 
import numpy as np
import healpy as hp
import pyfits as pf
import astropy as ap
import os


# In[43]:

cd ~/downloads


# In[44]:

file1 = "camb_cls_nside2.fits"
file2 = "camb_cls_nside4.fits"
file3 = "camb_cls_nside8.fits"


# In[45]:

ff1 = pf.open(file1) # open a FITS file
ff2 = pf.open(file2) 
ff3 = pf.open(file3) #type()=pyfits.hdu.hdulist.HDUList


# In[46]:

###Recall there are four columns: temp, E pol, B pol, grad-temp cross terms
##first two values are zero, i.e. monopole, dipole
cls1 = ff1[1].data  # actually second HDU, first is empty
cls2 = ff2[1].data 
cls3 = ff3[1].data 


# In[47]:

#XXX.field() references columns by 0-index
##field(0) is temperature values

cltemp1 = cls1.field(0) #all Cl scalar temp values put into ndarray
cltemp2 = cls2.field(0)
cltemp3 = cls3.field(0) #type()=numpy.ndarray

#len(cltemp1) = 97
#len(cltemp2) = 385
#len(cltemp3) = 1537


# In[48]:

##define ell values
#array from 0 to lmax, the size of map

ll_1 = np.arange(97) 
ll_2 = np.arange(385)
ll_3 = np.arange(1537)


# In[49]:

##P_0 is the monopole, P_1 is the dipole
##remove 0, 1

ll_nodipoles1 = np.delete(ll_1, [0,1]) #numpy.ndarray, [2,3,..,95,96], len()=95
ll_nodipoles2 = np.delete(ll_2, [0,1])
ll_nodipoles3 = np.delete(ll_3, [0,1])


# In[50]:

# First calculate the covariance matrix by the definition, i.e. 
# 
# C_ij = < \delta T_i \delta T_j > = 1/4pi \sum^{N_pix} p = 
# ( T_i(p) - mean(T_i) ) * ( T_j(p) - mean(T_j) )
#
#
# Larson, Weiland, Hinshaw, Bennett, 
# http://arxiv.org/pdf/1409.7718.pdf
#


# In[51]:

temp1 = "camb_nside2.fits" #CAMB simulated maps associated with scalar C_l values above
temp2 = "camb_nside4.fits"
temp3 = "camb_nside8.fits"


# In[52]:

tempmap1 = hp.mrdfits("camb_nside2.fits") #type() = list
tempmap2 = hp.mrdfits("camb_nside4.fits")
tempmap3 = hp.mrdfits("camb_nside8.fits")


# In[53]:

tempdata1 = tempmap1[0] #len()=48, type()=numpy.ndarray
tempdata2 = tempmap2[0] #len()=192
tempdata3 = tempmap3[0] #len()=768


# In[54]:

# First calculate the mean
#
# numpy.mean
#
# numpy.mean(a, axis=None, dtype=None, out=None, keepdims=False)
# Compute the arithmetic mean along the specified axis.


# Use either mean1 = np.mean(tempdata1) or mean1 = np.average(tempdata1)


# In[55]:

mean1 = np.average(tempdata1) #print mean1 = -3.92499e-07
mean2 = np.average(tempdata2) #print mean2 = -1.006e-07
mean3 = np.average(tempdata3) #print mean3 = -1.76321e-08


# In[56]:

Tpi1 = (tempdata1 - mean1) #type()=numpy.ndarray
Tpi2 = (tempdata2 - mean2)
Tpi3 = (tempdata3 - mean3)


# In[57]:

Tp1 = np.matrix(Tpi1) #Tp1.shape = (1, 48)
Tp2 = np.matrix(Tpi2) #Tp2.shape = (1, 192)
Tp3 = np.matrix(Tpi3) #Tp3.shape = (1, 768)


# In[58]:

transpose1 = Tp1.T #shape = (48, 1)
transpose2 = Tp2.T #shape = (192, 1)
transpose3 = Tp3.T #shape = (768, 1)


# In[59]:

cov1 = (1/(4*math.pi)) * ( transpose1 * Tp1 ) #covariance matrix 1, shape (48, 48)
cov2 = (1/(4*math.pi)) * ( transpose2 * Tp2 ) #covariance matrix 2, shape (192, 192)
cov3 = (1/(4*math.pi)) * ( transpose3 * Tp3 ) #covariance matrix 3, shape (768, 768)


# In[60]:

###Begin calculating S_ij piece by piece, in order to do the summation correctly. 
#
# S_ij = sum(2ell+1) C_l P_l(dotproductmatrix)


# In[61]:

ell_1 = np.arange(48)
ell_2 = np.arange(192)
ell_3 = np.arange(768)


# In[62]:

vecval_1 = hp.pix2vec(2, ell_1) #Nside = 2, type()=tuple
vecval_2 = hp.pix2vec(4, ell_2) #Nside = 4
vecval_3 = hp.pix2vec(8, ell_3) #Nside = 8
###will give three arrays
##arrays of all x values, all y values, all z values
##RING scheme default

#len()=3
#type()=tuple


# In[63]:

vecvalx1 = vecval_1[0] #shape (48,)
vecvaly1 = vecval_1[1]
vecvalz1 = vecval_1[2]

vecvalx2 = vecval_2[0] #shape (192,)
vecvaly2 = vecval_2[1]
vecvalz2 = vecval_2[2]

vecvalx3 = vecval_3[0] #shape (768,) 
vecvaly3 = vecval_3[1]
vecvalz3 = vecval_3[2]


# In[64]:

###First arrange arrays vertically
##numpy.vstack = Stack arrays in sequence vertically (row wise), input sequence of arrays

totalvecval_1 = np.vstack((vecvalx1, vecvaly1, vecvalz1)) #type()=numpy.ndarray
totalvecval_2 = np.vstack((vecvalx2, vecvaly2, vecvalz2))
totalvecval_3 = np.vstack((vecvalx3, vecvaly3, vecvalz3))


#type(totalvecval) = numpy.ndarray


# In[65]:

trans1 = totalvecval_1.T
trans2 = totalvecval_2.T
trans3 = totalvecval_3.T


# In[66]:

dotproductmatrix1 = trans1.dot(totalvecval_1)
dotproductmatrix2 = trans2.dot(totalvecval_2)
dotproductmatrix3 = trans3.dot(totalvecval_3)


# In[67]:

###We set C_l scalar TT values to lmax=2Npix
###Fix this by simply deleting high C_l values

#print np.arange(48,96)

#len(cltemp1) = 97


###our indices to remove, arr1, arr2, arr3

arr1 = np.arange(48,97,1) #numpy.ndarray such that [48,49,...95,96]
arr2 = np.arange(192,384)
arr3 = np.arange(768,1536)


# In[68]:

#We screwed up above. Write lmax = Npix. Remove excess C_l values. 

 #[x for i,x in enumerate(a) if i not in ind2remove]

#enumerate() is a build-in python function

#clscalar1 = np.delete(cltemp1, [48:96])
#clscalar2 = np.delete(cltemp2, [192:384])
#clscalar3 = np.delete(cltemp3, [768:1536])

newcls1 = [ x for i,x in enumerate(cltemp1) if i not in arr1]
newcls2 = [ x for i,x in enumerate(cltemp2) if i not in arr2]
newcls3 = [ x for i,x in enumerate(cltemp3) if i not in arr3]


# In[69]:

from scipy.special import eval_legendre


# In[70]:

###Begin calculating S_ij piece by piece, in order to do the summation correctly. 
#
# S_ij = sum(2ell+1) C_l P_l(dotproductmatrix)

summatrix1 = sum( [eval_legendre(i, dotproductmatrix1) for i in ell_1])
summatrix2 = sum( [eval_legendre(i, dotproductmatrix2) for i in ell_2])

#summatrix3 = sum( [eval_legendre(i, dotproductmatrix3) for i in ell_3])



# In[71]:

#matrix_total = 
#(1/(4*math.pi)) * sum((2 * ll + 1) * cltemp ) * eval_legendre(ll, matrix_dotprod)
#
# Begin with adding theoretical scalar C_l values
#
add_clvalues1 = sum([ i * summatrix1 for i in newcls1 ])
add_clvalues2 = sum([ i * summatrix2 for i in newcls2 ])

#add_clvalues3 = sum([ i * summatrix3 for i in newcls3 ])


# In[72]:

#matrix_total = 
#(1/(4*math.pi)) * np.sum((2*ll + 1) * cltemp ) * eval_legendre(ll, matrix_dotprod)

wholematrix1 = sum([((2 * i) + 1) * add_clvalues1 for i in newcls1 ])
wholematrix2 = sum([((2 * i) + 1) * add_clvalues2 for i in newcls2 ])

#wholematrix3 = sum([((2 * i) + 1) * add_clvalues3 for i in newcls3 ])


# In[73]:

covmatrix1 = (1/(4 * math.pi)) * wholematrix1 #covariance matrix for Nside=2
covmatrix2 = (1/(4 * math.pi)) * wholematrix2 #covariance matrix for Nside=4

#covmatrix3 = (1/(4 * math.pi)) * wholematrix3 #covariance matrix for Nside=8



# In[74]:

np.set_printoptions(threshold='nan')  
##Use this to print all values, disables corner printing


# In[75]:

#Compare with what we calculated above

#covmatrix1.shape = (48, 48), covmatrix1.size = 2304
#covmatrix2.shape = (192, 192), covmatrix2.size = 36864
#covmatrix3.shape = (768, 768), covmatrix3.size = 589824


# In[76]:

###Comparison 1: 
###Temp-temp covariance matrix, cov1
###Dot product-generated covariance matrix, covmatrix 1
#

print cov1
print "******"
print covmatrix1


# In[77]:

###Comparison 2: 
###Temp-temp covariance matrix, cov2
###Dot product-generated covariance matrix, covmatrix 2
#

print cov2
print "******"
print covmatrix2


# In[78]:

###Comparison 3: 
###Temp-temp covariance matrix, cov3
###Dot product-generated covariance matrix, covmatrix 3
#

#print cov3
#print "******"
#print covmatrix3


# In[ ]:



