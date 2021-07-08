import numpy as np
import numpy.ma as ma
from scipy.spatial.distance import cdist
import pandas as pd


def compute_dist(A,B):
	M,dim=A.shape
	N=B.shape[0]	
	dist=[]
	for K in range(M):
		C=np.ones((N,1)) * A[K,]
		D=(C-B)*(C-B)
		dist.append(D.min())	
	return max(dist)

# Quantile Hausdorff
def quantile_hausdorff(A,B,QUANTILE=0.75):
	"""
	dH(A, B) = max(h(A, B),h(B, A)),  
	where h(A, B) = max(min(d(a, b))),  
	and d(a, b) is a L2 norm. 
	dist_H = hausdorff( A, B ) 
	A: First point sets. 
	B: Second point sets. 
	************************************************
	** A and B may have different number of rows, **
	** but must have the same number of columns   **
	************************************************
	Hassan RADVAR-ESFAHLAN; Universit du Qubec; TS; Montral; CANADA 
	15.06.2010
	"""

	# remove nans from the inputs
	A=ma.compress_rows(ma.array(A,mask=np.isnan(A)))
	B=ma.compress_rows(ma.array(B,mask=np.isnan(B)))
		
	A[:,0]=ma.masked_inside(A[:,0],B[:,0].min()-1,B[:,0].max()+1)
	A[:,1]=ma.masked_inside(A[:,1],B[:,1].min()-1,B[:,1].max()+1)
	
	if A.shape[1] != B.shape[1]:
		print('WARNING: dimensionality must be the same')
		dist=[]
		return dist
	
	D=cdist(A,B)
	D=pd.Series(D.min(axis=1))
	
	return D.quantile(q=QUANTILE,interpolation='linear')

# Classic Hausdorff metric
def hausdorff(A,B):
	"""
	dH(A, B) = max(h(A, B),h(B, A)),  
	where h(A, B) = max(min(d(a, b))),  
	and d(a, b) is a L2 norm. 
	dist_H = hausdorff( A, B ) 
	A: First point sets. 
	B: Second point sets. 
	************************************************
	** A and B may have different number of rows, **
	** but must have the same number of columns   **
	************************************************
	Hassan RADVAR-ESFAHLAN; Universit du Qubec; TS; Montral; CANADA 
	15.06.2010
	"""

	# remove nans from the inputs
	A=ma.compress_rows(ma.array(A,mask=np.isnan(A)))
	B=ma.compress_rows(ma.array(B,mask=np.isnan(B)))
		
	A[:,0]=ma.masked_inside(A[:,0],B[:,0].min()-1,B[:,0].max()+1)
	A[:,1]=ma.masked_inside(A[:,1],B[:,1].min()-1,B[:,1].max()+1)
	
	if A.shape[1] != B.shape[1]:
		print('WARNING: dimensionality must be the same')
		dist=[]
		return dist
			
	return max(compute_dist(A,B),compute_dist(B,A))

# Modified Hausdorff metric
def mod_hausdorff(A,B):
	"""
	This function computes the Modified Hausdorff Distance (MHD) which is 
	proven to function better than the directed HD as per Dubuisson et al. 
	in the following work:
		
		M. P. Dubuisson and A. K. Jain. A Modified Hausdorff distance for object 
		matching. In ICPR94, pages A:566-568, Jerusalem, Israel, 1994.
		http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=576361
		
	The function computed the forward and reverse distances and outputs the 
	minimum of both.

	Format for calling function:

	MHD = mod_hausdorff(A,B)

	where
	MHD = Modified Hausdorff Distance.
	A -> Point set 1
	B -> Point set 2
	
	No. of samples of each point set may be different but the dimension of
	the points must be the same.
	"""

	# remove nans from the inputs
	A=ma.compress_rows(ma.array(A,mask=np.isnan(A)))
	B=ma.compress_rows(ma.array(B,mask=np.isnan(B)))
	
	A[:,0]=ma.masked_inside(A[:,0],B[:,0].min()-1,B[:,0].max()+1)
	A[:,1]=ma.masked_inside(A[:,1],B[:,1].min()-1,B[:,1].max()+1)
		
	# Check if the points have the same dimensions
	if A.shape[1] != B.shape[1]:
		error('The dimensions of points in the two sets are not equal');

	# Calculating the forward HD

	D = cdist(A,B)
	
	# Calculating the forward HD
	fhd = D.min(axis=1).sum() / A.shape[0]

	# Calculating the reverse HD
	rhd = D.min(axis=0).sum() / B.shape[0]

	# Find the minimum of fhd/rhd as the mod hausdorff dist
	return max(fhd,rhd)

