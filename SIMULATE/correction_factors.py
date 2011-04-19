# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from tables import *
from numpy import *
from pymc import invlogit, logit
from pymc.flib import logsum
from cf_helper import cfh, cfhs
# import magic_load
# magic_load.magic_load('age_pr_datasets','/home/anand/MAP/Age standardization/age_pr_datasets/__init__.py')
from age_pr_datasets import *
import auxiliary_data

from scipy import interpolate as interp

__all__ = ['two_ten_factors', 'age_corr_factors', 'known_age_corr_factors', 'known_age_corr_likelihoods_f', 'stochastic_known_age_corr_likelihoods', 'age_corr_likelihoods', 'age_corr_factors_from_limits']
proot = auxiliary_data.parameter_model
    
F_trace, P_trace = proot.col('F_pred'), proot.col('P_pred')

aroot = auxiliary_data.age_dist_model

S_trace = aroot.col('S_pred')

def two_ten_factors(N):
    """
    factors = two_ten_factors(N)
    
    Factors returns a length-N array of samples from the prior 
    distribution of r, where PfPR[2, 10) = r P'. P' is the 
    location-caused component.
    
    NB this function sets 'F' to 1; it assumes infections are 
    'detected' wp 1, or alternatively assumes we are interested 
    in computing the probability of being infected rather than 
    the probability an infection is detected.
    """
    factors = empty(N)
    p_indices = random.randint(P_trace.shape[0], size=N)
    S_indices = random.randint(S_trace.shape[0], size=N)
    
    P_samps = P_trace[p_indices]
    F_samps = F_trace[p_indices]
    S_samps = S_trace[S_indices,0,:]
                
    for j in xrange(N):
        P = P_samps[j,:]
        F = F_samps[j,:]
        S = S_samps[j,:]

        factors[j] = sum(S[2:10] * P[2:10]) / sum(S[2:10])
            
    return factors

def age_corr_factors_from_limits(a_lo, a_hi, N):
    p_indices = random.randint(P_trace.shape[0], size=N)
    S_indices = random.randint(S_trace.shape[0], size=N)
    
    P_samps = P_trace[p_indices]
    F_samps = F_trace[p_indices]
    S_samps = S_trace[S_indices,0,:]
    
    a_index_min = np.where(a<=a_lo)[0][-1]
    a_index_max = np.where(a>=a_hi)[0]
    if len(a_index_max)==0:
        a_index_max = len(a)-1
    else:
        a_index_max = a_index_max[0]
    
    factors = np.empty(N)
    for j in xrange(N):
        P = P_trace[p_indices[j],a_index_min:a_index_max+1]
        F = F_trace[p_indices[j],a_index_min:a_index_max+1]
        S = S_trace[S_indices[j],0,a_index_min:a_index_max+1]

        factors[j] = sum(S * P * F) / sum(S)
    return factors
    


def age_corr_factors(A, N):
    """
    factors = age_corr_factors(A,N)
    
    Factors returns a len(A) by N array of samples from the prior 
    distribution of r, where PfPR[a_min, a_max) = r P'. P' is the 
    location-caused component.
    """
    N_recs = len(A)
    factors = empty((N_recs, N))
    p_indices = random.randint(P_trace.shape[0], size=N)
    S_indices = random.randint(S_trace.shape[0], size=N)
    
    P_samps = P_trace[p_indices]
    F_samps = F_trace[p_indices]
    S_samps = S_trace[S_indices,0,:]
                
    for i in xrange(N_recs):
        
        a_index_min = np.where(a<=A.LOW_AGE[i])[0][-1]
        a_index_max = np.where(a>=A.UP_AGE[i])[0]
        if len(a_index_max)==0:
            a_index_max = len(a)-1
        else:
            a_index_max = a_index_max[0]
            
        for j in xrange(N):
            P = P_samps[j,a_index_min:a_index_max+1]
            F = F_samps[j,a_index_min:a_index_max+1]
            S = S_samps[j,a_index_min:a_index_max+1]
    
            factors[i,j] = sum(S * P * F) / sum(S)
            
    return factors

def known_age_corr_factors(A,N):
    factors = empty((len(A),N))
    p_indices = random.randint(P_trace.shape[0], size=N)
    S_indices = random.randint(S_trace.shape[0], size=N)
    
    P_samps = P_trace[p_indices]
    F_samps = F_trace[p_indices]
    # S_samps = S_trace[S_indices,0,:]
    
    for i in xrange(len(A)):
        factors[i,:] = P_samps[:,A[i]]*F_samps[:,A[i]]

    return factors

def outer_small(p1, fac_array, Ai, spi, posi, negi, likes_now):
    N = fac_array.shape[1]
    
    # For each sample from the posterior of r
    for j in xrange(N):
        # Compute actual PfPR[A[i]]
        p2 = fac_array[:,j][Ai]

        # Record likelihood
        # k1=np.add.outer(np.log(p1),np.log(p2))
        k1 = np.log(p1)*spi + np.dot(posi,np.log(p2))

        # k2 = log(1.-np.outer(p1,p2))
        k2 = cfh(p1,p2,negi)
        # k1 = np.dot(k1,posi)
        # k2 = np.dot(k2,negi)

        likes_now[j,:] = k1 + k2
        
    # Average log-likelihoods.
    return np.apply_along_axis(logsum,0,likes_now) - log(N)#(mean(likes_now,axis=0))

def outer_large(p1, fac_array, Ai, spi, posi, negi, likes_now):
    # TODO: Treat case of n sufficiently large as delta function or normal approx if n > npos > 0.
    # Remember, even if it's a delta function you have the distribution of the correction factors to deal with.
    # So spline up the correction factor histograms ahead of time, there will be thousands but it's OK.
    # If n = npos or n = 0, 
    N = fac_array.shape[1]
    likes_now *= 0
    
    Aset = set(Ai)
    pos = {}
    negs = {}
    sps = {}
    outs = {}
    for Au in Aset:
        where_this_age = where(Ai==Au)[0]
        pos = sum(posi[where_this_age])
        tot = len(where_this_age)
        neg = len(where_this_age) - pos
        
        p2 = fac_array[Au,:]
        k1 = np.add.outer(np.log(p2), np.log(p1))*pos
        k2 = log(1. - np.outer(p2, p1))*neg
        
        likes_now += k1 + k2
        
    # Average log-likelihoods.
    return np.apply_along_axis(logsum,0,likes_now) - log(N)#(mean(likes_now,axis=0))        
        

def known_age_corr_likelihoods_f(pos, A, fac_array, f_mesh, nug, type=None):
    """
    Computes spline representations over P_mesh for the likelihood 
    of N_pos | N_exam, A
    """

    # TODO: Optimize large-N case using CLT of some kind.

    # Allocate work and output arrays.
    N_recs = len(A)

    likelihoods = empty((N_recs, len(f_mesh)))
    likes_now = empty((fac_array.shape[1], len(f_mesh)), dtype=float128)
    splreps = []
    
    p1 = invlogit(f_mesh)
    
    # For each record
    for i in xrange(N_recs):
        posi = pos[i]
        Ai = A[i]
        spi = np.sum(posi)
        negi = 1.-posi

        if type is None:
            if len(Ai) < 100:
                fn = outer_small
            else:
                fn = outer_large
        elif type=='s':
            fn = outer_small
        else:
            fn = outer_large

        likelihoods[i,:] = fn(p1, fac_array, Ai, spi, posi, negi, likes_now)

        # Clean out occasional infinities on the edges.
        good_indices = where(1-isinf(likelihoods[i,:]))[0]

        # Compute spline representations.
        this_splrep = interp.splrep(x=f_mesh[good_indices], y=likelihoods[i,good_indices].squeeze())
        def this_fun(x, sp=this_splrep, Pml=f_mesh[good_indices].min(), Pmh=f_mesh[good_indices].max()):
            out = np.atleast_1d(interp.splev(x, sp))
            if np.any(x<Pml) or np.any(x>Pmh):
                out[np.where(x<Pml)] = -np.Inf
                out[np.where(x>Pmh)] = -np.Inf
            return out.reshape(np.shape(x))

        splreps.append(this_fun)        
    return splreps

def stochastic_known_age_corr_likelihoods(pos, A, fac_array):
    """
    Computes spline representations over P_mesh for the likelihood 
    of N_pos | N_exam, A
    """

    # Allocate work and output arrays.
    N_recs = len(A)
    N = fac_array.shape[1]

    funs = []

    import time

    # For each record
    for i in xrange(N_recs):
        posi = pos[i]
        Ai = A[i]
        spi = np.sum(posi)
        negi = 1.-posi

        j = np.random.randint(N)
        this_fac_array=fac_array[:,j]

        # Compute actual PfPR[A[i]]
        p2 = this_fac_array[Ai]
        
        p3 = np.dot(posi, np.log(p2))


        # Compute spline representations.
        def this_fun(x, p2=p2, p3=p3,negi=negi, posi=posi, Ai=Ai):
            p1 = np.log(invlogit(x))
            return p1*spi + p3 + cfh(p1,p2,negi)

        funs.append(this_fun)        
    return funs
    
# def known_age_corr_likelihoods(pos, A, fac_array, P_mesh):
#     """
#     Computes spline representations over P_mesh for the likelihood 
#     of N_pos | N_exam, A
#     """
#     
#     # Allocate work and output arrays.
#     N_recs = len(A)
#     N = fac_array.shape[1]
#     
#     likelihoods = empty((N_recs, len(P_mesh)))
#     likes_now = empty((N, len(P_mesh)), dtype=float128)
#     splreps = []
#     
#     # For each record
#     for i in xrange(N_recs):
#         Np = N_pos[i]
#         Ne = N_exam[i]
#         Ai = A[i]
# 
#         # For each sample from the posterior of r
#         for j in xrange(N):
#             # Compute actual PfPR[A[i]]
#             P_adjusted = fac_array[Ai,j] * P_mesh
#             # Record likelihood
#             likes_now[j,:] = Np * log(P_adjusted) + (Ne-Np) * log(1.-P_adjusted)
#             
#         # Average log-likelihoods.
#         likelihoods[i,:] = np.apply_along_axis(logsum,0,likes_now) - log(N)#(mean(likes_now,axis=0))
#         
#         # Clean out occasional infinities on the edges.
#         good_indices = where(1-isinf(likelihoods[i,:]))[0]
#         
#         # Compute spline representations.
#         this_splrep = interp.splrep(x=P_mesh[good_indices], y=likelihoods[i,good_indices].squeeze(), xb=0, xe=1)
#         def this_fun(x, sp=this_splrep, Pml=P_mesh[good_indices].min(), Pmh=P_mesh[good_indices].max()):
#             out = np.atleast_1d(interp.splev(x, sp))
#             out[np.where(x<Pml)] = -np.Inf
#             out[np.where(x>Pmh)] = -np.Inf
#             return out.reshape(np.shape(x))
#             
#         splreps.append(this_fun)
#         
#         # if np.any(this_splrep[1]>0):
#         
#     return splreps

def age_corr_likelihoods(A, N, P_mesh, datafile_name):
    """
    Returns samples from the log-likelihood p(N|A,P'=P_mesh,r), 
    where A is a MAP-style record array and r[i] is the factor 
    converting P' to PfPR[A[i]], and a spline representation of 
    the log marginal likelihood  p(N|A,P'=P_mesh) = 
    \int_r p(N|A,P'=P_mesh,r) p(r|A).
    """
    
    splreps_name = datafile_name + '_splreps.hdf5'
    if splreps_name in os.listdir('./'):
        splrep_file = openFile(splreps_name)
        splreps = splrep_file.root.splreps[:]
        likelihoods = splrep_file.root.likelihoods[:]
        splrep_file.close()
        return likelihoods, splreps
    
    print 'Recomputing likelihood spline representations...'
    # Call to age_corr_factors to get samples from the predictive distribution of r.
    factors = age_corr_factors(A,N)
    
    # Allocate work and output arrays.
    N_recs = len(A)
    likelihoods = empty((N_recs, len(P_mesh)))
    likes_now = empty((N, len(P_mesh)), dtype=float128)
    splreps = []
    indices = random.randint(factors.shape[1], size=N)

    splrep_file = openFile(splreps_name,'w')
    splrep_file.createVLArray(splrep_file.root,'splreps',ObjectAtom())
    splrep_file.createCArray(splrep_file.root, 'likelihoods', FloatAtom(), (N_recs, len(P_mesh)))

    # For each record
    for i in xrange(N_recs):
        N_pos = A.PF[i]
        N_exam = A.EXAMINED[i]

        # For each sample from the posterior of r
        for j in xrange(N):
            # Compute actual PfPR[A[i]]
            P_adjusted = factors[i,indices[j]] * P_mesh
            # Record likelihood
            likes_now[j,:] = N_pos * log(P_adjusted) + (N_exam - N_pos) * log(1.-P_adjusted)

        # Standardize, exponentiate and average log-likelihoods.
        likes_now = exp(likes_now - likes_now.max())
        likelihoods[i,:] = log(mean(likes_now,axis=0))
        
        # Clean out occasional infinities on the edges.
        good_indices = where(1-isinf(likelihoods[i,:]))
        
        # Compute spline representations.
        this_splrep = interp.splrep(x=P_mesh[good_indices], y=likelihoods[i,good_indices].squeeze(), xb=0, xe=1)
        splreps.append(this_splrep)
        splrep_file.root.splreps.append(this_splrep)

    print 'Done'
    print splrep_file.root.likelihoods.shape, (likelihoods - likelihoods.max()).shape
    splrep_file.root.likelihoods[:,:] = likelihoods - likelihoods.max()

    splrep_file.close()
    
    return likelihoods-likelihoods.max(), splreps
            
