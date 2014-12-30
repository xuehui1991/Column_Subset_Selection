#!/bin/bash
# linear least square solution by householder transformation, 1965 by Golub

from numpy import *;
import math;
def sign(v):
    if v > 0:   return 1
    else:   return -1

M = array([(1,2,3),(3,4,6),(7,9,4)]);

def HouseHolder_step(M,k):
    (r,c) = M.shape    
    if r <= k or c <= k:    return M
    
    sigmak = 0;
    for i in xrange(k,r):
        sigmak += M[i,k] * M[i,k];
    sigmak = math.sqrt(sigmak);
    
    betak  = sigmak * (sigmak + abs(M[k,k])); 
    betak  = 1.0/betak;

    u = zeros([r,1]);
    u[k][0] = sign(M[k,k]) * (sigmak + abs(M[k,k]));
    for i in xrange(k+1,r):
        u[i, 0] = M[i,k];

    uut = dot(u, transpose(u));
    M = M - betak * dot(uut,M);
    
    return M;

        
def HouseHolder(M):
    (r,c) = M.shape;
    k = min(r,c);    

    for i in xrange(k):
        M = HouseHolder_step(M,i);
    return M;
