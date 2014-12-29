#!/bin/python2
from numpy import *;

#def css(M,k):

def initial_stage(M, k):
    r,c    = M.shape;
    u,d,vt = linalg.svd(M); 
    p      = array([0 for i in xrange(c)]);
    for j in xrange(c):
        for i in xrange(k):
            p[j] += vt[i,j]*vt[i,j];
        p[j] /= k;
    
   
#def randomized_stage():

#def deterministic_stage(M):

