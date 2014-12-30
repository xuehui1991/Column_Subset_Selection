#!/bin/python2
from numpy import *;
from CSS_Util import *;

#def css(M,k):

def initial_stage(A, k):
    row,col = A.shape;
    u,d,vt  = linalg.svd(A); 
    p       = array([0.0 for i in xrange(col)]);
    for j in xrange(col):
        for i in xrange(k):
            p[j] +=  vt[i,j]*vt[i,j];
        p[j] /= k;
    
    Af        = linalg.norm(A,'fro');
    c0_square =  4 * 0.8 * 0.8 / Af / Af;
    print c0_square;
    c         =  1600 *  c0_square  *  k * math.log( 800 * c0_square * k ) / math.log(2);

    return u,d,vt,p,c;
   
   
def randomized_stage(A, k, p, c):
    trials = int(ceil(a));
    for t in xrange(trials):
        i = roulette_pick(p); 
       


#def deterministic_stage(M):

