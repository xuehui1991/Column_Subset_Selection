#!/bin/python2
import os;
import sys;

util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
rrqr_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/RankRevealing_QR_Factorization/MingGu_StrongRRQR";
sys.path.append(util_dir);
sys.path.append(rrqr_dir);

from numpy import *;
from Roulette import *;
from StrongRRQR import *;

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
    c         =  int(ceil(c));

    return u,d,vt,p,c;
   
   
def randomized_stage(A, k, vt, p, c):
    row,col = A.shape;
    S1 = array([[0.0 for j in xrange(c)] for i in xrange(col)]);
    D1 = array([[0.0 for j in xrange(c)] for i in xrange(col)]); 
         

    vkts1d1 = array([[0.0 for j in xrange(c)] for i in xrange(k)]);
    for t in xrange(c):
        i = roulette_pick(p); 
        S1[i,t] = 1.0;
        D1[t,t] = 1.0 / math.sqrt( c * p[i] );
        for r in xrange(k):
            vkts1d1[r,t] = vt[r,i] * D1[t,t];


    return vkts1d1, S1, D1;


def deterministic_stage(vts1d1,k):
    StrongRRQR.rrqr(vts1d1,k);
    

