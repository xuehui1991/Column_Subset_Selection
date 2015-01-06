#!/bin/python2
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
rrqr_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/RankRevealing_QR_Factorization/MingGu_StrongRRQR";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
rrqr_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../RankRevealing_QR_Factorization/MingGu_StrongRRQR";
sys.path.append(util_dir);
sys.path.append(rrqr_dir);
sys.path

from numpy import *;
import Roulette;
import StrongRRQR;


def initial_stage(A, k):
    row,col = A.shape;
    u,d,vt  = linalg.svd(A); 
    p       = array([0.0 for i in xrange(col)]);
    p1      = array([0.0 for i in xrange(col)]);
    for j in xrange(col):
        for i in xrange(k):
            p1[j] +=  vt[i,j]*vt[i,j];
        p1[j] /= k;

    d_lou_minus_k  = diag(d)[k:row,k:col];
    vt_lou_minus_k = vt[k:row,:];
    m              = dot(d_lou_minus_k, vt_lou_minus_k);
    m_f            = linalg.norm(m,'fro');
    m_f            = m_f * m_f;
    p2             = array([0.0 for i in xrange(col)]);
  
    for j in xrange(col):
        for i in xrange(k,row):
            p2[j] += m[i-k,j] * m[i-k,j];
        p2[j] /= m_f;
 
    for j in xrange(col):
        #p[j] = p1[j]/2 + p2[j]/2;    
        p[j]  = p1[j];

    Af        = linalg.norm(A,'fro');
    c0_square =  1.0
    c         =  1600 *  c0_square  *  k * math.log( 800 * c0_square * k ) / math.log(2);
    c         =  int(ceil(c));


    return u,d,vt,p,c;
   
   
def randomized_stage(A, k, vt, p, c):
    row,col = A.shape;
    S1 = array([0.0 for j in xrange(c)]); 
    D1 = array([0.0 for j in xrange(c)]);

    vkts1d1 = array([[0.0 for j in xrange(c)] for i in xrange(k)]);
    for t in xrange(c):
        i = Roulette.roulette_pick(p); 
        S1[t] = i;
        D1[t] = 1.0 / math.sqrt( c * p[i] );
        for r in xrange(k):
            vkts1d1[r,t] = vt[r,i] * D1[t];

    return vkts1d1, S1, D1;


def deterministic_stage(vts1d1,k):
    row,col      = vts1d1.shape;
    hat          = zeros([col,col]);
    hat[0:row,:] = vts1d1;
    R,S2         = StrongRRQR.rrqr(hat,k);

    return S2[0:k];


def css(A, k):
    m,n = A.shape;    

    u,d,vt,p,c      =  initial_stage(A, k);
    vkts1d1, s1, d1 =  randomized_stage(A, k, vt, p, c);
    s2              =  deterministic_stage(vkts1d1, k);
    
    C               =  zeros([1,k]);
    for i in xrange(k):
        C[i] = s1[s2[i]];
    return C;
                

