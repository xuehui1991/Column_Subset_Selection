#!/bin/python2
import sys;
import os;
util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
sys.path.append(util_dir);
from Matrix_Utils import *;
import Roulette;

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

testM = array([[1,2,3,4],[0,2,3,4],[0,0,3,4],[0,0,3,4]]);
testM = testM*5;
m,n   = testM.shape;
k     = 2;
u,d,vt,p,c = initial_stage(testM,k);
print "u:"
matrix_show(u);
print "d:";
matrix_show(d);
print "vt:";
matrix_show(vt);
print "p:";
matrix_show(p);
print "c:",c;

vkts1d1,s1,d1 = randomized_stage(testM, k, vt, p, c);
print "vkts1d1:";
matrix_show(vkts1d1);
print "s1";
matrix_show(s1);
print "d1:";
matrix_show(d1);

sigma = diag(d);
sigma_p_minus_k = sigma[k:m, k:n];
vt_p_minus_k    = vt[k:n,:];
m               = dot(sigma_p_minus_k, vt_p_minus_k);
value1          = linalg.norm(m,'fro');
print value1*value1;

sum_value = 0;
for i in xrange(10000000):
    vkts1d1,s1,d1 = randomized_stage(testM, k, vt, p, c);
    S1            = zeros([n,c]);
    D1            = diag(d1);
    for j in xrange(c):
        S1[s1[j],j] = 1.0;
    #print "p:";
    #matrix_show(p);
    #print "S1:";
    #matrix_show(S1);
    #print "D1:";
    #matrix_show(D1);
    tmp = dot(dot(m,S1),D1);
    v1  = linalg.norm(tmp,'fro');
    sum_value    += v1 * v1;
    print i+1,value1*value1, sum_value/(i+1);
