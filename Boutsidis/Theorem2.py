#!/bin/python

import os;
import sys;

util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
rrqr_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/RankRevealing_QR_Factorization/MingGu_StrongRRQR";
sys.path.append(util_dir);
sys.path.append(rrqr_dir);

import numpy as np;
np.random.seed(0);
from Roulette import *;
import random;
random.seed(0);
import math;

#This is a toy code for theorem2 in "An improved Approximation Algorithm for the Column Subset Selection Problem"

def exactly(A,c,p):
    m,n = A.shape;   
    pick_i = 0; 

    S = zeros([n,c]);
    for t in xrange(c):
        i = roulette_pick(p);
        #i = pick_i%n;
        S[i,t]  = 1/math.sqrt(c * p[i]);
        #pick_i += 1;
        #pick_i %= n;

    C = dot(A,S);
    return C;

def calculate_P(A):
    m,n = A.shape;
    A_f = np.linalg.norm(A, 'fro');
    p = array([0.0 for i in xrange(n)]);
    for j in xrange(n):
        for i in xrange(m):
            p[j] += A[i,j] * A[i,j];
        p[j] /= A_f * A_f;
    return p;

if __name__ == "__main__":
    sum_value  = 0;
    sum_value1 = 0;
    
    M  = np.random.random([200,200]);
    #M1 = zeros([200,200]);
    #for x in xrange(200):
    #    M1[x,x] = x;
 
    k = 50;
    print "start svd";
    u,d,v = linalg.svd(M);
    v     = copy(v[0:k,:]);
    print "M svd ends";
    #u1,d1,v1 = linalg.svd(M1);
    #v1    = copy(v1[0:k,:]);
    #print "M1 svd ends";

    print "calculate p";
    p     = calculate_P(v);
    print "p calculate ends";
    #p1    = calculate_P(v1);   
    #print "p1 calculate ends";

    for i in xrange(10000000):

        #c = 2 * 4* math.log(4)/math.log(2);
        #c = int(ceil(c));
        c = 1600 * k * math.log(800 * k) / math.log(2);
        print "c:",c;
        c = int(ceil(c));
        #c = 200000;         
 
        C = exactly(v,c,p);
        #C1 = exactly(v1,c,p1);
    
        m = dot(v,transpose(v)) - dot(C,transpose(C));
        tu,td,tv = linalg.svd(m);
        sum_value += td[0];

        #m = dot(v1,transpose(v1)) - dot(C1,transpose(C1));
        #tu1,td1,tv1 = linalg.svd(m);
        #sum_value1 += td1[0];
       
        print "i:",i,"\tc:",c, "\tvalue:",sum_value/(i+1); #"\tvalue1:",sum_value1/(i+1);
    
