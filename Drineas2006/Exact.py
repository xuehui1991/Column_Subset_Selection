#!/bin/python2
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(util_dir1);

from Matrix_Utils import *;
from Roulette     import *;
import math;
import numpy as np;

is_debug = False;

def eval(A, C):
    PC   = matrix_A_dot_PI(A, C);
    diff = A - np.dot(np.dot(PC, np.linalg.pinv(PC)), A);
    return np.linalg.norm(diff, 'fro');

def css(R, k, epsilon, delta):
    m,n = R.shape;
    if k > n:
        raise Exception("CSSP requires k <= n, but k=%d, n=%d"%(k,n));    

    num  = int(ceil(math.log(1/delta)));
    Cmin = [];
    Vmin = 10000000;
    c    = k * k * math.log( 1/delta) /epsilon / epsilon;
    for iter1 in xrange(num):

        if is_debug == True:
            print "iter %d >>>>>>>>",iter1;
    
        u,d,vt = np.linalg.svd(R);
        p      = array([0.0 for j in xrange(n)]);
        for j in xrange(n):
            p[j]  = np.linalg.norm(vt[0:k, j:j+1]);
            p[j] *= p[j];
            p[j] /= k;
        
        C     = [];
        while len(C) < c:
            i = roulette_pick(p);
            C.append(i);
        C  = np.array(C);
    
        value = eval(R, C);
        if value < Vmin:
            Vmin = value;
            Cmin = C;            


    return Cmin;
