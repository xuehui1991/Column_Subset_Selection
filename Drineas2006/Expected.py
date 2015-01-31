#!/bin/python2
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(util_dir1);

from Roulette     import *;
from Matrix_Utils import *;
import math;
import numpy as np;

def eval(A, C):
    PC   = matrix_A_dot_PI(A, C);
    diff = A - np.dot(np.dot(PC, np.linalg.pinv(PC)), A);
    return np.linalg.norm(diff, 'fro');

def css(R, k, epsilon, delta = 0.5): ##default delta = 0.5, num = 1;
    m,n = R.shape;
    if k > n:
        raise Exception("CSSP requires k <= n, but k=%d, n=%d"%(k,n));    

    num  = int(ceil(math.log(1/delta)));
    Cmin = [];
    Vmin = 10000000;
    c    = k * k * math.log( 1/delta ) /epsilon / epsilon;
    for iter1 in xrange(num):
    
        u,d,vt = np.linalg.svd(R);
        p      = array([0.0 for j in xrange(n)]);
        for j in xrange(n):
            p[j]  = np.linalg.norm(vt[0:k, j:j+1]);
            p[j] *= p[j];
            p[j] /= k;
        
        t     = 0;
        C     = [];
        Cdict = [0 for j in xrange(n)];
        while t < c:
            t += 1;
            i  = roulette_pick(p);
            if 1 == Cdict[i]: continue;
            C.append(i);
            Cdict[i] = 1;
        C = np.array(C);

        if 1 == num:
            Cmin = C;
        else:
            value = eval(R, C);
            if value < Vmin:
                Vmin = value;
                Cmin = C;         

    return Cmin;
