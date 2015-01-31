#!/bin/python2
import os;
import sys;

path      = os.path.split(os.path.realpath(__file__))[0];
sys.path.append(path + "/../utils/Python_Utils");
sys.path.append(path + "/../Python_Utils");

from Matrix_Utils import *;
from Roulette     import *;
import math;
import numpy as np;
import Logger;


def eval(A, C):
    PC   = matrix_A_dot_PI(A, C);
    diff = A - np.dot(np.dot(PC, np.linalg.pinv(PC)), A);
    return np.linalg.norm(diff, 'fro');

def css(R, k, epsilon, delta = 0.5):
    m,n = R.shape;
    if k > n:
        raise Exception("CSSP requires k <= n, but k=%d, n=%d"%(k,n));    

    num  = int(ceil(math.log(1/delta)));
    Cmin = [];
    Vmin = 10000000;
    c    = k * k * math.log( 1/delta) /epsilon / epsilon;
    for iter1 in xrange(num):

        Logger.instance.debug("iter %d >>>>>>>>",iter1);
    
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
    
        if 1 == num:
            Cmin = 1;
        else:
            value = eval(R, C);
            if value < Vmin:
                Vmin = value;
                Cmin = C;            

    return Cmin;
