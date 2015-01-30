#!/bin/python2
import os;
import sys;

path      = os.path.split(os.path.realpath(__file__))[0];
sys.path.append(path + "/../utils/Python_Utils");
sys.path.append(path + "/../Python_Utils");


from Roulette import *;
import numpy as np;
import Logger;

def css(R,k):
    m,n = R.shape;
    if k > n:
        raise Exception("Weibi Algorithm requires k <= n, but k=%d, n=%d"%(k,n));    

    Logger.instance.info("Algorithm starts");
    Logger.instance.info("SVD starts");
    u,d,vt = np.linalg.svd(R);
    vtk    = vt[0:k, :];
    Logger.instance.info("SVD ends");
    
    p      = array([0.0 for j in xrange(n)]);
    for j in xrange(n):
        p[j]  = np.linalg.norm(vt[0:k,j:j+1]);
        p[j] *= p[j];
        p[j] /= k;
    Logger.instance.info("Probabilities calculation ends");
     
    C     = [];
    Cdict = [0 for j in xrange(n)];
    while len(C) < k:
        i = roulette_pick(p);
        if 0 ==  Cdict[i]:
            Cdict[i] = 1;
            C.append(i);
    Logger.instance.info("Algorithm ends");

    return np.array(C); 
