#!/bin/python2
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(util_dir1);

from Roulette import *;
import numpy as np;

def css(R,k):
    m,n = R.shape;
    if k > n:
        raise Exception("Weibi Algorithm requires k <= n, but k=%d, n=%d"%(k,n));    

    u,d,vt = np.linalg.svd(R);
    vtk    = vt[0:k, :];
    
    p      = array([0.0 for j in xrange(n)]);
    for j in xrange(n):
        for i in xrange(k):
            p[j] += vtk[i,j] * vtk[i,j];
        p[j] /= k;
    
    
    C     = [];
    Cdict = [0 for j in xrange(n)];
    while len(C) < k:
        i = roulette_pick(p);
        if 0 ==  Cdict[i]:
            Cdict[i] = 1;
            C.append(i);
    return array(C); 
