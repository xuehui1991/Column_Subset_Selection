#!/bin/python
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(util_dir1);


import numpy      as np;
import algorithm1 as greedy;
from Matrix_Utils import *;

is_debug = False;

def css(A, k, delta = 0.1):
    greedy.is_debug = is_debug;
    m,n             = A.shape;
    if k > min(m,n):
        raise Exception("k must be <= min(m,n), but k=%d, m=%d, n=%d"%(k,m,n));

    u,d,vt = np.linalg.svd(A);
    target = u[:,0:k];
    for i in xrange(k):
        target[:,i:i+1] = d[i] * target[:,i:i+1];
 
    isPass, Lambda = greedy.greedy(A, target, delta);
    
    C   = [];
    for i in xrange(len(Lambda)):
        if True == Lambda[i]:
            C.append(i);    
    C   = np.array(C);
    return C;
