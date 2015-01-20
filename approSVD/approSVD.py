#!/bin/python
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(rrqr_dir);

import numpy      as np;
import algorithm1 as greedy;
from Matrix_Utils import *;


def css(A, k, delta=0.1):
    m,n    = A.shape;

    v,d,vt = np.linalg.svd(A);
    Ak     = np.dot(u[:,0:k], vt[0:k,:]);
    isPass, Lambda = greedy(A, Ak, delta);
    
    C   = [];
    for i in xrange(len(Lambda)):
        if True == Lambda[i]:
            C.append(i);    
    C   = np.array(C);
    return C;
