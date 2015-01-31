#!/bin/python
import os;
import sys;

path      = os.path.split(os.path.realpath(__file__))[0];
sys.path.append(path + "/../utils/Python_Utils");
sys.path.append(path + "/../Python_Utils");

from Matrix_Utils import *;
import numpy      as np;
import algorithm1 as greedy;
import Logger;

def css(A, k, epsilon = 0.1):
    m,n             = A.shape;
    if k > min(m,n):
        raise Exception("k must be <= min(m,n), but k=%d, m=%d, n=%d"%(k,m,n));

    Logger.instance.info("svd starts");
    u,d,vt = np.linalg.svd(A);
    target = u[:,0:k];
    for i in xrange(k):
        target[:,i:i+1] = d[i] * target[:,i:i+1];    
    Logger.instance.info("svd ends");

    Ak   = np.dot(target, vt[0:k,:]);
    norm = np.linalg.norm(A - Ak, 'fro'); 

    isPass, Lambda = greedy.greedy(A, target, epsilon * norm);
    C   = [];
    for i in xrange(len(Lambda)):
        if True == Lambda[i]:
            C.append(i);    
    C   = np.array(C);
    return C;
