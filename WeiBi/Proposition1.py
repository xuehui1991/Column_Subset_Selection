#!/bin/python2
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(util_dir1);

from Test_Utils import *;
from Matrix_Utils import *;
from Roulette import *;
from Float_Utils import *;
import numpy as np;
np.random.seed(0);

def vtk_c(R,k):
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

    return vtk, array(C); 


def test():
    test_start_show();
    k = 10;
    m = 100;
    n = 100;

    s = 0;
    for iter1 in xrange(1000):
        A = np.random.random([m,n]);
        #print "A:";
        #matrix_show(A);
        vtk, C = vtk_c(A, k);
        VtkC   = matrix_A_dot_PI(vtk, C);
        #print "VtkC:";
        #matrix_show(VtkC);

        det = np.linalg.det(VtkC);
        if not eq(det,0.0):
            s += 1;

        print "iter=%d,\tdet=%f\t,prob_of_fullrank=%f"%(iter1,det, s*1.0/(iter1+1));


    test_end_show();   

if __name__ == "__main__":
    test();    
