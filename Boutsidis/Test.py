#!/bin/python2
import os;
import sys;
util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
sys.path.append(util_dir);


import Boutsidis as Bs;
from Matrix_Utils import *;
from Test_Utils   import *;
import numpy as np;
np.random.seed(0);
import random;
random.seed(0);

if __name__ == "__main__":
    Bs.is_debug = True;

    test_start_show();
    testM = array([[1,2,3,4],[0,2,3,4],[0,0,3,4],[0,0,3,4]]);
    testM = np.random.random([100,100]);
    k     = 10;
    u,d,vt = linalg.svd(testM);
    print "testM:";
    matrix_show(testM);
    print "u:"
    matrix_show(u);
    print "d:";
    matrix_show(d);
    print "vt:";
    matrix_show(vt);

    print "Initialization Stage:";
    u,d,vt,p,c =  Bs.initial_stage(testM, k);
    print "p:";
    matrix_show(p);
    print "c:";
    print c;

    print "Randomized Stage:";
    vkts1d1,s1,d1 = Bs.randomized_stage(testM, k, vt, p, c);
    print "vkts1d1.shape";
    f = open("data","w");
    row,col = vkts1d1.shape;
    for i in xrange(row):
        for j in xrange(col):
            f.write("%f\t"%vkts1d1[i,j]);
        f.write("\n");
    f.close();
    print vkts1d1.shape;
    print "s1.shape";
    print s1.shape;
    print "d1.shape";
    print d1.shape;

    print "Deterministic Stage:";
    s2 =  Bs.deterministic_stage(vkts1d1, k);
    print "s2.shape:";
    matrix_show(s2);
    test_end_show();

    test_start_show();
    k = 10;
    A = np.copy(testM);
    C = Bs.css(A,k);
 
    Ak = matrix_Ak(A,k);
    print "A%d"%k;
    matrix_show(Ak);
    error = A - Ak;
    error_norm = np.linalg.norm(error, 'fro');
    print "|A-Ak|",error_norm;

    AC = matrix_A_dot_PI(A,C);
    PC = np.dot(AC,matrix_pinv(AC));
    error = A - dot(PC,A);
    error_norm = np.linalg.norm(error, 'fro');
    print "|A-PcA|F", error_norm;
    
    test_end_show();
