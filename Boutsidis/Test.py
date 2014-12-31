#!/bin/python2
import os;
import sys;
util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
sys.path.append(util_dir);


from numpy import *;
from Boutsidis import *;
from Matrix_Utils import *;
import numpy as np;


if __name__ == "__main__":
    testM = array([[1,2,3,4],[0,2,3,4],[0,0,3,4],[0,0,3,4]]);
    testM = testM*5;
    k     = 1;
    u,d,vt = linalg.svd(testM);
    print "testM:";
    print testM;
    print "u:"
    print u;
    print "d:";
    print d;
    print "vt:";
    print vt;

    print "Initialization Stage:";
    u,d,vt,p,c =  initial_stage(testM, 2);
    print "p:";
    print p;
    print "c:";
    print c;

    print "Randomized Stage:";
    vkts1d1,s1,d1 = randomized_stage(testM, k, vt, p, c);
    print "vkts1d1";
    matrix_show(vkts1d1);
    print "s1:";
    matrix_show(s1);
    print "d1:";
    matrix_show(d1);

    print "Deterministic Stage:";
    s2 =  deterministic_stage(vkts1d1, k);
    print "s2:";
    matrix_show(s2);
