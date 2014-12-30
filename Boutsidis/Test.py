#!/bin/python2
from numpy import *;
from Boutsidis import *;
import numpy as np;

if __name__ == "__main__":
    testM = array([[1,2,3,4],[0,2,3,4],[0,0,3,4],[0,0,3,4]]);
    testM = testM*5;
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
