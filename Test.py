#!/bin/python2
from numpy import *;
import numpy as np;

if __name__ == "__main__":
    testM = array([[1,2,3,4],[0,2,3,4],[0,0,3,4],[0,0,3,4]]);
    testM = testM/10.0
    u,d,vt = linalg.svd(testM);
    print "testM:";
    print testM;
    print "u:"
    print u;
    print "d:";
    print d;
    print "vt:";
    print vt;
