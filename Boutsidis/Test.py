#!/bin/python2
import os;
import sys;
util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
sys.path.append(util_dir);


from Boutsidis import *;
from Matrix_Utils import *;
import numpy as np;
np.random.seed(0);

if __name__ == "__main__":
    testM = array([[1,2,3,4],[0,2,3,4],[0,0,3,4],[0,0,3,4]]);
    testM = np.random.random([100,100]);
    k     = 2;
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
    u,d,vt,p,c =  initial_stage(testM, k);
    print "p:";
    print p;
    print "c:";
    print c;

    print "Randomized Stage:";
    vkts1d1,s1,d1 = randomized_stage(testM, k, vt, p, c);
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
    s2 =  deterministic_stage(vkts1d1, k);
    print "s2.shape:";
    matrix_show(s2);

    print "css >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
    C = css(testM,1);
    matrix_show(C);
        
