#!/bin/python2
import sys;
import os;
util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
sys.path.append(util_dir);
sys.path.append("./");
from Matrix_Utils import *;
from Boutsidis import *;
import Roulette;


testM = array([[1,2,3,4],[0,2,3,4],[0,0,3,4],[0,0,3,4]]);
testM = testM*5;
m,n   = testM.shape;
k     = 2;
u,d,vt,p,c = initial_stage(testM,k);


sum_value = 0;
for i in xrange(10000000):
    vkts1d1,s1,d1 = randomized_stage(testM, k, vt, p, c);
    u,d,v         = linalg.svd(vkts1d1);
    if d[k-1] > 0.5:
        sum_value    += 1.0;
    print i+1,d[k-1], sum_value/(i+1);
