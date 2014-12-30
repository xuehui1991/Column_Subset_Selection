#!/bin/python
from numpy import *;
import random;
import math;

#This is a toy code for theorem2 in "An improved Approximation Algorithm for the Column Subset Selection Problem"
def exactly(A,c):
    #p[0,n] = 1/n
    m,n = A.shape;
    p = array([1.0/n for i in xrange(n)]);
    
    S = zeros(n,c);
    for t in xrange(c):
        i = pick(p);
        S[i,t] = 1/math.sqrt(c * p[i]);

    C = dot(A,S);

def pick(p):
    v = random.random();
    s = 0;
    for i in xrange(len(p)):
        s += p[i];
        if s[i] >= v:
            return i;
    return len(p)-1;

if __name__ == "__main__":
    
