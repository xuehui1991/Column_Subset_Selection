#!/bin/python
from numpy import *;
import HouseHolder;
#Efficient Algorithm for Computing A Strong Rank-Relealing QR Factorization  850 page 2.RRQR Algorithm   

def test():
    M = array([(1,2,3),(2,2,4),(4,5,4)]);
    print M;

    m = copy(M);
    [r,p] = rrqr(m,0.1);
    print p;
    print r;

    testM = array([[1,0,1,1,0],[0,1,1,1,0],[1,1,1,0,0],[1,0,1,0,0],[0,0,0,0,1]]);
    print testM;
    [r,p] = rrqr(testM,0.001);
    print p;
    print r;

def rrqr(M,sigma):
    (r,c) = M.shape;
    if r < c:
        raise Exception("M.r(%d) < M.c(%d)"%(r,c))

    PI   = eye(c);
    PI1D = range(c);
    R    = M;
    
    S    = [0 for col in xrange(c)];
    for col in xrange(c):
        for row in xrange(r):
            S[col] += M[row,col] * M[row,col];

    for k in xrange(c):
        j,value = maxOf(S,k)
        if value < sigma:   break

        #swap j and k
        PI[PI1D[j],j] = 0;
        PI[PI1D[k],k] = 0;
        tmp        = PI1D[k];
        PI1D[k]    = PI1D[j];
        PI1D[j]    = tmp;  
        PI[PI1D[k],k] = 1;
        PI[PI1D[j],j] = 1;

        tmp  = S[j];
        S[j] = S[k];
        S[k] = tmp;

        R[:,[j,k]] = R[:,[k,j]];

        R = HouseHolder.HouseHolder_step(R, k);        
        for i in xrange(k,len(S)):
            S[i] -= R[k,i] * R[k,i];

    return R,PI

def maxOf(S,k):
    max_p = k;
    max_v = S[k];
    for i in xrange(k+1,len(S)):
        if S[i] > max_v:
            max_p = i;
            max_v = S[i];
    return max_p,max_v;

if __name__ == "__main__":
    test();
