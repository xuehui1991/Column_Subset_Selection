#/bin/python
from numpy import *;

def is_matrix_equals(M1,M2):
    if 2 == M1.ndim: 
        if M1.shape != M2.shape:
            return False;
        r,c = M1.shape;
        for i in xrange(r):
            for j in xrange(c):
                if abs(M1[i,j] - M2[i,j]) >= 1e-9:
                    return False;
        return True;
    elif 1 == M1.ndim:
        if M1.shape != M2.shape:
            return False;
        r = M1.shape[0];
        for i in xrange(r):
            if abs(M1[i] - M2[i]) >= 1e-9:
                return False;
        return True;

    else:
        raise Exception("equals_matrix function not support ndim = %d"%(M1.ndim));
    

def matrix_show(M):
    if 2 == M.ndim:
        r,c = M.shape;
        for i in xrange(r):
            for j in xrange(c):
                if M[i,j] >= 0:
                    print " %.3f\t"%M[i,j],
                else:
                    print "%.3f\t"%M[i,j],
            print "";
    elif 1 == M.ndim:
        l = len(M);
        for i in xrange(l):
            if M[i] >= 0:
                print " %.3f\t"%M[i],
            else:
                print "%.3f\t"%M[i],
        print "";
    else:
        raise Exception("show_matrix not support ndim=%d yet"%M.ndim);

