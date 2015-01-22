#!/bin/python2

import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(util_dir1);

from Float_Utils  import *;
from Matrix_Utils import *;
from Test_Utils   import *;
import numpy as np;
import math;

np.random.seed(0);
is_debug = False;

def choose(A, B, Lambda):
    am,an = A.shape;
    max_i = -1;
    max_v = -1;
    for i in xrange(an):
        if 1 == Lambda[i]:  continue;
        BtAi  = np.dot(np.transpose(B),A[:,i:i+1]);
        sigma = np.linalg.norm(BtAi); 
        if sigma > max_v:
            max_v = sigma;
            max_i = i;
    return max_i,max_v;

def choose_optimized(D):
    max_p = -1;
    max_v = -1;
    m,n   = D.shape;
    for j in xrange(n):
        value = 0.0;
        for i in xrange(m):
            value += D[i,j] * D[i,j];
        value = math.sqrt(value);
        if value > max_v:
            max_p = j;
            max_v = value;
    return max_p, max_v;

def normalize_column(A, Lambda = None,selected = None):
    m,n = A.shape;
    if selected == None:    selected = 1;
    for j in xrange(n):
        if None != Lambda and Lambda[j] != selected:    continue;
        s = 0.0; 
        for i in xrange(m):
            s += A[i,j] * A[i,j];
        if eq(s,0.0):   continue;
        s = math.sqrt(s);
        for i in xrange(m):
            A[i,j] = A[i,j] / s;

    return A;

def greedy(A, B, delta):
    am,an   = A.shape;
    bm,bn   = B.shape;

    if am != bm:
        raise Exception("The greedy algorithm requires A.row = B.row, \
                        but A.row = %d B.row = %d"%(am,bm));
    Lambda  = np.array([0 for j in xrange(an)]);
    is_pass = True;

    A = normalize_column(A);

    B_fro  = np.linalg.norm(B,'fro');
    if is_debug:
        iter1 = 0;
        print "iter=%d:>>>>>>>>>>>>>>>>>>>>>>>>>>>"%iter1
        print "B_fro=%f"%(B_fro);
        print "A"
        matrix_show(A);
        print "B"
        matrix_show(B);
        print "Lambda";
        matrix_show(Lambda);
        print ""


    while B_fro > delta:
        max_j,max_v = choose(A, B, Lambda);
        if -1 == max_j: 
            is_pass = False;
            break;

        Lambda[max_j] = 1;
        for j in xrange(an):
            if 0 != Lambda[j]: continue;
            A[:,j:j+1] -= np.dot(np.transpose(A[:,j:j+1]), A[:,max_j:max_j+1]) * A[:,max_j:max_j+1]  
        for j in xrange(bn):
            B[:,j:j+1] -= np.dot(np.transpose(B[:,j:j+1]), A[:,max_j:max_j+1]) * A[:,max_j:max_j+1];

                
        A = normalize_column(A, Lambda, 0);
        B_fro = np.linalg.norm(B,'fro');   
        if is_debug:
            iter1 += 1;
            print "iter=%d >>>>>>>>>>>>>>>>>>>>>>>>"%iter1;
            print "|BA|col = %f"%max_v;
            print "A";
            matrix_show(A);
            print "B";
            matrix_show(B);
            print "Lambda";
            matrix_show(Lambda);
            print "B_fro=%f\n"%(B_fro);
 
    if is_pass:
        return True, Lambda 
    else:
        return False, Lambda;

def check(A, B, Lambda, sigma, is_pass):
    print "check start>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    m,n = A.shape;
    num = 0;
    for j in xrange(len(Lambda)):
        num  += Lambda[j];

    C   = zeros([m,num]);
    c   = 0;
    for j in xrange(len(Lambda)):
        if 0 != Lambda[j]:
            C[:,c:c+1] = np.copy(A[:,j:j+1]);
            c += 1;

    X    = np.dot(matrix_pinv(C), B);
    diff = np.dot(C,X) - B;
    print "|CX-B|";
    matrix_show(diff);    
    norm_diff = np.linalg.norm(diff,'fro');
    print "norm_diff:",norm_diff

    if (norm_diff <= sigma and is_pass) or (norm_diff > sigma and False == is_pass): 
        print "check end"
        return True;
    else:
        print "check end"
        return False;

def test():
    global is_debug;
    is_debug = True;
    
    test_start_show();
    print "test choose";
    A = np.array([[1.0,2.0,3.0,4.0],[0.0,2.0,3.0,4.0],[0.0,0.0,3.0,4.0]]);
    B = np.copy(A);
    Lambda = array([0 for i in xrange(4)]);
    i,value = choose(A,B,Lambda);
    print "A";
    matrix_show(A);
    print "B";
    matrix_show(B);
    print "choice:",i;
    if i != 3:
        print "Test fails.";
        exit(1);  
    test_end_show();

    test_start_show();
    A = np.array([[1.0,2.0,3.0,4.0],[0,2,3,4],[0,0,3,4]]);
    B = np.copy(A);
    Ao = np.copy(A);
    Bo = np.copy(B);
    is_pass, Lambda = greedy(A,B,0.001);
    if not check(Ao,Bo, Lambda, 0.001, is_pass):
        print "Test fails. ";
        exit(1);
    test_end_show();

    
    test_start_show();
    A = np.array([[1.0,2.0,3.0,4.0],[0,2,3,4],[0,0,3,4]]);
    B = np.copy(A[:,1:2]);
    Ao = np.copy(A);
    Bo = np.copy(B);
    is_pass, Lambda = greedy(A,B,0.001);
    if not check(Ao,Bo, Lambda, 0.001, is_pass):
        print "Test fails. ";
        exit(1);
    test_end_show();


    for iter1 in xrange(1,15):
        test_start_show();
        A  = np.random.random([iter1,iter1%6+1]);
        B  = np.random.random([iter1,iter1%4+1]);
        Ao = np.copy(A);
        Bo = np.copy(B);
        is_pass, Lambda = greedy(A,B,0.001);
        if not check(Ao,Bo, Lambda, 0.001, is_pass):
            print "Test fails. ";
            exit(1);
        test_end_show();

if __name__ == "__main__":
    test();
