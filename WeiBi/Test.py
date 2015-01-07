#!/bin/python2
import os;
import sys;

util_dir  = os.path.split(os.path.realpath(__file__))[0] + "/../utils/Python_Utils";
util_dir1 = os.path.split(os.path.realpath(__file__))[0] + "/../Python_Utils";
sys.path.append(util_dir);
sys.path.append(util_dir1);

from Test_Utils import *;
from Matrix_Utils import *;
from WeiBi import *;
import numpy as np;

def test():
    test_start_show();
    k =2;
    A = np.array([[1,2,3],[0,2,3],[0,0,3]]);
    print "A";
    matrix_show(A);
    print "k=%d"%k;
    C = css(A,k);
    print "C:";
    matrix_show(C);

    Ak = matrix_Ak(A,k);
    print "A%d"%k;
    matrix_show(Ak);
    error = A - Ak;
    fro   = np.linalg.norm(error,'fro');
    print "fro |A-Ak|";
    print fro;    

    AC = matrix_A_dot_PI(A,C);
    PC = np.dot(AC,matrix_pinv(AC));
    error = A - dot(PC,A);
    fro   = np.linalg.norm(error,'fro');
    print "fro |A-PcA|";
    print fro;
    test_end_show();   

if __name__ == "__main__":
    test();    
