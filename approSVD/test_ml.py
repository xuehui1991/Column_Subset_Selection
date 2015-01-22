#!/bin/python2
import os;
import sys;
util_dir = os.path.split(os.path.realpath(__file__))[0]+"/../utils/Python_Utils";
sys.path.append(util_dir);


import approSVD  as approsvd;
from Matrix_Utils import *;
from Test_Utils   import *;
import numpy as np;
np.random.seed(0);
import random;
random.seed(0);

def readdata(filename = "/home/ll/multi_label_data/delicious/train.y"):
    f = open(filename, "r");
    A = [];
    for line in f:
        row  = [];
        eles = line.strip().split(" ");
        for i in xrange(len(eles)):
            row.append(float(eles[i])); #not row.append(int(eles[i]));
        A.append(row);
    return np.array(A);

if __name__ == "__main__":
    approsvd.is_debug = True;

    test_start_show();
    A = readdata(); 
    k = 100;
    C = approsvd.css(A,k);
 
    Ak = matrix_Ak(A,k);
    matrix_show(Ak);
    error = A - Ak;
    error_norm = np.linalg.norm(error, 'fro');
    print "|A-Ak|",error_norm;

    AC = matrix_A_dot_PI(A,C);
    PC = np.dot(AC,matrix_pinv(AC));
    error = A - dot(PC,A);
    error_norm = np.linalg.norm(error, 'fro');
    print "|A-PcA|F", error_norm;
    test_end_show();

