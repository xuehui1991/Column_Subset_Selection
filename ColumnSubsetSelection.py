#!/bin/python
import sys;
import os;
path      = os.path.split(os.path.realpath(__file__))[0];
sys.path.append(path + "/utils");
sys.path.append(path + "/utils/Python_Utils");
sys.path.append(path + "/utils/Python_Utils");
sys.path.append(path + "/WeiBi");
sys.path.append(path + "/Drineas2006");
sys.path.append(path + "/approSVD");
sys.path.append(path + "/Eval");

from Matrix_Utils import *;
from Test_Utils   import *;
import numpy as np;
np.random.seed(0);
import random;
random.seed(0);
import Logger;


type_algorithm_str   = [["" for i in xrange(1)],\
                        ["" for i in xrange(2)]];
type_algorithm_str[0][0] = "WeiBi et al, ICML2013, Efficient Multi-label Classification with Many Labels"
type_algorithm_str[1][0] = "Drineas et al, 2006, Subspace Sampling and Relative Error Approximation: Column-based Method, Expect Algorithm";
type_algorithm_str[1][1] = "A.Civril and M.Magdon-Ismail, 2012, Column Subset Selection via Sparse Aprroximation of SVD";

def printUsages():
    print "Usage: ColumnSubsetSelection.py [options] matrix_file output_file";
    print "options:"
    print "-t cssp_type : set type of column subset selection (default 0)";
    print "      0 -- NS (number of columns is specified)";
    print "      1 -- ES (epsilon is specified)";

    print "-a algorithm : set algorithm (default 0)"
    print "   for number of columns specified (NS) ";
    print "      0 -- %s"%type_algorithm_str[0][0];
    print "   for epsilon is specified (ES)";
    print "      0 -- %s"%type_algorithn_str[0][0];
    print "      1 -- %s"%type_algorithm_str[1][1];
    print "-k the number of columns selected";
    print "-e the epsilon";
    print "-log log_file (if not specified, write log information to stderr)";
    print "-level log_level (default the info level)";
    print "      debug -- the debug level";
    print "      info  -- the info level";

       

def parseParameter(args):
    if len(args) < 3: #at least 3 paramters: ColumnSubsetSelection matrix_file output_file
        printUsages();
        exit(1);

    matrix_filename = args[len(args) - 2];
    output_filename = args[len(args) - 1];
    
    dicts      = dict();
    dicts["t"] = 0;
    dicts["a"] = 0;
    i = 1;
    while i + 1 < len(args) - 2:
        if "-t" == args[i]:
            dicts["t"] = int(args[i+1]);
        elif "-a" == args[i]:
            dicts["a"] = int(args[i+1]);
        elif "-k" == args[i]:
            dicts["k"] = int(args[i+1]);
        elif "-e" == args[i]:
            dicts["e"] = float(args[i+1]);
        else:
            printUsages();
            exit(1);
        i += 2;

    if 0 == dicts["t"] and "k" not in dicts:
        Logger.instance.error("we must provide k for the NS type");
    if 1 == dicts["t"] and ("k" not in dicts or "e" not in dicts):
        Logger.instance.error("we must provide k and epsilon for the ES type");

    return matrix_filename, output_filename, dicts;

def read_matrix(filename):
    f = open(filename, "r");
    A = [];
    for line in f:
        eles = line.strip().split(" ");
        row  = map(float, eles);
        A.append(row);
    return np.array(A); 

def write_result(C, filename):
    f = None;
    try:f = open(filename, "w");
    except:
        logger.waring("fail to open %s as output_file, \
        store the result to ./result"%filename);
        f = open("./result", "w");

    for j in xrange(len(C)):
        f.write("%d "%C[j]);
    f.write("\n");


def css(A, opts):
    C = [];
    if   0  == opts["t"] and 0  == opts["a"]:
        import WeiBi;
        k  = opts["k"];
        C  = WeiBi.css(A, k);
    elif 1  == opts["t"] and 0  == opts["a"]:
        import Expected;
        k  = opts["k"];
        e  = opts["e"];
        C  = Expected.css(A, k, e);
    elif 1  == opts["t"] and 1  == opts["a"]:
        import approSVD;
        k  = opts["k"];
        e  = opts["e"];
        C  = approSVD.css(A, k, e);
    else:
        Logger.instance.error("Not supported type=%d and algorithm=%d"%(opts["t"], opts["a"]));
        exit(1);
    return C;

if __name__ == "__main__":
    matrix_filename, output_filename, opts = parseParameter(sys.argv);
    Logger.initlog(opts);
    Logger.instance.info("k = %d"%opts["k"]);
    if "e" in opts:
        Logger.instance.info("epsilon = %f"%opts["e"]);
    Logger.instance.info("matrix_file [ %s ]"%matrix_filename);
    Logger.instance.info("output_file [ %s ]"%output_filename);
    if not os.path.exists(matrix_filename):
        Logger.instance.error("matrix_file [ %s ] not exists"%matrix_filename);
        exit(1);


    Logger.instance.info("data reading starts");
    A = read_matrix(matrix_filename);
    Logger.instance.info("Column Subset Selection starts");
    C = css(A, opts);
    write_result(C, output_filename);
    Logger.instance.info("Column Subset Selection completes. ");
    
