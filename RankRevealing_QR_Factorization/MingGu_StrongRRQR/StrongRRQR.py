#/bin/python
import sys;
from numpy import *;
import numpy as np;
from math  import *;
sys.path.append("../Golub_RRQR");
import Golub_RRQR
import HouseHolder

np.random.seed(0);

Debug=False;

def equals_matrix(M1,M2):
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
    

def show_matrix(M):
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

def gt(f1,f2):
    return f1-f2 > 1e-6;
def lt(f1,f2):
    return f2-f1 > 1e-6;    

def check_final(R,k,f=1.414):
    r,c = R.shape;
    Ak = R[0:k,0:k];
    Bk = R[0:k,k:c];
    Ck = R[k:r,k:c];
    pkn = math.sqrt(1+ f*f*c*(c-k))
    isPass = True;

    print "final check starts";
    u,sigma_R,v = linalg.svd(R);
    u,sigma_Ak,v = linalg.svd(Ak);
    u,sigma_Ck,v = linalg.svd(Ck);
    for i in xrange(k):
        if gt(sigma_Ak[i], sigma_R[i]):
            print "sigma_Ak[%d] > sigma_R[%d]"%(i,i);
            isPass = False;
        if lt(sigma_Ak[i], sigma_R[i] / pkn):
            print "sigma_Ak[%d] < sigma_R[%d]/pkn"%(i,i);
            isPass = False;

    for j in xrange(c-k):
        if lt(sigma_Ck[j], sigma_R[k+j]):
            print "sigma_Ck[%d] < sigma_R[%d]"%(j,k+j);
            isPass = False;
        if gt(sigma_Ck[j], pkn * sigma_R[k+j]):
            print "sigma_Ck[%d] > pkn*sigma_R[%d]"%(j,k+j);
            isPass = False;

    if False == isPass: 
        print "k:",k;
        print "pkn:",pkn;
        print "sigma_R:";
        show_matrix(sigma_R);
        print "sigma_Ak:";
        show_matrix(sigma_Ak);
        print "sigma_Ck:";
        show_matrix(sigma_Ck);
    

    print "final check ends\n";

    return isPass;

def check_step(R, invA_B, omega, gamma, k):
    isPass = True;
    r,c = R.shape;
    Ak = R[0:k,0:k];
    Bk = R[0:k,k:c];
    Ck = R[k:r,k:c];
    
    print "check start: k = ",k;

    invAk = linalg.inv(Ak);
    check_omega = array([0.0 for i in xrange(r)]);
    for i in xrange(k):
        for j in xrange(k):
            check_omega[i] += invAk[i][j] * invAk[i][j];
        check_omega[i] = math.sqrt(check_omega[i]);
    if not equals_matrix(check_omega, omega):
        print "omega:";
        show_matrix(omega);
        print "check_omega:";
        show_matrix(check_omega);
        isPass  = False;

    a = check_gamma = array([0.0 for i in xrange(c)]);
    for j in xrange(k,c):
        for i in xrange(k,r):
            check_gamma[j] += R[i,j] * R[i,j];
        check_gamma[j] = math.sqrt(check_gamma[j]);
    if not equals_matrix(check_gamma, gamma):
        print "gamma:";
        show_matrix(gamma);
        print "check_gamma:";
        show_matrix(check_gamma);
        isPass = False;

    check_invA_B =  dot(linalg.inv(Ak),Bk);    
    if not equals_matrix(invA_B, check_invA_B):
        print "invA_B:";
        show_matrix(invA_B);
        print "check_invA_B:";
        show_matrix(check_invA_B);
        isPass = False;

    print "end check";
    print "";
    return isPass;

def show_step(R, invA_B, omega, gamma, k):
    print "R:"
    show_matrix(R);
    print "gamma:"
    show_matrix(gamma);
    print "omega:"
    show_matrix(omega); 
    print "";



def rrqr(M, k, f=1.414):
    R,PI = Golub_RRQR.rrqr(M,0.0001);
    m,n  = R.shape;
    Ak   = R[0:k,  0:k];
    Bk   = R[0:k,  k:n]; 
    invA = linalg.inv(Ak);
    omega = array([0.0 for i in xrange(m)]);
    for i in xrange(k):
        for j in xrange(k):
            omega[i] += invA[i][j] * invA[i][j];
        omega[i] = math.sqrt(omega[i]);
        
    gamma = array([0.0 for j in xrange(n)]);
    for j in xrange(k,n):
        for i in xrange(k,m):
            gamma[j] += R[i][j] * R[i][j];
        gamma[j] = math.sqrt(gamma[j]);

    invA_B = dot(invA, Bk);

    if True == Debug:
        print "start:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
        show_step(R, invA_B, omega, gamma, k);    
        check_step(R, invA_B, omega, gamma, k);

    flag,i,j = is_rho_less_f(invA_B, omega, gamma, k, f);
    while not flag:
        R, invA_B, omega, gamma =  \
        update_swap_k_kplusj(R, invA_B, omega, gamma, k, j);
        if True == Debug:
            print "swap_k_kplusj: j = ",j,"  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
            show_step(R, invA_B, omega, gamma, k);
            if not check_step(R, invA_B, omega, gamma, k):
                print "Test fails. Please debug this matrix:";
                show_matrix(M);
                exit(0);

        R, invA_B, omega, gamma =  \
        update_shift(R, invA_B, omega, gamma, i, k);
 
        if True == Debug:
            print "shift:i = ",i, "   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
            show_step(R, invA_B, omega, gamma,k);
            if not check_step(R, invA_B, omega, gamma, k):
                print "Test fails. Please debug this matrix:";
                show_matrix(M);
                exit(0);                

        R, invA_B, omega, gamma =  \
        update_swap_kminus1_k(R, invA_B, omega, gamma,k);           
        if True == Debug:
            print "swap_kminus1_k  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
            show_step(R, invA_B, omega, gamma, k);
            if not check_step(R, invA_B, omega, gamma, k):
                print "Test fails. Please debug this matrix:";
                show_matrix(M);
                exit(0);

        flag, i, j = is_rho_less_f(invA_B, omega, gamma, k ,f); 
    
    if True == Debug:
        if not check_final(R, k, f):
            print "Test fails. Please debug this matrix:";
            show_matrix(M);
            exit(0);
    return R;


def is_rho_less_f(invA_B, omega, gamma, k, f=1.414):    
    r,c = invA_B.shape;
    for i in xrange(r):
        for j in xrange(c):
            rho = invA_B[i][j] * invA_B[i][j] + \
                  gamma[j+k]*gamma[j+k] * omega[i] * omega[i];
            if rho > f:
                return False, i, j;

    return True,None,None;
        


def update_swap_k_kplusj(R, invA_B, omega, gamma, k, j): 
    #swap column k and column k + j;  
    row,col = R.shape; 
    if j == 0 or k+j >= col:
        return R, invA_B, omega, gamma,;

    tmp             = copy(R[:,k:k+1]);
    R[:,k:k+1]      = R[:,k+j:k+j+1];
    R[:,k+j:k+j+1]  = tmp;    
    ##employ Given Rotation to ensure Ck[:,0] = [notzero, 0,...0]
    ##Ck[:,0] = [notzero, 0, ..., 0] is important to Update_swap_kminus1_k
    for i in xrange(1,j+1):
        lou = math.sqrt(R[k,k]*R[k,k] + R[k+i,k]*R[k+i,k]);
        c   = R[k,k] / lou;
        s   = -1 * R[k+i,k] / lou;
        for p in xrange(k,col):
            up   = c*R[k,p] - s*R[k+i,p];
            down = s*R[k,p] + c*R[k+i,p];
            R[k,p]   = up;
            R[k+i,p] = down;
        R[k,k] = lou;
       
 
    tmp         = gamma[ k ];
    gamma[k]    = gamma[ k+j ];
    gamma[k+j]  = tmp;

    tmp             =  copy(invA_B[:,j:j+1]);
    invA_B[:,j:j+1] =  invA_B[:,0:1];
    invA_B[:,0:1]   =  tmp;
    
    return R, invA_B , omega, gamma;

def update_shift(R, invA_B, omega, gamma, i, k):
    #4.2
    m,n = R.shape;
    if i >= k:   return R, invA_B, omega, gamma
    tmp = copy(R[0:k,i:i+1]);
    for idx in xrange(i,k-1):
        R[0:k, idx:idx+1] = R[0:k, idx+1:idx+2];
    R[0:k,k-1:k] = tmp;

    
    ## Given Rotation
    for idx in xrange(i,k-1):
        r = math.sqrt(R[idx][idx] * R[idx][idx] + R[idx+1][idx]*R[idx+1][idx]);
        c = R[idx, idx] / r;
        s = -R[idx+1, idx] / r;
        for j in xrange(idx,n):
            up   = c * R[idx][j] - s * R[idx+1][j];
            down = s * R[idx][j] + c * R[idx+1][j];
            R[idx, j]   = up;
            R[idx+1, j] = down;      

    tmp = omega[i];
    for idx in xrange(i,k):
        omega[idx] = omega[idx+1];
    omega[k-1] = tmp;

    tmp = copy(invA_B[i:i+1,:]);
    for idx in xrange(i,k-1):
        invA_B[idx:idx+1,:] = invA_B[idx+1:idx+2,:];
    invA_B[k-1:k,:] = tmp;
    
    return R, invA_B, omega, gamma;

 
def update_swap_kminus1_k(R,invA_B, omega,gamma,k):
    row,col = R.shape;
    r  = R[k-1][k-1];
    print r;
    mu = R[k-1][k] / r;
    nu = R[k][k]   / r;
    lou = math.sqrt(mu * mu + nu * nu);

    test_code='''
    print "before update_swap_kminus1_k";
    print "R:";
    show_matrix(R);
    print "invA_B";
    show_matrix(invA_B);
    print "omega:";
    show_matrix(omega);
    print "gamma:";
    show_matrix(gamma);'''    


    invA_B_r,invA_B_c = invA_B.shape;
    Akminus1 = copy(R[0:k-1, 0:k-1]);
    b1       = copy(R[0:k-1, k-1:k]);
    u        = dot(linalg.inv(Akminus1),b1);
    u1       = invA_B[0:k-1,0:1];
    u2       = transpose(invA_B[k-1:k, 1:invA_B_c]);
    U        = invA_B[0:k-1, 1:invA_B_c];

    test_code='''
    print "Akminus1:";
    show_matrix(Akminus1);
    print "b1";
    show_matrix(b1);
    print "u:";
    show_matrix(u);
    print "u1:";
    show_matrix(u1);
    print "u2:";
    show_matrix(u2);
    print "U:";
    show_matrix(U);'''

    #R and gamma
    tmp             = copy(R[0:k-1, k-1:k])
    R[0:k-1, k-1:k] = R[0:k-1, k:k+1]
    R[0:k-1, k:k+1] = tmp; 
    R[k-1][k-1] = r * lou;
    R[k-1][k]   = r * mu / lou;
    R[k][k]     = r * nu / lou;
    gamma[k]    = abs(R[k,k]);

    for i in xrange(k+1,col):
        c1 = R[k-1][i] * mu / lou + R[k][i] * nu / lou;
        c2 = R[k-1][i] * nu / lou - R[k][i] * mu / lou;
        gamma[i]  = math.sqrt(gamma[i] * gamma[i] + \
                              c2 * c2 - R[k][i] * R[k][i]); 
        R[k-1][i] = c1;
        R[k][i]   = c2;

    
    #omega
    hatr = R[k-1][k-1];
    omega[k-1] = abs(1/hatr);
    for i in xrange(k-1):
        omega[i] = math.sqrt(omega[i]*omega[i] \
                   - u[i]*u[i]/r/r + (u1[i]+mu*u[i])*(u1[i]+mu*u[i])/hatr/hatr );

    #invA_B
    c1         = R[k-1:k, k+1:col];
    c2         = R[k:k+1, k+1:col];
    new_invA_B = copy(invA_B);
    invA_B_row, invA_B_col           = new_invA_B.shape; 
    new_invA_B[k-1][0]               = mu/lou/lou;
    new_invA_B[k-1:k, 1:invA_B_col]  = copy(R[k-1:k, k+1:col])/hatr;
    new_invA_B[0:k-1, 0:1         ]  = (nu*nu*u - mu*u1)/lou/lou;
    new_invA_B[0:k-1, 1:invA_B_col]  = U + (nu*dot(u,c2) - dot(u1,c1)) / hatr; 
    return R, new_invA_B, omega, gamma;    

if __name__ == "__main__":
    Debug = True;
    print "case:***************************************************************************************************************************************************"
    testM = array([[1,0,1,1],[0,1,1,1],[1,1,1,0],[1,0,1,0]]);
    R = rrqr(testM,3);


    print "case:***************************************************************************************************************************************************"
    testM[3,3] = 10;
    R = rrqr(testM,2);


    print "case: **************************************************************************************************************************************************"
    testM = array([[1,0,1,1,1],[0,1,1,1,1],[1,1,1,0,1],[1,0,1,0,1],[10,10,10,10,10]]);
    R = rrqr(testM,3);

    print "case: ***************************************************************************************************************************************************"
    R = rrqr(testM,1);

    for i in xrange(100):
        for k in xrange(1,6): 
            print "case:***************************************************************************************************************************************************";
            del testM;
            testM = np.random.rand(6,6);
            R = rrqr(testM,k);

    
    for i in xrange(100):
        for k in xrange(1,7):
            print "case:***************************************************************************************************************************************************";
            del testM;
            testM = np.random.rand(7,7);
            R = rrqr(testM,k);    
