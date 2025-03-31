#include <stdio.h>
#include <math.h>

void print_mat(double * mat) {
    printf("[ [ %5lg, %5lg, %5lg ], \n", mat[0], mat[3], mat[6]);
    printf("  [ %5lg, %5lg, %5lg ], \n", mat[1], mat[4], mat[7]);
    printf("  [ %5lg, %5lg, %5lg ] ]\n", mat[2], mat[5], mat[8]);
}


double * tot_qp;
double * tot_qw;
int tot_qn;

// (1-coef[2])*((1-coef[0]-coef[1])*x[3*ind[0]+k]+coef[0]*x[3*ind[1]+k]+coef[1]*x[3*ind[2]+k]) + 
//    coef[2] *((1-coef[0]-coef[1])*x[3*ind[3]+k]+coef[0]*x[3*ind[4]+k]+coef[1]*x[3*ind[5]+k]);
void DP(size_t* ind, double* coef, double* x, double* M) {
    for (int k=0;k<3;k++) {
        M[k+0*3] = (1-coef[2])*(-x[3*ind[0]+k]+x[3*ind[1]+k]) + 
            coef[2]*(-x[3*ind[3]+k]+x[3*ind[4]+k]);
        M[k+1*3] = (1-coef[2])*(-x[3*ind[0]+k]+x[3*ind[2]+k]) + 
            coef[2]*(-x[3*ind[3]+k]+x[3*ind[5]+k]);
        M[k+2*3] = -((1-coef[0]-coef[1])*x[3*ind[0]+k]+coef[0]*x[3*ind[1]+k]+coef[1]*x[3*ind[2]+k]) + 
            ((1-coef[0]-coef[1])*x[3*ind[3]+k]+coef[0]*x[3*ind[4]+k]+coef[1]*x[3*ind[5]+k]);
    }
}

void inv3(double* m ,double* inv, double* det_ptr) {
    double det = m[0+3* 0] * (m[1+3* 1] * m[2+3* 2] - m[2+3* 1] * m[1+3* 2]) -
             m[0+3* 1] * (m[1+3* 0] * m[2+3* 2] - m[1+3* 2] * m[2+3* 0]) +
             m[0+3* 2] * (m[1+3* 0] * m[2+3* 1] - m[1+3* 1] * m[2+3* 0]);
    double invdet = 1 / det;
    inv[0+3* 0] = (m[1+3* 1] * m[2+3* 2] - m[2+3* 1] * m[1+3* 2]) * invdet;
    inv[0+3* 1] = (m[0+3* 2] * m[2+3* 1] - m[0+3* 1] * m[2+3* 2]) * invdet;
    inv[0+3* 2] = (m[0+3* 1] * m[1+3* 2] - m[0+3* 2] * m[1+3* 1]) * invdet;
    inv[1+3* 0] = (m[1+3* 2] * m[2+3* 0] - m[1+3* 0] * m[2+3* 2]) * invdet;
    inv[1+3* 1] = (m[0+3* 0] * m[2+3* 2] - m[0+3* 2] * m[2+3* 0]) * invdet;
    inv[1+3* 2] = (m[1+3* 0] * m[0+3* 2] - m[0+3* 0] * m[1+3* 2]) * invdet;
    inv[2+3* 0] = (m[1+3* 0] * m[2+3* 1] - m[2+3* 0] * m[1+3* 1]) * invdet;
    inv[2+3* 1] = (m[2+3* 0] * m[0+3* 1] - m[0+3* 0] * m[2+3* 1]) * invdet;
    inv[2+3* 2] = (m[0+3* 0] * m[1+3* 1] - m[1+3* 0] * m[0+3* 1]) * invdet;
    det_ptr[0] = det;
}

void EPS(size_t* ind, double* coef, double* x0, double* x1, double* eps, double* det_ptr) {
    double M0[9];
    DP(ind, coef, x0, M0);
    double M0inv[9];
    inv3(M0, M0inv, det_ptr);
    double M1[9];
    DP(ind, coef, x1, M1);
    double M[9];
    // $AD II-LOOP
    for (int i=0;i<3;i++) {
        // $AD II-LOOP
        for (int j=0;j<3;j++) {
            M[i+j*3] = 0;
            // $AD II-LOOP
            for (int k=0;k<3;k++) M[i+j*3] += M1[i+k*3] * M0inv[k+j*3];
        }
    }
    // $AD II-LOOP
    for (int i=0;i<3;i++) {
        // $AD II-LOOP
        for (int j=0;j<3;j++) {
            eps[i+j*3] = 0;
            // $AD II-LOOP
            for (int k=0;k<3;k++) eps[i+j*3] += M[i+k*3] * M[j+k*3];
        }
        eps[i+i*3] -= 1;
    }
}

double HookEnergy(double lam, double gam, double* eps) {
    double a=0,b=0;
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) a += eps[i+3*j]*eps[i+3*j];
        b += eps[i+i*3]*eps[i+i*3];
    }
    return lam*a + gam*b;
}

double Energy(double lam, double gam, size_t* ind, double* coef, double* x0, double* x1) {
    double eps[9];
    double det;
    EPS(ind, coef, x0, x1, eps, &det);
    return HookEnergy(lam, gam, eps)*det;
}

double ElementEnergy(double lam, double gam, size_t* ind, double* x0, double* x1) {
    double energy = 0;
    // $AD II-LOOP
    for (int i=0; i<tot_qn; i++) {
        double coef[3];
        coef[0] = tot_qp[0+i*3];
        coef[1] = tot_qp[1+i*3];
        coef[2] = tot_qp[2+i*3];
        energy += tot_qw[i] * Energy(lam, gam, ind, coef, x0, x1);
    }
    return energy;
}

double TotalEnergy(double lam, double gam, size_t* ind, size_t ind_n, double* x0, double* x1) {
    double energy = 0;
    // $AD II-LOOP
    for(size_t i=0; i<ind_n; i++) {
        size_t el_ind[6];
        for (int j=0; j<6; j++) el_ind[j] = ind[j+6*i];
        energy += ElementEnergy(lam, gam, el_ind, x0, x1);
    }
    return energy;
}

