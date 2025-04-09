#include <stdio.h>
#include <math.h>

void print_mat(double * mat) {
    printf("[ [ %5lg, %5lg, %5lg ], \n", mat[0], mat[3], mat[6]);
    printf("  [ %5lg, %5lg, %5lg ], \n", mat[1], mat[4], mat[7]);
    printf("  [ %5lg, %5lg, %5lg ] ]\n", mat[2], mat[5], mat[8]);
}


extern double * tot_qp;
extern double * tot_qw;
extern int tot_qn;
extern size_t * ind;
extern size_t el_n;
extern double * x0;

void Interpolate(size_t el, double coord[3], double* x, double v[3]) {
    // $AD II-LOOP
    for (int k=0;k<3;k++) {
        v[k] = (1-coord[2])*((1-coord[0]-coord[1])*x[3*ind[0+6*el]+k]+coord[0]*x[3*ind[1+6*el]+k]+coord[1]*x[3*ind[2+6*el]+k]) + 
                coord[2] *((1-coord[0]-coord[1])*x[3*ind[3+6*el]+k]+coord[0]*x[3*ind[4+6*el]+k]+coord[1]*x[3*ind[5+6*el]+k]);
    }
}

void DP(size_t el, double coord[3], double* x, double M[9]) {
    // $AD II-LOOP
    for (int k=0;k<3;k++) {
        M[k+0*3] = (1-coord[2])*(-x[3*ind[0+6*el]+k]+x[3*ind[1+6*el]+k]) + 
            coord[2]*(-x[3*ind[3+6*el]+k]+x[3*ind[4+6*el]+k]);
        M[k+1*3] = (1-coord[2])*(-x[3*ind[0+6*el]+k]+x[3*ind[2+6*el]+k]) + 
            coord[2]*(-x[3*ind[3+6*el]+k]+x[3*ind[5+6*el]+k]);
        M[k+2*3] = -((1-coord[0]-coord[1])*x[3*ind[0+6*el]+k]+coord[0]*x[3*ind[1+6*el]+k]+coord[1]*x[3*ind[2+6*el]+k]) + 
            ((1-coord[0]-coord[1])*x[3*ind[3+6*el]+k]+coord[0]*x[3*ind[4+6*el]+k]+coord[1]*x[3*ind[5+6*el]+k]);
    }
}

void inv3(double m[9] ,double inv[9], double det_ptr[1]) {
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

void EPS(size_t el, double coord[3], double* x1, double eps[9], double det[1]) {
    double M0[9];
    double M0inv[9];
    double M1[9];
    double M[9];
    DP(el, coord, x0, M0);
    inv3(M0, M0inv, det);
    DP(el, coord, x1, M1);
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

double Energy(size_t el, double coord[3], double coef[3], double* x1, double* v1) {
    double eps[9];
    double det[1];
    double v[3];
    EPS(el, coord, x1, eps, det);
    Interpolate(el, coord, v1, v);
    double a=0,b=0,c=0;
    // $AD II-LOOP
    for (int i=0;i<3;i++) {
        // $AD II-LOOP
        for (int j=0;j<3;j++) a += eps[i+3*j]*eps[i+3*j];
        b += eps[i+i*3]*eps[i+i*3];
        c += v[i]*v[i];
    }
    return (coef[0]*a+coef[1]*b+coef[2]*c)*det[0];
}

double ElementEnergy(size_t el, double coef[3], double* x1, double* v1) {
    double energy = 0;
    // $AD II-LOOP
    for (int i=0; i<tot_qn; i++) {
        double coord[3];
        coord[0] = tot_qp[0+i*3];
        coord[1] = tot_qp[1+i*3];
        coord[2] = tot_qp[2+i*3];
        energy += tot_qw[i] * Energy(el, coord, coef, x1, v1);
    }
    return energy;
}

double TotalEnergy(double coef[3], double* x1, double* v1) {
    double energy = 0;
    // $AD II-LOOP
    for(size_t el=0; el<el_n; el++) {
        energy += ElementEnergy(el, coef, x1, v1);
    }
    return energy;
}

