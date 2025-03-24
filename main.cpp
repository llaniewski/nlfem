#include <stdio.h>

void DP(size_t* ind, size_t el, double* coef, double* x, double* M) {
    for (int k=0;k<3;k++) {
        M[k+0*3] = (1-coef[2])*(-x[3*ind[6*el+0]+k]+x[3*ind[6*el+1]+k]) + 
            coef[2]*(-x[3*ind[6*el+3]+k]+x[3*ind[6*el+4]+k]);
        M[k+1*3] = (1-coef[2])*(-x[3*ind[6*el+0]+k]+x[3*ind[6*el+2]+k]) + 
            coef[2]*(-x[3*ind[6*el+3]+k]+x[3*ind[6*el+5]+k]);
        M[k+2*3] = -((1-coef[0]-coef[1])*x[3*ind[6*el+0]+k]+coef[0]*x[3*ind[6*el+1]+k]+coef[1]*x[3*ind[6*el+2]+k]) + 
            ((1-coef[0]-coef[1])*x[3*ind[6*el+3]+k]+coef[0]*x[3*ind[6*el+4]+k]+coef[1]*x[3*ind[6*el+5]+k]);
    }
}

void inv3(double* m ,double* inv) {
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
}

void EPS(size_t* ind, size_t el, double* coef, double* x0, double* x1, double* eps) {
    double M0[9];
    DP(ind, el,coef, x0, M0);
    double M0inv[9];
    inv3(M0,M0inv);
    double M1[9];
    DP(ind, el,coef, x1, M1);
    double M[9];
    for (int i=0;i<3;i++) {
        
        for (int j=0;j<3;j++) {
            M[i+j*3] = 0;
            for (int k=0;k<3;k++) M[i+j*3] += M0inv[i+k*3]*M1[k+j*3];
        }
    }
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            eps[i+j*3] = 0;
            for (int k=0;k<3;k++) eps[i+j*3] += M[i+k*3]*M[j+k*3];
        }
        eps[i+i*3] -= 1;
    }
}

void Hook(double lam, double gam, double* eps, double* sigma) {
    double tr = 0;
    for (int k=0;k<3;k++) tr += eps[k+3*k];
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) sigma[i+j*3] = lam*eps[i+3*j];
        sigma[i+i*3] += gam;
    }
}

double HookEnergy(double lam, double gam, double* eps) {
    double a=0,b=0;
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) a += eps[i+3*j];
        b += eps[i+i*3];
    }
    return lam*a + gam*b;
}

double Energy(double lam, double gam, size_t* ind, size_t el, double* coef, double* x0, double* x1) {
    double eps[9];
    EPS(ind, el, coef, x0, x1, eps);
    return HookEnergy(lam, gam, eps);
}

double ElementEnergy(double lam, double gam, size_t* ind, size_t el, double* x0, double* x1) {
    double energy = 0;
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) if (i+j<2) {
            for(int k=0; k<2; k++) {
                double coef[3];
                coef[0] = i;
                coef[1] = j;
                coef[2] = k;
                energy += Energy(lam, gam, ind, el, coef, x0, x1);
            }
        }
    }
    return energy;
}

double TotalEnergy(size_t *ind, size_t el_n, double lam, double gam, double* x0, double* x1) {
    double energy = 0;
    for(size_t el=0; el<el_n; el++) {
        energy += ElementEnergy(lam, gam, ind, el, x0, x1);
    }
    return energy;
}



int main () {

    return 0;
}