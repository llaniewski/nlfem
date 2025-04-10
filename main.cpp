#include <stdio.h>
#include <math.h>

double * tot_qp;
double * tot_qw;
int tot_qn;
size_t * ind;
size_t el_n;
double * x0;
double lam;
double gam;
double rho;

extern "C" {
    double TotalEnergy(double coef[3], double* x1, double* v1);
}

int main () {
    int tri_qn = 6;
    double tri_qw[] = { 0.223381589678011,
                        0.223381589678011,
                        0.223381589678011,
                        0.109951743655322,
                        0.109951743655322,
                        0.109951743655322};
    double tri_qp1[]= {	0.445948490915965,
                        0.445948490915965,
                        0.108103018168070,
                        0.091576213509771,
                        0.091576213509771,
                        0.816847572980459};
    double tri_qp2[]= {	0.108103018168070,
                        0.445948490915965,
                        0.445948490915965,
                        0.816847572980459,
                        0.091576213509771,
                        0.091576213509771};
    int int_qn = 3;
    double int_qp[] = { -sqrt(3.0/5.0),
                        0.0,
                        sqrt(3.0/5.0)};
    double int_qw[] = { 5.0/18.0,
                        8.0/18.0,
                        5.0/18.0 };

    tot_qn = tri_qn * int_qn;
    tot_qp = new double[3*tot_qn];
    tot_qw = new double[tot_qn];
    for (int i=0;i<tri_qn;i++){
        for (int j=0;j<int_qn;j++){
            tot_qp[0+3*(i+tri_qn*j)] = tri_qp1[i];
            tot_qp[1+3*(i+tri_qn*j)] = tri_qp2[i];
            tot_qp[2+3*(i+tri_qn*j)] = int_qp[j];
            tot_qw[i+tri_qn*j] = tri_qw[i] * int_qw[j]/2;
        }
    }


    int mx = 10;
    int my = 10;
    int pnt_n = mx*my*2;
    el_n = (mx-1)*(my-1)*2;
    x0 = new double[pnt_n*3];
    double * x = new double[pnt_n*3];
    double * v = new double[pnt_n*3];
    ind = new size_t[el_n*6];
    double LX = 1;
    double LY = 1;
    double LZ = 1;
    for (int i=0;i<mx;i++){
        for (int j=0;j<my;j++){
            for (int k=0;k<2;k++){
                x0[0+3*(i+mx*(j+my*k))] = LX*i/(mx-1);
                x0[1+3*(i+mx*(j+my*k))] = LY*j/(my-1);
                x0[2+3*(i+mx*(j+my*k))] = LZ*k;
            }
        }
    }
    for (int i=0;i<3*mx*my*2;i++) x[i] = x0[i];
    for (int i=0;i<3*mx*my*2;i++) v[i] = 1.0;
    double a = 0;
    for (int i=0;i<pnt_n;i++) {
        x[0+3*i] = x0[0+3*i]*cos(a) - x0[1+3*i]*sin(a);
        x[1+3*i] = x0[0+3*i]*sin(a) + x0[1+3*i]*cos(a);
        x[2+3*i] = x0[2+3*i]*1.0;
    }
    for (int i=0;i<mx-1;i++){
        for (int j=0;j<my-1;j++){
            ind[0+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*(j+my*0);
            ind[1+6*(i+(mx-1)*(j+(my-1)*0))] = (i+1)+mx*(j+my*0);
            ind[2+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*((j+1)+my*0);
            ind[3+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*(j+my*1);
            ind[4+6*(i+(mx-1)*(j+(my-1)*0))] = (i+1)+mx*(j+my*1);
            ind[5+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*((j+1)+my*1);
            ind[0+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*((j+1)+my*0);
            ind[1+6*(i+(mx-1)*(j+(my-1)*1))] = i+mx*((j+1)+my*0);
            ind[2+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*(j+my*0);
            ind[3+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*((j+1)+my*1);
            ind[4+6*(i+(mx-1)*(j+(my-1)*1))] = i+mx*((j+1)+my*1);
            ind[5+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*(j+my*1);
        }
    }
    
    lam = 1;
    gam = 1;
    rho = 1;
    double coef[3] = {lam,gam,rho};
    double energy = TotalEnergy(coef, x,v);
    printf("%lg\n",energy);

    return 0;
}