#include <stdio.h>
#include <math.h>
#include <functional>
#include <memory>
#include <string>

double * tot_qp;
double * tot_qw;
int tot_qn;
size_t * ind;
size_t el_n;
double * x0;
double lam;
double gam;
double rho;

typedef std::vector<double> vec;

extern "C" {
    double TotalEnergy(double coef[3], const double* x1, const double* v1);
    void TotalEnergy_b(double coef[3], const double *x1, double *x1b, const double *v1, double *v1b, double TotalEnergyb);
    void TotalEnergy_b_d(double coef[3], const double *x1, const double *x1d, double *x1b, 
        double *x1bd, const double *v1, const double *v1d, double *v1b, double *v1bd, 
        double TotalEnergyb);
}

void write_vtu(char* filename, vec& points, vec& v, std::vector<size_t>& triangles);

double skal(const vec& a, const vec& b) {
    size_t n = a.size();
    double sum = 0;
    for (int i=0;i<n;i++) sum += a[i]*b[i];
    return sum;
}

int Solve_(std::function<void(const vec&, vec&)> mult, const vec& b, vec& x, int max_iter=100, double eps=1e-6) {
    size_t n = b.size();
    static vec r; r.resize(n);
    static vec Ar; Ar.resize(n);
    static vec p; p.resize(n);
    static vec Ap; Ap.resize(n);
    x.resize(n);
    for (int iter=0; iter<max_iter; iter++) {
        mult(x,r);
        for (size_t i=0;i<n;i++) r[i] = b[i] - r[i];
        double res = sqrt(skal(r,r));
        printf("linear %5d res=%lg\n",iter, res);
        if (res < eps) return iter;
        if (iter>0) {
            //mult(r,Ar);
            double beta = skal(Ap,r) / skal(Ap,p);
            //double beta = 0;
            for (size_t i=0;i<n;i++) p[i] = r[i] - beta*p[i];
        } else {
            for (size_t i=0;i<n;i++) p[i] = r[i];
        }
        mult(p,Ap);
        double alpha = skal(p,r) / skal(Ap,p);
        for(size_t i=0;i<n;i++) x[i] = x[i] + alpha*p[i];
//        for(size_t i=0;i<n;i++) r[i] = r[i] - alpha*Ap[i];
    }
    return max_iter;
}



int Solve(std::function<void(const vec&, vec&)> mult, const vec& b, vec& x, int max_iter=100, double eps=1e-6) {
    size_t n = b.size();
    static vec r; r.resize(n);
    static vec p; p.resize(n);
    static vec Ap; Ap.resize(n);
    static vec q; q.resize(n);
    static vec Aq; Aq.resize(n);
    double res;
    for (int iter = 0; iter < max_iter; iter++) {
        mult(x, r);
        for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
        res = sqrt(skal(r, r));
        //res_draw(res);
        
        if (res < eps) {
            printf("linear %5d iterations converged (%lg)\n", iter, res);
            return iter;
        }
        for (int i = 0; i < n; i++) {
            p[i] = r[i];
        }
        if (iter > 0) {
            mult(p, Ap);
            mult(q, Aq);
            double beta = skal(Aq, p) / skal(Aq, q);
            for (int i = 0; i < n; i++) p[i] = p[i] - q[i] * beta;
        }
        mult(p, Ap);
        //printf("|p|^2=%lg\n",skal(p,p));
        //printf("|Ap|^2=%lg\n",skal(Ap,Ap));
        double alpha = skal(p, r) / skal(Ap, p);
        //printf("alpha=%lg\n",alpha);
        for (int i = 0; i < n; i++) x[i] = x[i] + p[i]*alpha;
        for (int i = 0; i < n; i++) q[i] = p[i];
    }
    printf("linear unconverged final residual=%lg\n", res);
    return max_iter;
}


void GetMatrix(size_t n, std::function<void(const vec&, vec&)> mult, vec& M) {
    static vec x; x.resize(n);
    static vec r; r.resize(n);
    M.resize(n*n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) x[j] = 0;
        x[i] = 1;
        mult(x, r);
        for (int j = 0; j < n; j++) M[i+j*n] = r[j];
    }
}

void PrintMatrix(size_t n, vec& M) {
    M.resize(n*n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a = M[i+j*n];
            if (fabs(a) < 1e-10) a = 0;
            printf("%3.3lg", a);
            if (j<n-1) printf(" ");
        }
        printf("\n");
    }
}


void ODE(double dt, vec& x1, vec& v1, vec& x2, vec& v2) {
    size_t n = x1.size();
    static vec Lx1;
    static vec Lv1;
    static vec Lx2;
    static vec Lv2;
    static vec Lx2d;
    static vec Lv2d;
    static vec x2d;
    static vec r;
    static vec dr;
    Lx1.resize(n);
    Lv1.resize(n);
    Lx2.resize(n);
    Lv2.resize(n);
    Lx2d.resize(n);
    Lv2d.resize(n);
    x2d.resize(n);
    r.resize(n);
    dr.resize(n);
    
    x2.resize(n);
    v2.resize(n);

    double coef[3] = {lam,gam,-rho};
    TotalEnergy_b(coef, x1.data(), Lx1.data(), v1.data(), Lv1.data(), 1);
    for (int i=0;i<n;i++) x2[i] = x1[i];
    for (int i=0;i<n;i++) v2[i] = v1[i];
    for (int titer=0; titer<30; titer++) {
        for (int iter=0; iter<30; iter++) {
            for (int i=0;i<n;i++) Lv2[i] = 0;
            for (int i=0;i<n;i++) Lx2[i] = 0;
            for (int i=0;i<n;i++) x2[i] = x1[i] + (v1[i] + v2[i])/2*dt;
            TotalEnergy_b(coef, x2.data(), Lx2.data(), v2.data(), Lv2.data(), 1);
            for (int i=0;i<n;i++) r[i] = (Lv2[i]-Lv1[i])/dt - (Lx1[i]+Lv1[i])/2;
            double res = sqrt(skal(r,r));
            printf("nonlin %5d res=%lg\n",iter, res);
            std::function<void(const vec&, vec&)> mult = [&coef, &x2, &v2, &dt, &n](const vec& v2d, vec& rd){
                for (int i = 0; i < n; ++i) Lx2d[i] = 0;
                for (int i = 0; i < n; ++i) Lv2d[i] = 0;
                // for (int i = 0; i < n; ++i) x2d[i] = 0;
                for (int i = 0; i < n; ++i) x2d[i] = dt*v2d[i]/2;
                TotalEnergy_b_d(coef, x2.data(), x2d.data(), Lx2.data(), Lx2d.data(), v2.data(), v2d.data(), Lv2.data(), Lv2d.data(), 1.0);
                // for (int i = 0; i < n; ++i) rd[i] = Lv2d[i];
                for (int i = 0; i < n; ++i) rd[i] = (0-Lv2d[i])/dt - (0+Lx2d[i])/2;
                
                // // for (int i = 0; i < n; ++i) rd[i] = 0;
                // for (int i = 0; i < n-1; ++i) {
                //     double a = v2d[i+1] - v2d[i];
                //     rd[i] += -a;
                //     rd[i+1] += a;
                // }
                // // rd[0] = v2d[0];
                // rd[n-1] = v2d[n-1];
                //for (int i = 0; i < n; ++i) rd[i] = v2d[i];
            };
            //vec M;
            //GetMatrix(n,mult,M);
            //PrintMatrix(n,M);
            Solve(mult, r, dr);
            for (int i = 0; i < n; ++i) v2[i] = v2[i] + dr[i];
        }
        for (int i=0;i<n;i++) x1[i] = x2[i];
        for (int i=0;i<n;i++) v1[i] = v2[i];    
    }
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


    std::string outpath = "output/";
    std::string name = "box";

    int mx = 3;
    int my = 3;
    int pnt_n = mx*my*2;
    el_n = (mx-1)*(my-1)*2;

    vec points;
    points.resize(pnt_n*3);
    double LX = 1;
    double LY = 1;
    double LZ = 1;
    for (int i=0;i<mx;i++){
        for (int j=0;j<my;j++){
            for (int k=0;k<2;k++){
                points[0+3*(i+mx*(j+my*k))] = LX*i/(mx-1);
                points[1+3*(i+mx*(j+my*k))] = LY*j/(my-1);
                points[2+3*(i+mx*(j+my*k))] = LZ*k;
            }
        }
    }

    std::vector<size_t> triangles;
    triangles.resize(el_n*6);
    for (int i=0;i<mx-1;i++){
        for (int j=0;j<my-1;j++){
            triangles[0+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*(j+my*0);
            triangles[1+6*(i+(mx-1)*(j+(my-1)*0))] = (i+1)+mx*(j+my*0);
            triangles[2+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*((j+1)+my*0);
            triangles[3+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*(j+my*1);
            triangles[4+6*(i+(mx-1)*(j+(my-1)*0))] = (i+1)+mx*(j+my*1);
            triangles[5+6*(i+(mx-1)*(j+(my-1)*0))] = i+mx*((j+1)+my*1);
            triangles[0+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*((j+1)+my*0);
            triangles[1+6*(i+(mx-1)*(j+(my-1)*1))] = i+mx*((j+1)+my*0);
            triangles[2+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*(j+my*0);
            triangles[3+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*((j+1)+my*1);
            triangles[4+6*(i+(mx-1)*(j+(my-1)*1))] = i+mx*((j+1)+my*1);
            triangles[5+6*(i+(mx-1)*(j+(my-1)*1))] = (i+1)+mx*(j+my*1);
        }
    }
    


    x0 = points.data();
    ind = triangles.data();

    vec x(pnt_n*3);
    vec v(pnt_n*3);
    for (int i=0;i<3*mx*my*2;i++) x[i] = points[i];
    for (int i=0;i<3*mx*my*2;i++) v[i] = 1.0;

    double a = 0;
    for (int i=0;i<pnt_n;i++) {
        x[0+3*i] = points[0+3*i]*cos(a) - points[1+3*i]*sin(a);
        x[1+3*i] = points[0+3*i]*sin(a) + points[1+3*i]*cos(a);
        x[2+3*i] = points[2+3*i]*1.0;
    }
    
    lam = 1;
    gam = 1;
    rho = 1;
    double coef[3] = {lam,gam,-rho};
    double energy = TotalEnergy(coef, x.data(),v.data());
    printf("E = %lg\n",energy);

    vec nx(pnt_n*3);
    vec nv(pnt_n*3);
    double dt = 0.01;
    int iter = 0;
    for (double t=0; t<0.1;t+=dt) {
        ODE(dt, x, v, nx, nv);
        char str[1024];
        sprintf(str, "%s%s_%08d.vtu", outpath.c_str(), name.c_str(), iter);
        write_vtu(str, x, v, triangles);
        for (int i=0;i<pnt_n*3;i++) {
            x[i] = nx[i];
            v[i] = nv[i];
        }
        iter++;
    }
    // size_t n = x.size();
    // vec b; b.resize(n);
    // for (int i = 0; i < n; ++i) {
    //     x[i] = 0;
    //     b[i] = 1;
    // }

    // Solve([](const vec& x, vec& r){
    //     size_t n = x.size();
    //     for (int i = 0; i < n; ++i) r[i] = 0;
    //     for (int i = 0; i < n-1; ++i) {
    //         double a = x[i+1] - x[i];
    //         r[i] += -a;
    //         r[i+1] += a;
    //     }
    //     r[0] = x[0];
    //     r[n-1] = x[n-1];
    // }, b, x);


    return 0;
}