#include <stdio.h>
#include <math.h>

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

#define n 5

void f (double a, double t, double *x, double *s);
void f (double a, double t, double *x, double *s) { // s --- вектор из производных
    
    s[0] = x[1];
    s[1] = x[3] - x[0] * exp(-a * t);
    s[2] = x[3] * exp(-a * t);
    s[3] = -x[2];
    s[4] = x[3] * x[3];
}

double lambda(double a, double t);
double lambda(double a, double t) {
   
    double l;
    l = exp(-3*a*t)*(exp(3*a*t)+exp(2*a*t)*sqrt(5*exp(2*a*t)-8*exp(a*t)+4))/4;
    return l;
}

double hnew (double eps, double err, double h);
double hnew (double eps, double err, double h) {
    
    double h_new;
    double chi;
    
    if (fabs(err) < 1e-15)
        chi = 10;
    else {
        chi = pow(eps / err, 1.0 / 6);
        if (chi > 10) chi = 10;
        if (chi < 0.1) chi = 0.1;
    }
    h_new = 0.95 * h * chi;
    
    return h_new;
}

double norma (double a0, double b0, double s0, double s1);
double norma (double a0, double b0, double s0, double s1) {
    
    double norma;
    norma = sqrt((s0-a0)*(s0-a0) + (s1-b0)*(s1-b0));
    return norma;
}

double Runge_Kutta(double a, double t0, void(*f)(double, double, double*, double*), double x[n], double x1[n], double h);
double Runge_Kutta(double a, double t0, void(*f)(double, double, double*, double*), double x[n], double x1[n], double h) {
    
    double x_tmp[n];
    double k1[n], k2[n], k3[n], k4[n], k5[n], k6[n], k7[n];
    double err = 0;
    
    static const double r_45 = 1.0 / 45.0;
    static const double r_8_9 = 8.0 / 9.0;
    static const double r_6561 = 1.0 / 6561.0;
    static const double r_167904 = 1.0 / 167904.0;
    static const double r_142464 = 1.0 / 142464.0;
    static const double r_21369600 = 1.0 / 21369600.0;
    
    for (int i = 0; i < n; i++)
        x_tmp[i] = x[i];
    (*f)(a, t0, x_tmp, k1);
    for (int i = 0; i < n; i++)
        x_tmp[i] = x[i] + 0.2 * h * k1[i];
    (*f)(a, t0 + 0.2 * h, x_tmp, k2);
    for (int i = 0; i < n; i++)
        x_tmp[i] = x[i] + h * (0.075 * k1[i] + 0.225 * k2[i]);
    (*f)(a, t0 + 0.3 * h, x_tmp, k3);
    for (int i = 0; i < n; i++)
        x_tmp[i] = x[i] + h * r_45 * (44.0 * k1[i] - 168.0 * k2[i] + 160 * k3[i]);
    (*f)(a, t0 + 0.8 * h, x_tmp, k4);
    for (int i = 0; i < n; i++)
        x_tmp[i] = x[i] + r_6561 * h * (19372.0 * k1[i] - 76080.0 * k2[i] + 64448.0 * k3[i] - 1908.0 * k4[i]);
    (*f)(a, t0 + r_8_9 * h, x_tmp, k5);
    for (int i = 0; i < n; i++)
        x_tmp[i] = x[i] + r_167904 * h * (477901.0 * k1[i] - 1806240.0 * k2[i]
        + 1495424.0 * k3[i] + 46746.0 * k4[i] - 45927.0 * k5[i]);
    (*f)(a, t0 + h, x_tmp, k6);
    for (int i = 0; i < n; i++)
        x_tmp[i] = x[i] + r_142464 * h * (12985.0 * k1[i] + 64000.0 * k3[i] + 92750.0 * k4[i] - 45927.0 * k5[i] + 18656.0 * k6[i]);
    (*f)(a, t0 + h, x_tmp, k7);
    for (int i = 0; i < n; i++)
    {
        x1[i] = x[i] + r_21369600 * h * (1921409.0 * k1[i] + 9690880.0 * k3[i] + 13122270.0 * k4[i] - 5802111.0 * k5[i] + 1902912.0 * k6[i] + 534240.0 * k7[i]);
        err += fabs(r_21369600 * (26341.0 * k1[i] - 90880.0 * k3[i] + 790230.0 * k4[i] - 1086939.0 * k5[i] + 895488.0 * k6[i] - 534240.0 * k7[i]));
    }
    return err;
}

void GlobalError (double a, void(*f)(double, double, double*,double*), double x0[n], double x1[n], double h, double t_max, double eps, double *delta);
void GlobalError (double a, void(*f)(double, double, double*,double*), double x0[n], double x1[n], double h, double t_max, double eps, double *delta) {
    
    double err;
    double delta1;
    double t;
    double h_new;
    double x_prev[n];
    double x_tmp[n];
    int at = 0;
    
    h_new = h;
    t = 0;
    delta1 = 0;
    for (int i = 0; i < n; i ++) {
        x_tmp[i] = x0[i];
    }
    while (t < t_max) {
        for (int i = 0; i < n; i ++)
            x_prev[i] = x_tmp[i];
        err = Runge_Kutta(a, t, f, x_prev, x_tmp, h_new);
        while (err > eps) {
            h_new = min(hnew(eps, err, h_new), t_max - t);
            err = Runge_Kutta(a, t, f, x_prev, x_tmp, h_new);
        }
        delta1 = err + delta1 * exp(h_new * (lambda(a, t+h_new) - lambda(a, t))); // интеграл вычисляем по формуле левых прямоугольников
        t += h_new;
        h_new = min(hnew(eps, err, h_new), t_max - t);
        at++;
    }
    *delta = delta1;
    for (int i = 0; i < n; i ++)
        x1[i] = x_tmp[i];
}

void Inverse_matr (double res[2], double a0, double b0, double eps, double eps0, double h, double tmax, double a); //обращает матрицу
void Inverse_matr (double res[2], double a0, double b0, double eps, double eps0, double h, double tmax, double a) {
    
    double x[2][2];
    double x1[n];
    double x0[n];
    double x_tmp1[n];
    double x_tmp2[n];
    double delta;
    double tmp;
    
    x0[0] = a0;
    x0[1] = 0;
    x0[2] = 0;
    x0[3] = b0;
    x0[4] = 0;
    
    GlobalError(a,f, x0, x1, h, tmax, eps, &delta);
    
    x0[0] = a0+eps0;
    x0[1] = 0;
    x0[2] = 0;
    x0[3] = b0;
    x0[4] = 0;
    
    GlobalError(a,f, x0, x_tmp1, h, tmax, eps, &delta);
    
    x0[0] = a0;
    x0[1] = 0;
    x0[2] = 0;
    x0[3] = b0+eps0;
    x0[4] = 0;
    
    GlobalError(a,f, x0, x_tmp2, h, tmax, eps, &delta);
    
    tmp = (x_tmp1[0] - x1[0])*(x_tmp2[1] - x1[1])/(eps0*eps0) - (x_tmp1[1] - x1[1])*(x_tmp2[0] - x1[0])/(eps0*eps0);
    
    x[0][0] = (x_tmp2[1] - x1[1])/(eps0 * tmp);
    x[0][1] = -(x_tmp2[0] - x1[0])/(eps0 * tmp);
    x[1][0] = -(x_tmp1[1] - x1[1])/(eps0 * tmp);
    x[1][1] = (x_tmp1[0] - x1[0])/(eps0 * tmp);
    
    res[0] = x[0][0]*x1[0] + x[0][1]*(x1[1] + M_PI_2);
    res[1] = x[1][0]*x1[0] + x[1][1]*(x1[1] + M_PI_2);
    
}

void Newton (double *a00, double *b00, double eps, double eps0, double h, double tmax, double a);
void Newton (double *a00, double *b00, double eps, double eps0, double h, double tmax, double a) {
    
    double s0, s1;
    double eps1 = 1e-11;
    double x[2];
    int at = 0;
    double a0, b0;
    a0 = *a00;
    b0 = *b00;
    do {
        s0 = a0;
        s1 = b0;
        Inverse_matr(x, s0, s1, eps, eps0, h, tmax, a);
        a0 = s0 - x[0];
        b0 = s1 - x[1];
        at++;
    } while (norma(a0,b0,s0,s1) > eps1);
    *a00 = a0;
    *b00 = b0;
}

int main() {
    double h = 0.001;
    double a0 = 1.3;
    double b0 = -1;
    double tmax = M_PI_2;
    double a = 10;
    double x0[n];
    double x1[n];
    double x2[n];
    double x3[n];
    double x4[n];
    double delta = 0;
    
    Newton (&a0, &b0, 1e-11, 1e-11, h, tmax, a);
    printf("a0 = %.12lf\n", a0);
    printf("b0 = %.12lf\n", b0);
    
    x0[0] = a0;
    x0[1] = 0;
    x0[2] = 0;
    x0[3] = b0;
    x0[4] = 0;
    GlobalError(a,f, x0, x1, h, tmax, 1e-7, &delta);
    GlobalError(a,f, x0, x2, h, tmax, 1e-9, &delta);
    printf("delta = %.12lf\n", delta);
    printf("7-9 = %.12lf\n", x1[4]-x2[4]);
    GlobalError(a,f, x0, x3, h, tmax, 1e-9, &delta);
    GlobalError(a,f, x0, x4, h, tmax, 1e-11, &delta);
    printf("9-11 = %.12lf\n", x3[4]-x4[4]);
    return 0;
}
