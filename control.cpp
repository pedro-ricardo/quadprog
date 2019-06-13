#include <stdio.h>

extern "C" {
    void quadprog (const int *n, const int *qeq, const double *vH ,const double *f, const double *vAeq, const double *beq, 
    const double *lb, const double *ub, double *x, const double *tol, double *res, int *err ) ;
    
    double qp_res_correction(const int *n, const double *f, double *res);
}
//g++ control.cpp quad_prog.o solve_qp.o util.o -lgfortran -lblas
main()
{
    int n = 3;
    int qeq = 2;
    double H[] = {
        2, 0, 0,
        0, 0, 0,
        0, 0, 0
    };
    double f[] = {-68.72, 0, 0};
    double Aeq[] = {
        1       , -1       , -1,
        3.4144e3, -3.0951e3, -2.7316e3
    };
    double beq[] = {0, 21800  
    };
    double lb[] = {1.944, 0, 1.944};
    double ub[] = {73.16, 49.30, 33.33};
    double x[] = {0,0,0};
    double tol = 5;
    double res = 0;
    int err = 0;
    double res_c = -1;
    
    quadprog(&n, &qeq, H, f, Aeq, beq, lb, ub, x, &tol, &res, &err);
    
    
    for(int i = 0; i < n; i++)
    {
        printf("%f\n", x[i]);
    }
    
    res_c = qp_res_correction(&n,f,&res);
    
    printf("%e\n",res_c);
    
    return 0;
}