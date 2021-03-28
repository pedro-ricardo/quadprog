// This Header file contains the interface to use the quadprog function
// in C and C++ language
#ifdef __cplusplus
extern "C" {
#endif

    // Quadratic Programming as in MATLAB sintax.
    void quadprog (const int n, const int q, const int qeq,
        const double *vH ,const double *f, const double *vAeq, const double *beq, 
        double **lb, double **ub, double *x,
        const double *tol, double *res, int *err );

    // True residual from equation.
    double qp_res_correction(const int n, const double *f, double *res);

#ifdef __cplusplus
}
#endif
