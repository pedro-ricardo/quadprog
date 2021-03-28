#include <stdio.h>
#include <math.h>
#include "../src/quad_prog.h"

int main() {
    // Sizes
    int n = 3;
    int q = 2;
    int qeq = 2;

    double tol = 1.0e-3;
    double res = 0;
    int err = 0;
    double res_c = -1;

    // Matrix Set
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
    double beq[] = {0, 21800};

    double lb[] = {1.944, 0, 1.944};
    double ub[] = {73.16, 49.30, 33.33};
    double x[] = {0,0,0};

    double *plb = lb;
    double *pub = ub;

    // Calculate
    quadprog(n, q, qeq, H, f, Aeq, beq, &plb, &pub, x, &tol, &res, &err);
    res_c = qp_res_correction(n,f,&res);
    
    // Print Results
    printf("\n Calculation Done\n");
    printf(" Error code :%2d\n",err);
    printf(" Eq. Residue:%10.3e\n",res_c);
    
    printf("\n Solution:\n");
    for(int i = 0; i < n; i++) {
        printf(" %8.3f\n", x[i]);
    }

    return 0;
}
