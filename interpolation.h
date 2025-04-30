#ifndef INTERPOLATION_H
#define INTERPOLATION_H

void calculate_newton_coefficients(int n, double *x, double *f_x, double *c);

double evaluate_newton_polynomial(double t, int n, double *x, double *c);

// Function to build cubic spline interpolation
// Parameters:
// n - number of interpolation points
// x - array of interpolation points (length n)
// f - array of function values at interpolation points (length n)
// a - array to store spline coefficients (length 4*n)
// Second derivatives at boundaries are computed from the function
void build_cubic_spline(int n, double *x, double *f, double *a, double (*func)(double));

// Function to evaluate cubic spline at a given point
// Parameters:
// x_val - point to evaluate
// a - left end of interval
// b - right end of interval
// n - number of interpolation points
// x - array of interpolation points (length n)
// coef - array of coefficients from build_cubic_spline
double eval_cubic_spline(double x_val, double a, double b, int n, double *x, double *coef);

#endif // INTERPOLATION_H 