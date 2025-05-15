#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#include <vector>

// Newton polynomial interpolation
void constructNewtonPolynomial(int n, const double* x, const double* f, double* a, double* work);
double evaluateNewtonPolynomial(double x, double a, double b, int n, const double* points, const double* coeffs);

class AkimaSpline {
public:
    // Construct Akima spline from points and values
    static void constructSpline(int n, const double* x, const double* f, double* coeffs);

    // Evaluate spline at point x
    static double evaluateSpline(double x, double a, double b, int n,
                               const double* points, const double* coeffs);

private:
    // Helper functions for Akima's method
    static double calculateSlope(double x0, double x1, double y0, double y1);
    static double calculateWeight(double s1, double s2);
    static void calculateExtraPoints(const double* x, const double* f, int n,
                                   std::vector<double>& extX, std::vector<double>& extF);

    // Calculate cubic polynomial coefficients for each segment
    static void calculateSegmentCoefficients(double x0, double x1, double y0, double y1,
                                           double m0, double m1, double* coeffs);

    // Find segment containing x
    static int findSegment(double x, const double* points, int n);
};

#endif // APPROXIMATION_H
