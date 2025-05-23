#include "approximation.h"
#include <cmath>
#include <algorithm>

// Newton polynomial interpolation implementation
void constructNewtonPolynomial(int n, const double* x, const double* f, double* a, double* work) {
    // work array is used to store divided differences
    // First column of work array stores function values
    for (int i = 0; i < n; ++i) {
        work[i] = f[i];
    }

    // Calculate divided differences
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            work[i + j * n] = (work[i + 1 + (j - 1) * n] - work[i + (j - 1) * n]) / (x[i + j] - x[i]);
        }
    }

    // Store coefficients in a
    for (int i = 0; i < n; ++i) {
        a[i] = work[i * n];
    }
}

double evaluateNewtonPolynomial(double x, double a, double b, int n, const double* points, const double* coeffs) {
    // Validate input range
    if (x < a || x > b) {
        return 0.0;  // Return 0 for points outside the interpolation range
    }

    double result = coeffs[0];
    double product = 1.0;

    for (int i = 1; i < n; ++i) {
        product *= (x - points[i - 1]);
        result += coeffs[i] * product;
    }

    return result;
}

void AkimaSpline::constructSpline(int n, const double* x, const double* f, double* coeffs) {
    if (n < 2) return;

    // Create extended points array with extra points for boundary conditions
    std::vector<double> extX, extF;
    calculateExtraPoints(x, f, n, extX, extF);

    // Calculate slopes between points
    std::vector<double> slopes(n + 3);
    for (int i = 0; i < n + 3; ++i) {
        slopes[i] = calculateSlope(extX[i], extX[i+1], extF[i], extF[i+1]);
    }

    // Calculate weights and tangents at each point
    std::vector<double> tangents(n);
    for (int i = 0; i < n; ++i) {
        double w1 = calculateWeight(slopes[i], slopes[i+1]);
        double w2 = calculateWeight(slopes[i+2], slopes[i+3]);

        // Improve numerical stability with a more robust comparison
        if (std::abs(w1 + w2) < 1e-10) {
            // If weights are near zero, use the average of slopes as fallback
            tangents[i] = 0.5 * (slopes[i+1] + slopes[i+2]);
        } else {
            tangents[i] = (w1 * slopes[i+2] + w2 * slopes[i+1]) / (w1 + w2);
        }
    }

    // Calculate coefficients for each segment
    for (int i = 0; i < n - 1; ++i) {
        calculateSegmentCoefficients(x[i], x[i+1], f[i], f[i+1],
                                   tangents[i], tangents[i+1],
                                   coeffs + i * 4);
    }
}

double AkimaSpline::evaluateSpline(double x, double a, double b, int n,
                                 const double* points, const double* coeffs) {
    // Validate input range
    if (x < a || x > b) {
        return 0.0;  // Return 0 for points outside the interpolation range
    }
    // Find the segment containing x
    int segment = findSegment(x, points, n);
    if (segment < 0) return 0.0;

    // Get coefficients for this segment
    const double* segCoeffs = coeffs + segment * 4;

    // Calculate relative x
    double dx = x - points[segment];

    // Evaluate cubic polynomial
    return segCoeffs[0] + dx * (segCoeffs[1] + dx * (segCoeffs[2] + dx * segCoeffs[3]));
}

double AkimaSpline::calculateSlope(double x0, double x1, double y0, double y1) {
    double dx = x1 - x0;
    if (std::abs(dx) < 1e-10) return 0.0;
    return (y1 - y0) / dx;
}

double AkimaSpline::calculateWeight(double s1, double s2) {
    double diff = std::abs(s2 - s1);
    if (diff < 1e-10) return 1.0;
    return 1.0 / diff;
}

void AkimaSpline::calculateExtraPoints(const double* x, const double* f, int n,
                                     std::vector<double>& extX, std::vector<double>& extF) {
    extX.resize(n + 4);
    extF.resize(n + 4);

    // Copy original points
    std::copy(x, x + n, extX.begin() + 2);
    std::copy(f, f + n, extF.begin() + 2);

    // Add extra points at boundaries
    double h = x[1] - x[0];
    extX[1] = x[0];
    extX[0] = x[0] - h;
    extF[1] = f[0];
    extF[0] = f[0] - (f[1] - f[0]);

    h = x[n-1] - x[n-2];
    extX[n+2] = x[n-1];
    extX[n+3] = x[n-1] + h;
    extF[n+2] = f[n-1];
    extF[n+3] = f[n-1] + (f[n-1] - f[n-2]);
}

void AkimaSpline::calculateSegmentCoefficients(double x0, double x1, double y0, double y1,
                                             double m0, double m1, double* coeffs) {
    double dx = x1 - x0;
    double dy = y1 - y0;

    coeffs[0] = y0;  // a
    coeffs[1] = m0;  // b

    // Calculate c and d coefficients
    double t = (dy / dx - m0) / dx;
    coeffs[2] = (3*t - (m1 - m0)/(dx)) / dx;  // c
    coeffs[3] = ((m1 - m0)/(dx) - 2*t) / dx;  // d
}

int AkimaSpline::findSegment(double x, const double* points, int n) {
    if (x < points[0] || x > points[n-1]) return -1;

    // Binary search for the segment
    int left = 0, right = n - 2;
    while (left <= right) {
        int mid = (left + right) / 2;
        if (x >= points[mid] && x <= points[mid+1]) return mid;
        if (x < points[mid]) right = mid - 1;
        else left = mid + 1;
    }

    return n - 2;  // Last segment
}
