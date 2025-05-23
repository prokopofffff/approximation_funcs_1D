#include "interpolation.h"
#include <cmath>
#include <algorithm>
#include <iostream>

Interpolation::Interpolation(double a, double b, int n, int k)
    : a(a), b(b), n(n), currentFunction(k), currentMethod(NEWTON) {
    points.resize(n);
    values.resize(n);
    newtonCoeffs.resize(n);
    akimaCoeffs.resize(4 * (n - 1));  // 4 coefficients per segment
    workArray.resize(n * n);

    setFunction(k);
    initializePoints();
    updateCoefficients();
}

Interpolation::~Interpolation() {
    // Nothing to do here as vectors will be automatically destroyed
}

void Interpolation::setFunction(int k) {
    currentFunction = k;
    updateCoefficients();
}

double Interpolation::calculateFunction(double x) const {
    switch (currentFunction) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return x * x;
        case 3: return x * x * x;
        case 4: return x * x * x * x;
        case 5: return exp(x);
        case 6: {
            // Protect against potential overflow in intermediate calculations
            double denom = 25.0 * x * x + 1.0;
            if (std::abs(denom) < 1e-10) return 1e10; // Return a large value instead of infinity
            return 1.0 / denom;
        }
        default: return 0.0;
    }
}

void Interpolation::initializePoints() {
    // Optimized equidistant points calculation
    double h = (b - a) / (n - 1);
    for (int i = 0; i < n - 1; ++i) {
        points[i] = a + i * h;
        values[i] = calculateFunction(points[i]);
    }
    originalMiddleValue = values[n / 2];
    points[n - 1] = b;
    values[n - 1] = calculateFunction(b);
}

void Interpolation::updateCoefficients(int skipIndex) {
    // Update function values at points
    for (int i = 0; i < n; ++i) {
        if (i == skipIndex) continue;
        values[i] = calculateFunction(points[i]);
        if(i == n / 2) {
            originalMiddleValue = values[i];
        }
    }

    // Update Newton coefficients
    if(n <= 50) {
        constructNewtonPolynomial(n, points.data(), values.data(), newtonCoeffs.data(), workArray.data());
    }

    // Update Akima coefficients
    AkimaSpline::constructSpline(n, points.data(), values.data(), akimaCoeffs.data());
}

double Interpolation::calculateInterpolation(double x, Method method) const {
    switch (method) {
        case NEWTON:
            return evaluateNewtonPolynomial(x, a, b, n, points.data(), newtonCoeffs.data());
        case AKIMA:
            return AkimaSpline::evaluateSpline(x, a, b, n, points.data(), akimaCoeffs.data());
        default:
            return 0.0;
    }
}

void Interpolation::setN(int newN) {
    if (newN < 2) return;

    n = newN;
    points.resize(n);
    values.resize(n);
    newtonCoeffs.resize(std::min(n, 50));
    akimaCoeffs.resize(4 * (n - 1));
    workArray.resize(std::min(n * n, 50 * 50));

    initializePoints();
    updateCoefficients();
}

void Interpolation::setValue(size_t index, double value) {
    if (index >= values.size()) return;

    values[index] = value;
    updateCoefficients(index);
}

void Interpolation::evaluateMultiple(const std::vector<double>& x, std::vector<double>& y, Method method) const {
    size_t numPoints = x.size();
    y.resize(numPoints);

    for (size_t i = 0; i < numPoints; ++i) {
        y[i] = calculateInterpolation(x[i], method);
    }
}

double Interpolation::factorial(int n) {
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}




