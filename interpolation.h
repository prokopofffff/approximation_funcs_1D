#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <cmath>
#include "approximation.h"

class Interpolation {
public:
    enum Method {
        NEWTON = 0,
        AKIMA = 1,
        BOTH = 2,
        ERROR = 3
    };

    Interpolation(double a, double b, int n, int k);
    ~Interpolation();

    // Set function by index k
    void setFunction(int k);

    // Calculate function value at point x
    double calculateFunction(double x) const;

    // Calculate interpolation value at point x
    double calculateInterpolation(double x, Method method = NEWTON) const;

    // Get interpolation points
    const std::vector<double>& getPoints() const { return points; }
    const std::vector<double>& getValues() const { return values; }
    const std::vector<double>& getNewtonCoeffs() const { return newtonCoeffs; }
    double getOriginalMiddleValue() const { return originalMiddleValue; }

    // Get function range
    double getA() const { return a; }
    double getB() const { return b; }

    // Get current function index
    int getCurrentFunction() const { return currentFunction; }

    // Get/set number of points
    int getN() const { return n; }
    void setN(int newN);

    // Set value at specific point
    void setValue(size_t index, double value);

    // Optimized evaluation for multiple points
    void evaluateMultiple(const std::vector<double>& x, std::vector<double>& y, Method method = NEWTON) const;

    Method getCurrentMethod() const { return currentMethod; }
    void setCurrentMethod(Method m) { currentMethod = m; }

    void setInterval(double newA, double newB) {
        a = newA;
        b = newB;
        initializePoints();
        updateCoefficients();
    }

private:
    double a, b;
    int n;
    int currentFunction;
    Method currentMethod = NEWTON;
    std::vector<double> points;
    std::vector<double> values;
    std::vector<double> newtonCoeffs;
    std::vector<double> akimaCoeffs;
    std::vector<double> workArray;
    double originalMiddleValue;
    // Initialize interpolation points and values
    void initializePoints();

    // Update polynomial coefficients
    void updateCoefficients(int skipIndex = -1);

    // Calculate factorial
    double factorial(int n);
};

#endif // INTERPOLATION_H
