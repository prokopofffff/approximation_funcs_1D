#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <QWidget>
#include <QPainter>
#include <functional>

// Function type for calculating function values
using FunctionCalculator = std::function<double(double)>;

// Structure to hold scaling information
struct ScalingInfo {
    double scaleX;
    double scaleY;
    double offsetX;
    double offsetY;
    double maxAbsValue;
};

// Function to draw a function on a QWidget
// Parameters:
//   painter - QPainter object to draw with
//   widget - QWidget to draw on
//   a - left end of the segment
//   b - right end of the segment
//   func - function to calculate values
// Returns:
//   ScalingInfo containing scaling parameters and max absolute value
ScalingInfo drawFunction(QPainter& painter, QWidget* widget, double a, double b, FunctionCalculator func);

// Function to draw coordinate axes
// Parameters:
//   painter - QPainter object to draw with
//   widget - QWidget to draw on
//   scaling - scaling information
//   a - left end of the segment
//   b - right end of the segment
void drawAxes(QPainter& painter, QWidget* widget, const ScalingInfo& scaling, double a, double b);

#endif // GRAPHICS_H 