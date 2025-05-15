#include "graphics.h"
#include <QFont>
#include <QPen>
#include <QPainterPath>
#include <cmath>
#include <iostream>

ScalingInfo drawFunction(QPainter& painter, QWidget* widget, double a, double b, FunctionCalculator func) {
    ScalingInfo scaling;
    
    // Calculate the maximum absolute value of the function
    const int numPoints = 10;
    double maxAbsValue = 0.0;
    double step = (b - a) / numPoints;
    
    for (int i = 0; i <= numPoints; ++i) {
        double x = a + i * step;
        double y = func(x);
        maxAbsValue = std::max(maxAbsValue, std::abs(y));
    }
    
    // Add 10% margin to the maximum value
    maxAbsValue *= 1.0;
    
    // Calculate scaling factors
    int widgetWidth = widget->width();
    int widgetHeight = widget->height();
    
    // Leave margins for axes and labels
    int margin = 0;
    scaling.scaleX = (widgetWidth - 2 * margin) / (b - a);
    scaling.scaleY = (widgetHeight - 2 * margin) / (2 * maxAbsValue);
    
    // Calculate offsets
    scaling.offsetX = margin;
    scaling.offsetY = widgetHeight / 2;
    
    // Store maximum absolute value
    scaling.maxAbsValue = maxAbsValue;
    
    // Draw the function
    QPainterPath path;
    bool first = true;
    
    for (int i = 0; i <= numPoints; ++i) {
        double x = a + i * step;
        double y = func(x);
        
        int screenX = static_cast<int>((x - a) * scaling.scaleX + scaling.offsetX);
        int screenY = static_cast<int>(widgetHeight - (y + maxAbsValue) * scaling.scaleY - scaling.offsetY);
        
        if (first) {
            path.moveTo(screenX, screenY);
            first = false;
        } else {
            path.lineTo(screenX, screenY);
        }
    }
    
    painter.setPen(QPen(Qt::blue, 2));
    painter.drawPath(path);
    
    // Output maximum absolute value
    std::cout << "Maximum absolute value of the function: " << maxAbsValue << std::endl;
    
    return scaling;
}

void drawAxes(QPainter& painter, QWidget* widget, const ScalingInfo& scaling, double a, double b) {
    int w = widget->width();
    int h = widget->height();
    painter.setPen(QPen(Qt::black, 1));

    // Pre-calculate scaling factors
    double scaleX = scaling.scaleX;
    double scaleY = scaling.scaleY;
    double offsetX = scaling.offsetX;
    double offsetY = scaling.offsetY;

    // Compute minY and maxY in screen coordinates
    double minY = (h - offsetY) / scaleY;
    double maxY = (0 - offsetY) / scaleY;
    if (minY > maxY) std::swap(minY, maxY);

    // Draw Y axis at x=0 if 0 is in [a, b]
    if (a <= 0 && b >= 0) {
        int zeroX = static_cast<int>(0 * scaleX + offsetX);
        painter.drawLine(zeroX, 0, zeroX, h);
    }
    // Draw X axis at y=0 if 0 is in [minY, maxY]
    if (minY <= 0 && maxY >= 0) {
        int zeroY = h - static_cast<int>(0 * scaleY + offsetY);
        painter.drawLine(0, zeroY, w, zeroY);
    }

    // Draw ticks and labels as before...
    const int numTicks = 5;
    double step = (b - a) / (numTicks - 1);
    QFont font = painter.font();
    font.setPointSize(8);
    painter.setFont(font);
    QFontMetrics fm(font);
    for (int i = 0; i < numTicks; ++i) {
        double x = a + i * step;
        int screenX = static_cast<int>(x * scaleX + offsetX);
        int zeroY = h - static_cast<int>(0 * scaleY + offsetY);
        painter.drawLine(screenX, zeroY - 5, screenX, zeroY + 5);
        QString label = QString::number(x, 'g', 3);
        int labelWidth = fm.horizontalAdvance(label);
        painter.drawText(screenX - labelWidth/2, zeroY + 20, label);
    }
    // Y-axis ticks and labels
    double yMin = minY, yMax = maxY;
    double yStep = (yMax - yMin) / (numTicks - 1);
    for (int i = 0; i < numTicks; ++i) {
        double y = yMin + i * yStep;
        int zeroX = static_cast<int>(0 * scaleX + offsetX);
        int screenY = h - static_cast<int>(y * scaleY + offsetY);
        painter.drawLine(zeroX - 5, screenY, zeroX + 5, screenY);
        QString label = QString::number(y, 'g', 3);
        int labelWidth = fm.horizontalAdvance(label);
        painter.drawText(zeroX - labelWidth - 10, screenY + fm.height()/4, label);
    }
} 