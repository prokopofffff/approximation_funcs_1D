#include "window.h"
#include <QPainter>
#include <QMessageBox>
#include <QKeyEvent>
#include <cmath>
#include <cstdlib>
#include <vector>

Window::Window(QWidget *parent)
    : QWidget(parent), interpolation(nullptr),
      showOriginalFunction(true), showInterpolatedFunction(true), showError(false),
      currentScale(0), perturbationCount(0) {
    setMinimumSize(900, 700);
    resize(900, 700);
    setFocusPolicy(Qt::StrongFocus); // Enable keyboard focus
    setAttribute(Qt::WA_OpaquePaintEvent); // Optimize painting

    // Initialize with default values
    double a = -1.0;
    double b = 1.0;
    int n = 10;
    int k = 0;

    origA = a;
    origB = b;
    interpolation = new Interpolation(a, b, n, k);

    // Initialize scaling info
    scaling.scaleX = 1.0;
    scaling.scaleY = 1.0;
    scaling.offsetX = 1.0;
    scaling.offsetY = 1.0;
    scaling.maxAbsValue = 1.0;

    // Calculate initial scaling
    updateScaling();
}

Window::~Window() {
    delete interpolation;
}

bool Window::parse_command_line(int argc, char *argv[]) {
    double a = -1.0;  // default left bound
    double b = 1.0;   // default right bound
    int n = 10;       // default number of points
    int k = 0;        // default function index

    // If arguments are provided, use them
    if (argc == 5) {
        a = atof(argv[1]);
        b = atof(argv[2]);
        n = atoi(argv[3]);
        k = atoi(argv[4]);
    }

    // Validate parameters
    if (a >= b || n < 2 || k < 0 || k > 6) {
        return true;
    }

    interpolation = new Interpolation(a, b, n, k);
    return false;
}

void Window::set_func_by_id() {
    if (interpolation) {
        interpolation->setFunction(interpolation->getCurrentFunction());
        update();
    }
}

void Window::change_func() {
    if (interpolation) {
        int current = interpolation->getCurrentFunction();
        current = (current + 1) % 7;
        interpolation->setFunction(current);
        update();
    }
}

void Window::change_display_mode() {
    if (interpolation) {
        static int displayState = 0;
        static int akimaPoints = interpolation->getN();  // Store Akima's point count
        displayState = (displayState + 1) % 4;

        switch (displayState) {
            case 0:  // Newton only
                interpolation->setCurrentMethod(Interpolation::NEWTON);
                interpolation->setN(std::min(80, interpolation->getN()));  // Limit to 80 for Newton
                showInterpolatedFunction = true;
                showError = false;
                break;
            case 1:  // Akima only
                interpolation->setCurrentMethod(Interpolation::AKIMA);
                interpolation->setN(akimaPoints);  // Restore Akima's point count
                showInterpolatedFunction = true;
                showError = false;
                break;
            case 2:  // Both methods
                interpolation->setCurrentMethod(Interpolation::BOTH);
                akimaPoints = interpolation->getN();  // Store current point count for Akima
                interpolation->setN(std::min(80, interpolation->getN()));  // Limit to 80 for Newton
                showInterpolatedFunction = true;
                showError = false;
                break;
            case 3:  // Both methods errors
                interpolation->setCurrentMethod(Interpolation::ERROR);
                showInterpolatedFunction = false;
                showError = true;
                break;
        }
        update();
    }
}

void Window::increase_scale() {
    currentScale++;
    double scaleFactor = std::pow(2.0, currentScale);
    double newA = origA / scaleFactor;
    double newB = origB / scaleFactor;
    interpolation->setInterval(newA, newB);
    update();
}

void Window::decrease_scale() {
    currentScale--;
    double scaleFactor = std::pow(2.0, currentScale);
    double newA = origA / scaleFactor;
    double newB = origB / scaleFactor;
    interpolation->setInterval(newA, newB);
    update();
}

void Window::increase_points() {
    if (interpolation) {
        int maxPoints = (interpolation->getCurrentMethod() == Interpolation::NEWTON) ? 50 : 5000;
        if (interpolation->getN() < maxPoints) {
            interpolation->setN(interpolation->getN() * 2);
            updateInterpolation();
        }
    }
}

void Window::decrease_points() {
    if (interpolation && interpolation->getN() > 2) {
        interpolation->setN(interpolation->getN() / 2);
        updateInterpolation();
    }
}

void Window::increase_perturbation() {
    perturbationCount++;
    updateInterpolation();
}

void Window::decrease_perturbation() {
    perturbationCount--;
    updateInterpolation();
}

void Window::keyPressEvent(QKeyEvent *event) {
    switch (event->key()) {
        case Qt::Key_0:
            change_func();
            break;
        case Qt::Key_1:
            change_display_mode();
            break;
        case Qt::Key_2:
            increase_scale();
            break;
        case Qt::Key_3:
            decrease_scale();
            break;
        case Qt::Key_4:
            increase_points();
            break;
        case Qt::Key_5:
            decrease_points();
            break;
        case Qt::Key_6:
            increase_perturbation();
            break;
        case Qt::Key_7:
            decrease_perturbation();
            break;
        default:
            QWidget::keyPressEvent(event);
    }
}

void Window::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event);
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, false);
    painter.fillRect(rect(), Qt::white);
    if (interpolation) {
        drawGraph(painter);
        updateDisplayInfo(painter);
    }
}

void Window::resizeEvent(QResizeEvent *event) {
    Q_UNUSED(event);
    updateScaling();
    update();
}

void Window::drawGraph(QPainter& painter) {
    painter.setRenderHint(QPainter::Antialiasing, false);
    double a = interpolation->getA();
    double b = interpolation->getB();
    const int numPoints = std::min(200, width());
    std::vector<double> x(numPoints);
    std::vector<double> y(numPoints);
    double step = (b - a) / (numPoints - 1);
    double maxAbsValue = 0.0;
    for (int i = 0; i < numPoints; ++i) {
        x[i] = a + i * step;
        y[i] = interpolation->calculateFunction(x[i]);
        maxAbsValue = std::max(maxAbsValue, std::abs(y[i]));
    }
    scaling.maxAbsValue = maxAbsValue;
    updateScaling();
    if (showOriginalFunction) {
        drawFunctionFromPoints(painter, x, y, Qt::blue);
    }
    if (showInterpolatedFunction) {
        if (interpolation->getCurrentMethod() == Interpolation::NEWTON) {
            interpolation->evaluateMultiple(x, y, Interpolation::NEWTON);
            drawFunctionFromPoints(painter, x, y, Qt::red);
        } else if (interpolation->getCurrentMethod() == Interpolation::AKIMA) {
            interpolation->evaluateMultiple(x, y, Interpolation::AKIMA);
            drawFunctionFromPoints(painter, x, y, Qt::green);
        } else if (interpolation->getCurrentMethod() == Interpolation::BOTH) {
            // Draw both methods
            interpolation->evaluateMultiple(x, y, Interpolation::NEWTON);
            drawFunctionFromPoints(painter, x, y, Qt::red);
            interpolation->evaluateMultiple(x, y, Interpolation::AKIMA);
            drawFunctionFromPoints(painter, x, y, Qt::green);
        }
    }
    if (showError) {
        // Calculate and draw Newton error only if points <= 80
        if (interpolation->getN() <= 80) {
            #pragma omp parallel for
            for (int i = 0; i < numPoints; ++i) {
                y[i] = interpolation->calculateInterpolation(x[i], Interpolation::NEWTON) - interpolation->calculateFunction(x[i]);
            }
            drawFunctionFromPoints(painter, x, y, Qt::magenta);
        }

        // Calculate and draw Akima error
        #pragma omp parallel for
        for (int i = 0; i < numPoints; ++i) {
            y[i] = interpolation->calculateInterpolation(x[i], Interpolation::AKIMA) - interpolation->calculateFunction(x[i]);
        }
        drawFunctionFromPoints(painter, x, y, Qt::cyan);
    }
    drawAxes(painter, this, scaling, a, b);
}

void Window::drawFunctionFromPoints(QPainter& painter, const std::vector<double>& x,
                                  const std::vector<double>& y, const QColor& color) {
    double scaleX = scaling.scaleX;
    double scaleY = scaling.scaleY;
    double offsetX = scaling.offsetX;
    double offsetY = scaling.offsetY;
    int h = height();
    painter.setPen(QPen(color, 2));
    int prevX = 0, prevY = 0;
    bool first = true;
    for (size_t i = 0; i < x.size(); ++i) {
        int screenX = static_cast<int>(x[i] * scaleX + offsetX);
        int screenY = h - static_cast<int>(y[i] * scaleY + offsetY);
        if (first) {
            prevX = screenX;
            prevY = screenY;
            first = false;
        } else {
            painter.drawLine(prevX, prevY, screenX, screenY);
            prevX = screenX;
            prevY = screenY;
        }
    }
}

void Window::updateDisplayInfo(QPainter& painter) {
    painter.setFont(QFont("Arial", 10, QFont::Bold));
    int y = 20;
    // Display function description
    QString funcDesc = getFunctionDescription();
    painter.drawText(10, y, funcDesc);
    y += 20;
    // Display interval
    if (interpolation) {
        double a = interpolation->getA();
        double b = interpolation->getB();
        QString intervalText = QString("Interval: [%1, %2]").arg(a).arg(b);
        painter.drawText(10, y, intervalText);
        y += 20;
    }
    // Display current approximation method
    QString methodText = QString("Method: %1").arg(
        interpolation->getCurrentMethod() == Interpolation::NEWTON ? "Newton" :
        interpolation->getCurrentMethod() == Interpolation::AKIMA ? "Akima" :
        interpolation->getCurrentMethod() == Interpolation::BOTH ? "Both" : "Error");
    painter.drawText(10, y, methodText);
    y += 20;
    // Display number of points
    if (interpolation) {
        QString pointsText = QString("Points: %1").arg(interpolation->getN());
        painter.drawText(10, y, pointsText);
        y += 20;
    }
    // Display current scale
    QString scaleText = QString("Scale: 2^%1").arg(currentScale);
    painter.drawText(10, y, scaleText);
    y += 20;
    // Display perturbation
    QString pertText = QString("Perturbation: %1").arg(perturbationCount);
    painter.drawText(10, y, pertText);
    y += 20;

    // Display plot legends
    if (showOriginalFunction) {
        painter.fillRect(10, y, 15, 15, Qt::blue);
        painter.drawText(30, y + 12, "Original function");
        y += 20;
    }
    if (showInterpolatedFunction) {
        if (interpolation->getCurrentMethod() == Interpolation::NEWTON) {
            painter.fillRect(10, y, 15, 15, Qt::red);
            painter.drawText(30, y + 12, "Newton");
        } else if (interpolation->getCurrentMethod() == Interpolation::AKIMA) {
            painter.fillRect(10, y, 15, 15, Qt::green);
            painter.drawText(30, y + 12, "Akima");
        } else if (interpolation->getCurrentMethod() == Interpolation::BOTH) {
            painter.fillRect(10, y, 15, 15, Qt::red);
            painter.drawText(30, y + 12, "Newton");
            y += 20;
            painter.fillRect(10, y, 15, 15, Qt::green);
            painter.drawText(30, y + 12, "Akima");
        }
    }
    if (showError) {
        painter.fillRect(10, y, 15, 15, Qt::magenta);
        painter.drawText(30, y + 12, "Newton error");
        y += 20;
        painter.fillRect(10, y, 15, 15, Qt::cyan);
        painter.drawText(30, y + 12, "Akima error");
    }
}

QString Window::getFunctionDescription() const {
    if (!interpolation) return "";

    int k = interpolation->getCurrentFunction();
    QString desc;
    switch (k) {
        case 0: desc = "f(x) = 1"; break;
        case 1: desc = "f(x) = x"; break;
        case 2: desc = "f(x) = x^2"; break;
        case 3: desc = "f(x) = x^3"; break;
        case 4: desc = "f(x) = x^4"; break;
        case 5: desc = "f(x) = e^x"; break;
        case 6: desc = "f(x) = 1/(25x^2 + 1)"; break;
        default: desc = "Unknown function"; break;
    }
    return QString("k = %1, %2").arg(k).arg(desc);
}

void Window::updateInterpolation() {
    if (!interpolation) return;

    // Apply perturbation if needed
    if (perturbationCount != 0) {
        double maxValue = 0.0;
        const auto& points = interpolation->getPoints();
        const auto& values = interpolation->getValues();

        // Find maximum absolute value
        for (size_t i = 0; i < points.size(); ++i) {
            maxValue = std::max(maxValue, std::abs(values[i]));
        }

        // Apply perturbation to middle point
        size_t middleIndex = points.size() / 2;
        double perturbation = perturbationCount * 0.1 * maxValue;
        interpolation->setValue(middleIndex, values[middleIndex] + perturbation);
    }

    update();
}

void Window::updateScaling() {
    if (!interpolation) return;
    int w = width();
    int h = height();
    int margin = 50;
    double a = interpolation->getA(); // current interval
    double b = interpolation->getB(); // current interval
    // Sample the function to find min and max Y
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    int numPoints = std::min(200, w);
    double step = (b - a) / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i) {
        double x = a + i * step;
        double y = interpolation->calculateFunction(x);
        minY = std::min(minY, y);
        maxY = std::max(maxY, y);
    }
    if (maxY - minY < 1e-12) {
        maxY += 1.0;
        minY -= 1.0;
    }
    scaling.maxAbsValue = std::max(std::abs(minY), std::abs(maxY));
    scaling.scaleX = (w - 2 * margin) / (b - a); // fill width
    scaling.scaleY = (h - 2 * margin) / (maxY - minY);
    scaling.offsetX = margin - a * scaling.scaleX;
    scaling.offsetY = margin - minY * scaling.scaleY;
}
