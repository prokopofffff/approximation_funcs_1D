#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QPainter>
#include <QMouseEvent>
#include <QPainterPath>
#include "interpolation.h"
#include "graphics.h"

class Window : public QWidget {
    Q_OBJECT

public:
    Window(QWidget *parent = nullptr);
    ~Window();

    bool parse_command_line(int argc, char *argv[]);
    void set_func_by_id();

public slots:
    void change_func();  // Key 0
    void change_display_mode();  // Key 1
    void increase_scale();  // Key 2
    void decrease_scale();  // Key 3
    void increase_points();  // Key 4
    void decrease_points();  // Key 5
    void increase_perturbation();  // Key 6
    void decrease_perturbation();  // Key 7

protected:
    void paintEvent(QPaintEvent *event) override;
    void resizeEvent(QResizeEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

private:
    Interpolation *interpolation;
    bool showOriginalFunction;
    bool showInterpolatedFunction;
    bool showError;
    int currentScale;
    int perturbationCount;
    ScalingInfo scaling;
    double origA, origB;

    void drawGraph(QPainter& painter);
    void updateDisplayInfo(QPainter& painter);
    QString getFunctionDescription() const;
    void updateInterpolation();
    void updateScaling();
    void drawFunctionFromPoints(QPainter& painter, const std::vector<double>& x,
                              const std::vector<double>& y, const QColor& color);
};

#endif // WINDOW_H
