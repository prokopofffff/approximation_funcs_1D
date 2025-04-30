#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets>
#include "interpolation.h"

enum DisplayMode {
  // ORIGINAL_AND_CHEBYSHEV, // Show original and Chebyshev interpolation
  ORIGINAL_AND_NEWTON,    // Show original and Newton interpolation
  ORIGINAL_AND_SPLINE,    // Show original and spline interpolation
  ORIGINAL_AND_BOTH,      // Show original and both interpolations
  ERROR_ONLY,             // Show errors of both methods
  DISPLAY_MODE_COUNT      // Number of display modes
};

class Window : public QWidget
{
  Q_OBJECT

public:
  Window (QWidget *parent);
  ~Window();

  QSize minimumSizeHint () const override;
  QSize sizeHint () const override;

  int parse_command_line (int argc, char *argv[]);
  QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);
  
  static Window* active_instance;
  
  void build_interpolation_data();
  
  // double eval_chebyshev(double x_val);
  double eval_newton(double x_val);
  double eval_spline(double x_val);
  
  // double chebyshev_error(double x_val);
  double newton_error(double x_val);
  double spline_error(double x_val);
  
  double (*f) (double);
  
  void set_func_by_id();
  
  // Геттер для получения текущего идентификатора функции
  int get_func_id() const { return func_id; }

public slots:
  void change_func ();
  void change_display_mode();

protected:
  void paintEvent (QPaintEvent *event) override;
  void keyPressEvent (QKeyEvent *event) override;
  void resizeEvent(QResizeEvent *event) override;

private:
  int func_id;
  const char *f_name;
  double a;
  double b;
  int n;

  int scale_factor; 
  double perturb_factor;
  int perturb_point;

  DisplayMode display_mode;
  double *x_points;    
  double *f_values;
  double *f_perturbed;
  // double *chebyshev_coefs;
  double *newton_coefs;
  double *spline_coefs;
  

  double *buffer_x;
  double *buffer_y;
  // double *buffer_cheb;
  double *buffer_newton;
  double *buffer_spline;
  bool need_recalc;
  
  void draw_axes(QPainter &painter, double min_y, double max_y);
  void draw_original_function(QPainter &painter, double min_y, double max_y);
  void draw_chebyshev(QPainter &painter, double min_y, double max_y);
  void draw_newton(QPainter &painter, double min_y, double max_y);
  void draw_spline(QPainter &painter, double min_y, double max_y);
  void draw_errors(QPainter &painter, double min_y, double max_y);
  void draw_info_text(QPainter &painter);
  
  double get_f_value(int idx);
  double get_max_abs_value();
  
  void calculate_residuals(double scaled_a, double scaled_b, double delta_x, double &min_y, double &max_y);
};

#endif
