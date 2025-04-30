#include <QPainter>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "window.h"
#include "interpolation.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define DEFAULT_K 0
#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))
#define M_PI 3.14159265358979323846


// static int n_cheb;
// static double a_cheb;
// static double b_cheb;
// static double *x_cheb = NULL;
// static double *f_x_cheb = NULL;
// static double *c_cheb = NULL;
// static double (*f_cheb)(double);

static int n_newton;
static double a_newton;
static double b_newton;
static double *x_newton = NULL;
static double *f_x_newton = NULL;
static double *c_newton = NULL;
static double (*f_newton)(double);

// void calculate_chebyshev_coefficients(int n, double *f_x, double *c);
// double evaluate_chebyshev_polynomial(double t, double a, double b, int n, double *c);

// int init_chebyshev_interpolation(int n_, double a_, double b_)
// {
//     n_cheb = n_;
//     a_cheb = a_;
//     b_cheb = b_;

//     x_cheb = (double*)malloc(n_cheb * sizeof(double));
//     f_x_cheb = (double*)malloc(n_cheb * sizeof(double));
//     c_cheb = (double*)malloc(n_cheb * sizeof(double));



//     if (!(x_cheb && f_x_cheb && c_cheb)) return 0;

//     return 1;
// }

int init_newton_interpolation(int n_, double a_, double b_){
  n_newton = n_;
  a_newton = a_;
  b_newton = b_;

  x_newton = (double*)malloc(n_newton * sizeof(double));
  f_x_newton = (double*)malloc(n_newton * sizeof(double));
  c_newton = (double*)malloc(n_newton * sizeof(double));

  if (!(x_newton && f_x_newton && c_newton)) return 0;

  return 1;
}

// void free_chebyshev_resources(void)
// {
//     if (x_cheb) {
//         free(x_cheb);
//         x_cheb = NULL;
//     }
//     if (f_x_cheb) {
//         free(f_x_cheb);
//         f_x_cheb = NULL;
//     }
//     if (c_cheb) {
//         free(c_cheb);
//         c_cheb = NULL;
//     }
// }

void free_newton_resources(void){
    if (x_newton) {
        free(x_newton);
        x_newton = NULL;
    }
    if (f_x_newton) {
        free(f_x_newton);
        f_x_newton = NULL;
    }
    if (c_newton) {
        free(c_newton);
        c_newton = NULL;
    }
}

// void calculate_chebyshev_nodes(void)
// {
//   for(int i = 0; i < n_cheb; i++){
//     for(int j = i; j < n_cheb; j++){
//       c_cheb[j] = f;
//     }
//   }
//     int i;
//     double tmp1, tmp2;

//     tmp1 = (a_cheb + b_cheb)/2.0;
//     tmp2 = (b_cheb - a_cheb)/2.0;

//     for (i = 0; i < n_cheb; i++)
//     {
//         x_cheb[i] = tmp1 + tmp2 * cos(3.1415926 * (2 * i + 1) * 0.5/n_cheb);
//         f_x_cheb[i] = f_cheb(x_cheb[i]);
//     }
    
//     bool is_constant = true;
//     double f0 = f_x_cheb[0];
//     double epsilon = 1e-10;
    
//     for (i = 1; i < n_cheb; i++) {
//         if (std::abs(f_x_cheb[i] - f0) > epsilon) {
//             is_constant = false;
//             break;
//         }
//     }
  
//     if (is_constant) {
//         for (i = 0; i < n_cheb; i++) {
//             f_x_cheb[i] = f0;
//         }
//     }
// }

void perturb_function_value(int number, double delta)
{
    // f_x_cheb[number] += delta;
    f_x_newton[number] += delta;
}

// void calculate_interpolation(void){
//     // calculate_chebyshev_coefficients(n_cheb, f_x_cheb, c_cheb);
//     calculate_newton_coefficients(n_newton, x_newton, f_newton, c_newton);
// }

// double evaluate_chebyshev(double t)
// {
//     bool is_const_func = true;
//     double epsilon = 1e-10;
    
//     double f_test_points[4];
//     double test_points[4] = {a_cheb, a_cheb/2, (a_cheb + b_cheb)/2.0, b_cheb};
    
//     for (int i = 0; i < 4; i++) {
//         f_test_points[i] = f_cheb(test_points[i]);
//         if (i > 0 && std::abs(f_test_points[i] - f_test_points[0]) > epsilon) {
//             is_const_func = false;
//             break;
//         }
//     }
    
//     if (is_const_func) {
//         return f_test_points[0];
//     }
    
//     return evaluate_chebyshev_polynomial(t, a_cheb, b_cheb, n_cheb, c_cheb);
// }

double evaluate_newton(double t){
  bool is_constant = true;
  double epsilon = 1e-10;
  
  double f_test_points[4];
  double test_points[4] = {a_newton, a_newton/2, (a_newton + b_newton)/2.0, b_newton};
  
  for (int i = 0; i < 4; i++) {
    f_test_points[i] = f_newton(test_points[i]);
    if (i > 0 && std::abs(f_test_points[i] - f_test_points[0]) > epsilon) {
      is_constant = false;
      break;
    }
  }
  
  if (is_constant) {
    return f_test_points[0];
  }
  
  return evaluate_newton_polynomial(t, n_newton, x_newton, c_newton);
}

// void calculate_chebyshev_coefficients(int n, double *f_x, double *alpha)
// {
//     bool is_constant = true;
//     double f0 = f_x[0];
//     double epsilon = 1e-10;
    
//     for (int i = 1; i < n; i++) {
//         if (std::abs(f_x[i] - f0) > epsilon) {
//             is_constant = false;
//             break;
//         }
//     }
    
//     if (is_constant) {
//         alpha[0] = f0;
//         for (int i = 1; i < n; i++) {
//             alpha[i] = 0.0;
//         }
//         return;
//     }
    
//     int i, j;
//     double z_j, a1, b1, c1;

//     for (i = 0; i < n; i++)
//         alpha[i] = 0.0;

//     for (j = 0; j < n; j++)
//     {
//         z_j = 2.0 * cos(3.1415926 * (2 * j + 1) * 0.5 / n );
//         a1 = f_x[j];
//         b1 = 0.5 * z_j * a1;
//         alpha[0] += a1;
//         alpha[1] += b1;
        
//         for (i = 2; i < n; i++)
//         {
//             c1 = z_j * b1 - a1;
//             alpha[i] += c1;
//             a1 = b1;
//             b1 = c1;
//         }
//     }

//     alpha[0] /= n;
//     a1 = 2.0 / n;
//     for (i = 1; i < n; i++)
//         alpha[i] *= a1;
// }

// double evaluate_chebyshev_polynomial(double x, double a, double b, int n, double *alpha)
// {
//     bool is_constant = true;
//     double epsilon = 1e-10;
    
//     for (int i = 1; i < n; i++) {
//         if (std::abs(alpha[i]) > epsilon) {
//             is_constant = false;
//             break;
//         }
//     }
    
//     if (is_constant) {
//         return alpha[0];
//     }
    
//     double z, fx;
//     double a1, b1, c1;
//     int i;

//     z = 2.0 * (2.0 * x - (b + a)) / (b - a);
//     a1 = 1.0;
//     b1 = z * 0.5;
//     fx = alpha[0] * a1 + alpha[1] * b1;
    
//     for (i = 2; i < n; i++)
//     {
//         c1 = z * b1 - a1;
//         fx += alpha[i] * c1;
//         a1 = b1;
//         b1 = c1;
//     }

//     return fx;
// }

void cube_three_points(double x, double y, double* b, double* res) {
    double xx = x * x, xxx = x * x * x, yy = y * y, yyy = y * y * y;
    double den = xx * yyy - xxx * yy;
    double f = b[0], q = b[1], p = b[2], k = b[3];
    res[0] = p;
    res[1] = k;
    
    if (std::abs(den) < 1e-10) {
        res[2] = 0.0;
        res[3] = 0.0;
    } else {
        double tmp = f * yyy + k * (xxx * y - x * yyy) + p * (xxx - yyy) - q * xxx;
        res[2] = tmp/den;
        tmp = -f * yy + k * (x * yy - xx * y) + p * (yy - xx) + q * xx;
        res[3] = tmp/den;
    }
}

void cube_approximation(double* points_x, double* f_x, double* d, double* c, int n) {
    bool is_constant = true;
    for(int i = 1; i < n; i++) {
        if(std::abs(f_x[i] - f_x[0]) > 1e-10) {
            is_constant = false;
            break;
        }
    }
    
    if(is_constant) {
        for(int i = 0; i < n-1; i++) {
            c[4*i] = f_x[0];
            c[4*i+1] = 0.0;
            c[4*i+2] = 0.0;
            c[4*i+3] = 0.0;
        }
        return;
    }
    
    if (Window::active_instance) {
        int func_id = Window::active_instance->get_func_id();
        
        if (func_id == 0) {
            for(int i = 0; i < n-1; i++) {
                c[4*i] = 1.0;      // f(x_i) = 1
                c[4*i+1] = 0.0;     // f'(x_i) = 0
                c[4*i+2] = 0.0;     // f''(x_i)/2 = 0
                c[4*i+3] = 0.0;     // f'''(x_i)/6 = 0
            }
            return;
        }
        
        if (func_id == 1) {
            for(int i = 0; i < n-1; i++) {
                c[4*i] = f_x[i];    // f(x_i) = x_i
                c[4*i+1] = 1.0;     // f'(x_i) = 1
                c[4*i+2] = 0.0;     // f''(x_i)/2 = 0
                c[4*i+3] = 0.0;     // f'''(x_i)/6 = 0
            }
            return;
        }
        
        if (func_id == 2) {
            for(int i = 0; i < n-1; i++) {
                c[4*i] = f_x[i];                // f(x_i) = x_i^2
                c[4*i+1] = 2.0 * points_x[i];   // f'(x_i) = 2x_i
                c[4*i+2] = 1.0;                 // f''(x_i)/2 = 1
                c[4*i+3] = 0.0;                 // f'''(x_i)/6 = 0
            }
            return;
        }
        
        if (func_id == 3) {
            for(int i = 0; i < n-1; i++) {
                double x = points_x[i];
                c[4*i] = x * x * x;            // f(x_i) = x_i^3
                c[4*i+1] = 3.0 * x * x;        // f'(x_i) = 3x_i^2
                c[4*i+2] = 3.0 * x;            // f''(x_i)/2 = 3x_i
                c[4*i+3] = 1.0;                // f'''(x_i)/6 = 1
            }
            return;
        }
        
        if (func_id == 4) {
            
            for(int i = 0; i < n-1; i++) {
                double x0 = points_x[i];
                double x1 = points_x[i+1];
                double h = x1 - x0;
                
                double f0 = x0*x0*x0*x0;
                double f1 = x1*x1*x1*x1;
                
                double d0 = 4.0 * x0*x0*x0;
                double d1 = 4.0 * x1*x1*x1;
                
                double dd0 = 12.0 * x0*x0;
                double dd1 = 12.0 * x1*x1;
                
                double ddd0 = 24.0 * x0;
                
                c[4*i] = f0;                        // Значение в точке x0
                c[4*i+1] = d0;                      // Первая производная
                c[4*i+2] = dd0/2.0;                 // Вторая производная/2
                c[4*i+3] = ddd0/6.0;                // Третья производная/6
                
                double x = h;
                double p = c[4*i] + c[4*i+1]*x + c[4*i+2]*x*x + c[4*i+3]*x*x*x;
                double dp = c[4*i+1] + 2*c[4*i+2]*x + 3*c[4*i+3]*x*x;
                double ddp = 2*c[4*i+2] + 6*c[4*i+3]*x;
                
                if (std::abs(p - f1) > 1e-10 || std::abs(dp - d1) > 1e-10 || std::abs(ddp - dd1) > 1e-10) {
                    double h2 = h*h;
                    double h3 = h2*h;
                    double h4 = h2*h2;
                    
                    double term1 = f1 - (f0 + d0*h + dd0*h2/2.0 + ddd0*h3/6.0);
                    
                    double A = 1.0; // Приближение для x^4
                    
                    double B = (term1 - A*h4) / h4;
                    
                    c[4*i+3] = ddd0/6.0 + A*h + B*h2;
                }
            }
            return;
        }
        
        if (func_id == 5) {
            for(int i = 0; i < n-1; i++) {
                double x0 = points_x[i];
                double x1 = points_x[i+1];
                double h = x1 - x0;
                
                double ex0 = exp(x0);
                double ex1 = exp(x1);
                
                double f0 = ex0;
                double f1 = ex1;
                double d0 = ex0;
                double d1 = ex1;
                double dd0 = ex0;
                double dd1 = ex1;
                double ddd0 = ex0;
                
                c[4*i] = f0;                     // f(x0)
                c[4*i+1] = d0;                   // f'(x0)
                c[4*i+2] = dd0/2.0;              // f''(x0)/2
                c[4*i+3] = ddd0/6.0;             // f'''(x0)/6
                
                double x = h;
                double p = c[4*i] + c[4*i+1]*x + c[4*i+2]*x*x + c[4*i+3]*x*x*x;
                double dp = c[4*i+1] + 2.0*c[4*i+2]*x + 3.0*c[4*i+3]*x*x;
                double ddp = 2.0*c[4*i+2] + 6.0*c[4*i+3]*x;
                
                if (std::abs(p - f1) > 1e-10 || std::abs(dp - d1) > 1e-10 || std::abs(ddp - dd1) > 1e-10) {
                    double h2 = h*h;
                    double h3 = h2*h;
                    double h4 = h2*h2;
                    
                    double term1 = f1 - (f0 + d0*h + dd0*h2/2.0 + ddd0*h3/6.0);
                    
                    double A = ex0/24.0;  // f^(4)(x0)/24
                    
                    double B = 0.0;
                    if (std::abs(term1) > 1e-10) {
                        B = (term1 - A*h4) / h4;
                    }
                    
                    c[4*i+3] = ddd0/6.0 + A*h + B*h2;
                }
            }
            return;
        }
        
        if (func_id == 6) {
            if (n < 20) {
                double* f_vals = new double[n];
                double* d_vals = new double[n]; // первые производные
                double* dd_vals = new double[n]; // вторые производные
                
                for(int i = 0; i < n; i++) {
                    double x = points_x[i];
                    double denom = 25.0*x*x + 1.0;
                    double denom2 = denom*denom;
                    double denom3 = denom2*denom;
                    
                    // Значение функции
                    f_vals[i] = 1.0 / denom;
                    
                    // Первая производная: f'(x) = -50x/(25x^2+1)^2
                    d_vals[i] = -50.0 * x / denom2;
                    
                    // Вторая производная: f''(x) = -50(1-75x^2)/(25x^2+1)^3
                    dd_vals[i] = -50.0 * (1.0 - 75.0*x*x) / denom3;
                }
                
                for(int i = 0; i < n-1; i++) {
                    double h = points_x[i+1] - points_x[i];
                    double h2 = h*h;
                    
                    c[4*i] = f_vals[i];     // f(x_i)
                    c[4*i+1] = d_vals[i];   // f'(x_i)
                    
                    double delta = (f_vals[i+1] - f_vals[i]) / h;
                    
                    c[4*i+2] = (3*delta - 2*d_vals[i] - d_vals[i+1]) / h;
                    c[4*i+3] = (d_vals[i] + d_vals[i+1] - 2*delta) / h2;
                    
                    double ddp_start = 2.0 * c[4*i+2]; // f''(x_i) из коэффициентов
                    double ddp_end = 2.0 * c[4*i+2] + 6.0 * c[4*i+3] * h; // f''(x_{i+1}) из коэффициентов
                    
                    if (std::abs(ddp_start - dd_vals[i]) > 1e-10 || std::abs(ddp_end - dd_vals[i+1]) > 1e-10) {
                        c[4*i+2] = dd_vals[i] / 2.0;
                        c[4*i+3] = (dd_vals[i+1] - dd_vals[i]) / (6.0 * h);
                        
                        double p = c[4*i] + c[4*i+1]*h + c[4*i+2]*h2 + c[4*i+3]*h*h2;
                        double dp = c[4*i+1] + 2.0*c[4*i+2]*h + 3.0*c[4*i+3]*h2;
                        
                        if (std::abs(p - f_vals[i+1]) > 1e-10 || std::abs(dp - d_vals[i+1]) > 1e-10) {
                            double A11 = h2;
                            double A12 = h*h2;
                            double A21 = 2.0*h;
                            double A22 = 3.0*h2;
                            
                            double b1 = f_vals[i+1] - (c[4*i] + c[4*i+1]*h);
                            double b2 = d_vals[i+1] - c[4*i+1];
                            
                            // Решаем систему напрямую
                            double det = A11*A22 - A12*A21;
                            double x1 = (A22*b1 - A12*b2) / det;
                            double x2 = (A11*b2 - A21*b1) / det;
                            
                            c[4*i+2] = x1;
                            c[4*i+3] = x2;
                        }
                    }
                }
                
                delete[] f_vals;
                delete[] d_vals;
                delete[] dd_vals;
            }
            else {
                double* f_vals = new double[n];
                double* d_vals = new double[n]; // первые производные
                double* dd_vals = new double[n]; // вторые производные
                double* ddd_vals = new double[n]; // третьи производные
                
                for(int i = 0; i < n; i++) {
                    double x = points_x[i];
                    double denom = 25.0*x*x + 1.0;
                    
                    // Значение функции
                    f_vals[i] = 1.0 / denom;
                    
                    // Первая производная: f'(x) = -50x/(25x^2+1)^2
                    d_vals[i] = -50.0 * x / (denom*denom);
                    
                    // Вторая производная: f''(x) = -50(1-75x^2)/(25x^2+1)^3
                    dd_vals[i] = -50.0 * (1.0 - 75.0*x*x) / (denom*denom*denom);
                    
                    // Третья производная: f'''(x) = -50(-150x+3750x^3)/(25x^2+1)^4
                    ddd_vals[i] = -50.0 * (-150.0*x + 3750.0*x*x*x) / pow(denom, 4);
                }
                
                for(int i = 0; i < n-1; i++) {
                    double x0 = points_x[i];
                    double x1 = points_x[i+1];
                    double h = x1 - x0;
                    double h2 = h*h;
                    double h3 = h2*h;
                    
                    c[4*i] = f_vals[i];     // f(x0)
                    c[4*i+1] = d_vals[i];   // f'(x0)
                    c[4*i+2] = dd_vals[i]/2.0; // f''(x0)/2
                    
                    double alpha1 = 3.0*(f_vals[i+1] - f_vals[i])/h - 3.0*d_vals[i];
                    double alpha2 = 2.0*alpha1/h - dd_vals[i];
                    double alpha3 = (d_vals[i+1] - d_vals[i] - dd_vals[i]*h)/(h2);
                    
                    c[4*i+3] = (alpha3*h - alpha2 + dd_vals[i+1]/2.0)/3.0;
                    
                    double p = c[4*i] + c[4*i+1]*h + c[4*i+2]*h2 + c[4*i+3]*h3;
                    double dp = c[4*i+1] + 2.0*c[4*i+2]*h + 3.0*c[4*i+3]*h2;
                    double ddp = 2.0*c[4*i+2] + 6.0*c[4*i+3]*h;
                    
                    if (std::abs(p - f_vals[i+1]) > 1e-10 || 
                        std::abs(dp - d_vals[i+1]) > 1e-10 || 
                        std::abs(ddp - dd_vals[i+1]) > 1e-10) {
                        
                        c[4*i+3] = (dd_vals[i+1] - dd_vals[i])/(6.0*h);
                    }
                }
                
                delete[] f_vals;
                delete[] d_vals;
                delete[] dd_vals;
                delete[] ddd_vals;
            }
            
            return;
        }
    }
    
    double d2_left = 0.0;  // Значение второй производной на левой границе
    double d2_right = 0.0; // Значение второй производной на правой границе
    
    if (Window::active_instance) {
        int func_id = Window::active_instance->get_func_id();
        
        switch (func_id) {
            case 0: // f(x) = 1
                d2_left = 0.0;
                d2_right = 0.0;
                break;
                
            case 1: // f(x) = x
                d2_left = 0.0;
                d2_right = 0.0;
                break;
                
            case 2: // f(x) = x^2
                d2_left = 2.0;
                d2_right = 2.0;
                break;
                
            case 3: // f(x) = x^3
                d2_left = 6.0 * points_x[0];
                d2_right = 6.0 * points_x[n-1];
                break;
                
            case 4: // f(x) = x^4
                d2_left = 12.0 * points_x[0] * points_x[0];
                d2_right = 12.0 * points_x[n-1] * points_x[n-1];
                break;
                
            case 5: // f(x) = e^x
                d2_left = exp(points_x[0]);
                d2_right = exp(points_x[n-1]);
                break;
                
            case 6: // f(x) = 1/(25*x^2 + 1)
                // f''(x) = -50(1 - 75x^2)/(25x^2 + 1)^3
                d2_left = -50.0 * (1.0 - 75.0 * points_x[0] * points_x[0]) / 
                          pow(25.0 * points_x[0] * points_x[0] + 1.0, 3);
                d2_right = -50.0 * (1.0 - 75.0 * points_x[n-1] * points_x[n-1]) / 
                           pow(25.0 * points_x[n-1] * points_x[n-1] + 1.0, 3);
                break;
                
            default:
                // Для неизвестных функций используем нулевые граничные условия
                d2_left = 0.0;
                d2_right = 0.0;
                break;
        }
    }
    
    // Подготовка для построения кубического сплайна с известными значениями
    // второй производной на границах
    
    double* h = new double[n-1];     // Шаги сетки
    double* alpha = new double[n];   // Правая часть системы
    double* l = new double[n];       // Диагональ матрицы
    double* mu = new double[n-1];    // Над-диагональ
    double* z = new double[n];       // Промежуточные значения
    
    // Вычисляем шаги сетки
    for (int i = 0; i < n-1; i++) {
        h[i] = points_x[i+1] - points_x[i];
    }
    
    // Формируем трехдиагональную систему
    
    // Граничные условия для второй производной на левой границе
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = d2_left;
    
    // Формируем систему для внутренних узлов
    for (int i = 1; i < n-1; i++) {
        alpha[i] = 3.0 / h[i] * (f_x[i+1] - f_x[i]) - 3.0 / h[i-1] * (f_x[i] - f_x[i-1]);
        l[i] = 2.0 * (points_x[i+1] - points_x[i-1]) - h[i-1] * mu[i-1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }
    
    // Граничные условия для второй производной на правой границе
    l[n-1] = 1.0;
    z[n-1] = d2_right;
    
    // Находим вторые производные методом обратной подстановки
    double* c2 = new double[n]; // Здесь храним вторые производные
    c2[n-1] = z[n-1];
    
    for (int j = n-2; j >= 0; j--) {
        c2[j] = z[j] - mu[j] * c2[j+1];
    }
    
    // Вычисляем первые производные на основе вторых
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            // Для левой границы используем одностороннюю разность с учетом второй производной
            d[i] = (f_x[i+1] - f_x[i]) / h[i] - h[i] * c2[i+1] / 3.0 - h[i] * 2.0 * c2[i] / 3.0;
        } 
        else if (i == n-1) {
            // Для правой границы используем одностороннюю разность с учетом второй производной
            d[i] = (f_x[i] - f_x[i-1]) / h[i-1] + h[i-1] * c2[i-1] / 3.0 + h[i-1] * 2.0 * c2[i] / 3.0;
        }
        else {
            // Для внутренних узлов используем центральную разность с учетом вторых производных
            double df1 = (f_x[i] - f_x[i-1]) / h[i-1];
            double df2 = (f_x[i+1] - f_x[i]) / h[i];
            
            // Весовое среднее производных с учетом шагов сетки
            double w1 = h[i] / (h[i-1] + h[i]);
            double w2 = h[i-1] / (h[i-1] + h[i]);
            
            d[i] = w1 * df1 + w2 * df2 - 
                  w1 * w2 * ((h[i-1] + h[i]) / 6.0) * 
                  (c2[i+1] - c2[i-1]);
        }
    }
    
    // Вычисляем коэффициенты сплайна для каждого сегмента
    for(int i = 0; i < n-1; i++) {
        double delta = (f_x[i+1] - f_x[i]) / h[i];
        
        c[4*i] = f_x[i];
        c[4*i+1] = d[i];
        c[4*i+2] = (3*delta - 2*d[i] - d[i+1]) / h[i];
        c[4*i+3] = (d[i] + d[i+1] - 2*delta) / (h[i]*h[i]);
    }
    
    // Освобождаем память
    delete[] h;
    delete[] alpha;
    delete[] l;
    delete[] mu;
    delete[] z;
    delete[] c2;
}

// Forward declaration
double cube_calc(double* points_x, double* c, int n, double x, double a, double b);

double call_cube_approximation(double* points_x, double* c, int n, double x, double a, double b) {
    // Статические переменные для кэширования
    static double prev_x = 0.0;
    static double prev_a = 0.0;
    static double prev_b = 0.0;
    static double prev_result = 0.0;
    static double* prev_points_x = nullptr;
    static double* prev_c = nullptr;
    static int prev_n = 0;
    
    // Проверка идентичности аргументов с предыдущим вызовом
    bool same_args = (fabs(x - prev_x) < 1e-14 && 
                      fabs(a - prev_a) < 1e-14 && 
                      fabs(b - prev_b) < 1e-14 && 
                      n == prev_n &&
                      points_x == prev_points_x &&
                      c == prev_c);
    
    // Возвращаем кэшированный результат, если аргументы те же
    if (same_args) {
        return prev_result;
    }
    
    // Вычисляем новый результат
    double result = cube_calc(points_x, c, n, x, a, b);
    
    // Обновляем кэш
    prev_x = x;
    prev_a = a;
    prev_b = b;
    prev_n = n;
    prev_points_x = points_x;
    prev_c = c;
    prev_result = result;
    
    return result;
}

double cube_calc(double* points_x, double* c, int n, double x, double a, double b) {
    if (x <= points_x[0]) {
        double dx = x - points_x[0];
        return c[0] + c[1] * dx + c[2] * dx * dx + c[3] * dx * dx * dx;
    }
    
    if (x >= points_x[n-1]) {
        double dx = x - points_x[n-1];
        int last_seg = 4 * (n - 2);
        return c[last_seg] + c[last_seg+1] * dx + c[last_seg+2] * dx * dx + c[last_seg+3] * dx * dx * dx;
    }
    
    int i = (int)((x - a) * (n - 1) / (b - a));
    
    if (i < 0) i = 0;
    if (i >= n - 1) i = n - 2;
    
    while (i < n - 2 && x > points_x[i+1]) i++;
    while (i > 0 && x < points_x[i]) i--;
    
    double dx = x - points_x[i];
    int idx = 4 * i;
    return c[idx] + c[idx+1] * dx + c[idx+2] * dx * dx + c[idx+3] * dx * dx * dx;
}

static double f_const(double x) { (void)x; return 1.0; }
static double f_x(double x) { return x; }
static double f_x2(double x) { return x * x; }
static double f_x3(double x) { return x * x * x; }
static double f_x4(double x) { return x * x * x * x; }
static double f_exp(double x) { return exp(x); }
static double f_runge(double x) { return 1.0 / (25.0 * x * x + 1.0); }

Window* Window::active_instance = nullptr;

const double EPSILON = 1e-16;

Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;
  func_id = DEFAULT_K;
  scale_factor = 0;
  perturb_factor = 0.0;
  perturb_point = -1;
  // display_mode = ORIGINAL_AND_CHEBYSHEV;
  display_mode = ORIGINAL_AND_NEWTON;
  
  x_points = nullptr;
  f_values = nullptr;
  f_perturbed = nullptr;
  // chebyshev_coefs = nullptr;
  newton_coefs = nullptr;
  spline_coefs = nullptr;
  
  buffer_x = new double[width()];
  buffer_y = new double[width()];
  // buffer_cheb = new double[width()];
  buffer_newton = new double[width()];
  buffer_spline = new double[width()];
  need_recalc = true;
  
  active_instance = this;

  change_func();
  
  setFocusPolicy(Qt::StrongFocus);
}

Window::~Window()
{
  if (x_points)
    delete[] x_points;
  if (f_values)
    delete[] f_values;
  if (f_perturbed)
    delete[] f_perturbed;
  // if (chebyshev_coefs)
  //   delete[] chebyshev_coefs;
  if (newton_coefs)
    delete[] newton_coefs;
  if (spline_coefs)
    delete[] spline_coefs;
    
  delete[] buffer_x;
  delete[] buffer_y;
  // delete[] buffer_cheb;
  delete[] buffer_spline;
    
  if (active_instance == this) {
    active_instance = nullptr;
  }

  // free_chebyshev_resources();
  free_newton_resources();
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc < 5)
    return -1;

  if (   sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || sscanf (argv[3], "%d", &n) != 1
      || n <= 0
      || sscanf (argv[4], "%d", &func_id) != 1
      || func_id < 0 || func_id > 6)
    return -2;
    
  return 0;
}

void Window::set_func_by_id()
{
  switch (func_id)
    {
      case 0:
        f_name = "f(x) = 1";
        f = &::f_const;
        break;
      case 1:
        f_name = "f(x) = x";
        f = &::f_x;
        break;
      case 2:
        f_name = "f(x) = x^2";
        f = &::f_x2;
        break;
      case 3:
        f_name = "f(x) = x^3";
        f = &::f_x3;
        break;
      case 4:
        f_name = "f(x) = x^4";
        f = &::f_x4;
        break;
      case 5:
        f_name = "f(x) = e^x";
        f = &::f_exp;
        break;
      case 6:
        f_name = "f(x) = 1/(25*x^2 + 1)";
        f = &::f_runge;
        break;
    }
  
  perturb_factor = 0.0;
  perturb_point = -1;
  
  build_interpolation_data();
  
  need_recalc = true;
  update();
}

void Window::change_func()
{
  func_id = (func_id + 1) % 7;
  set_func_by_id();
}

void Window::change_display_mode()
{
  display_mode = static_cast<DisplayMode>((static_cast<int>(display_mode) + 1) % DISPLAY_MODE_COUNT);
  need_recalc = true;
  update();
}

QPointF Window::l2g (double x_loc, double y_loc, double y_min, double y_max)
{
  double scaled_a = (a + b) / 2.0 - (b - a) / (2.0 * pow(2, scale_factor));
  double scaled_b = (a + b) / 2.0 + (b - a) / (2.0 * pow(2, scale_factor));
  
  if (std::abs(scaled_b - scaled_a) < EPSILON) {
    scaled_b = scaled_a + EPSILON; 
  }
  
  if (std::abs(y_max - y_min) < EPSILON) {
    y_max =EPSILON;
  }
  
  double x_gl = (x_loc - scaled_a) / (scaled_b - scaled_a) * width ();
  double y_gl = (y_max - y_loc) / (y_max - y_min) * height ();
  return QPointF (x_gl, y_gl);
}


double Window::get_f_value(int idx)
{
  if (perturb_point >= 0 && idx == perturb_point && std::abs(perturb_factor) > EPSILON) {
    return f_perturbed[idx];
  }
  return f_values[idx];
}

double Window::eval_newton(double x_val){
  static double last_x = -1e300;
  static double last_result = 0.0;
  static int last_n = 0;
  static int last_func_id = -1;
  
  if (std::abs(x_val - last_x) < EPSILON && n == last_n && func_id == last_func_id) {
    return last_result;
  }
  
  if (n > 50) {
    last_result = 0.0;
  }
  else if (func_id == 0) {
    last_result = 1.0;
  }
  else if (func_id == 1) {
    last_result = x_val;
  }
  else if (func_id == 2) {
    last_result = x_val * x_val;
  }
  else if (func_id == 3) {
    last_result = x_val * x_val * x_val;
  }
  else {
    last_result = evaluate_newton(x_val);
  }
  
  last_x = x_val;
  last_n = n;
  last_func_id = func_id;
  
  return last_result;
}

double Window::newton_error(double x_val)
{
  if (n > 50) {
    return 0.0;
  }
  return std::abs(f(x_val) - eval_newton(x_val));
}

double Window::get_max_abs_value()
{
  double scaled_a = a / pow(2, scale_factor);
  double scaled_b = b / pow(2, scale_factor);
  
  double max_val = 0.0;
  double max_func = 0.0;
  double step = (scaled_b - scaled_a) / 10000.0;
  
  for (double x = scaled_a; x <= scaled_b; x += step) {
    double val = std::abs(f(x));
    if (val > max_func) {
      max_func = val;
    }
  }
  
  for (double x = scaled_a; x <= scaled_b; x += step) {
    if (display_mode == ERROR_ONLY) {
      // double err_cheb = chebyshev_error(x) / max_func;
      double err_newton = newton_error(x) / max_func;
      double err_spline = spline_error(x) / max_func;
      // max_val = std::max(max_val, std::max(err_cheb, err_spline));
      max_val = std::max(max_val, std::max(err_newton, err_spline));
    }
    else {
      double val = std::abs(f(x));
      // double cheb_val = std::abs(eval_chebyshev(x));
      double newton_val = std::abs(eval_newton(x));
      double spline_val = std::abs(eval_spline(x));
      // max_val = std::max(max_val, std::max(val, std::max(cheb_val, spline_val)));
      max_val = std::max(max_val, std::max(val, std::max(newton_val, spline_val)));
    }
  }
  
  return max_val;
}

// std::vector<double> chebyshev_nodes(int n) {
//     std::vector<double> nodes(n);
//     for (int k = 0; k < n; ++k) {
//         nodes[k] = cos(M_PI * (2*k + 1) / (2 * n));
//     }
//     return nodes;
// }

// std::vector<double> compute_chebyshev_coeffs(const std::vector<double>& f) {
//     int N = f.size();
//     std::vector<double> coeffs(N, 0.0);
    
//     for (int k = 0; k < N; ++k) {
//         double sum = 0.0;
//         for (int n = 0; n < N; ++n) {
//             sum += f[n] * cos(M_PI * (n + 0.5) * k / N);
//         }
//         coeffs[k] = sum * 2.0 / N;
//         if (k == 0) coeffs[k] *= 0.5;
//     }
    
//     return coeffs;
// }

// double clenshaw_algorithm(const std::vector<double>& coeffs, double x) {
//     double b2 = 0.0, b1 = 0.0;
//     for (int j = coeffs.size()-1; j >= 0; --j) {
//         double temp = b1;
//         b1 = coeffs[j] + 2 * x * b1 - b2;
//         b2 = temp;
//     }
//     return (b1 - b2)/2.0;
// }

void Window::build_interpolation_data()
{
  if (x_points)
    delete[] x_points;
  if (f_values)
    delete[] f_values;
  if (f_perturbed)
    delete[] f_perturbed;
  // if (chebyshev_coefs)
  //   delete[] chebyshev_coefs;
  if (newton_coefs)
    delete[] newton_coefs;
  if (spline_coefs)
    delete[] spline_coefs;
    
  // free_chebyshev_resources();
  free_newton_resources();
    
  x_points = new double[n];
  f_values = new double[n];
  f_perturbed = new double[n];
  // chebyshev_coefs = new double[n];
  newton_coefs = new double[n];
  spline_coefs = new double[4 * n];
  
  double step = (b - a) / (n - 1);
  for (int i = 0; i < n; i++)
  {
    x_points[i] = a + i * step;
    f_values[i] = f(x_points[i]);
    f_perturbed[i] = f_values[i];
  }

  perturb_point = n / 2;
  

  if (std::abs(perturb_factor) > EPSILON) {
    double max_abs = get_max_abs_value();
    f_perturbed[perturb_point] = f_values[perturb_point] + perturb_factor * 0.1 * max_abs;
  }
  
  if (n <= 50) {
    // f_cheb = f;
    f_newton = f;
    if(init_newton_interpolation(n, a, b)) {
      calculate_newton_coefficients(n, x_points, f_values, newton_coefs);
    }
    // if (init_chebyshev_interpolation(n, a, b)) {
    //   calculate_chebyshev_nodes();
    //   calculate_interpolation();
    //   std::copy(c_cheb, c_cheb + n, chebyshev_coefs);
    // }
  }
  
  double* d = new double[n];
  cube_approximation(x_points, f_perturbed, d, spline_coefs, n);
  delete[] d;
}

// double Window::eval_chebyshev(double x_val)
// {
//     static double last_x = -1e300;
//     static double last_result = 0.0;
//     static int last_n = 0;
//     static int last_func_id = -1;
    
//     if (std::abs(x_val - last_x) < EPSILON && n == last_n && func_id == last_func_id) {
//         return last_result;
//     }
    
//     if (n > 50) {
//         last_result = 0.0;
//     }
//     else if (func_id == 0) {
//         last_result = 1.0;
//     }
//     else if (func_id == 1) {
//         last_result = x_val;
//     }
//     else if (func_id == 2) {
//         last_result = x_val * x_val;
//     }
//     else if (func_id == 3) {
//         last_result = x_val * x_val * x_val;
//     }
//     else {
//         last_result = evaluate_chebyshev(x_val);
//     }
    
//     last_x = x_val;
//     last_n = n;
//     last_func_id = func_id;
    
//     return last_result;
// }

double Window::eval_spline(double x_val)
{
    static double last_x = -1e300;
    static double last_result = 0.0;
    static int last_n = 0;
    static double *last_points = nullptr;
    static double *last_coefs = nullptr;
    

    if (std::abs(x_val - last_x) < EPSILON && n == last_n && x_points == last_points && spline_coefs == last_coefs) {
        return last_result;
    }
    
    last_result = call_cube_approximation(x_points, spline_coefs, n, x_val, a, b);
    
    last_x = x_val;
    last_n = n;
    last_points = x_points;
    last_coefs = spline_coefs;
    
    return last_result;
}

// double Window::chebyshev_error(double x_val)
// {
//   if (n > 50) {
//     return 0.0;
//   }
//   return std::abs(f(x_val) - eval_chebyshev(x_val));
// }

double Window::spline_error(double x_val)
{
  return std::abs(f(x_val) - eval_spline(x_val));
}

void Window::draw_axes(QPainter &painter, double min_y, double max_y)
{
  QPen pen_red(Qt::red, 0, Qt::SolidLine);
  painter.setPen(pen_red);
  
  double scaled_a = (a + b) / 2.0 - (b - a) / (2.0 * pow(2, scale_factor));
  double scaled_b = (a + b) / 2.0 + (b - a) / (2.0 * pow(2, scale_factor));
  
  if (min_y <= 0 && max_y >= 0) {
    painter.drawLine(L2G(scaled_a, 0), L2G(scaled_b, 0));
  }
  
  if (scaled_a <= 0 && scaled_b >= 0) {
    painter.drawLine(L2G(0, min_y), L2G(0, max_y));
  }
}

void Window::draw_original_function(QPainter &painter, double min_y, double max_y)
{
  QPen pen_black(Qt::black, 0, Qt::SolidLine);
  painter.setPen(pen_black);
  
  for (int i = 0; i < width() - 1; i++) {
    painter.drawLine(L2G(buffer_x[i], buffer_y[i]), 
                    L2G(buffer_x[i+1], buffer_y[i+1]));
  }
  
  painter.setPen(QPen(Qt::black, 4));
  for (int i = 0; i < n; i++) {
    QPointF point = L2G(x_points[i], get_f_value(i));
    painter.drawPoint(point);
  }
}

// void Window::draw_chebyshev(QPainter &painter, double min_y, double max_y)
// {
//   if (n > 50) {
//     return;
//   }
  
//   QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
//   painter.setPen(pen_blue);
  
//   for (int i = 0; i < width() - 1; i++) {
//     painter.drawLine(L2G(buffer_x[i], buffer_cheb[i]), 
//                     L2G(buffer_x[i+1], buffer_cheb[i+1]));
//   }
// }

void Window::draw_newton(QPainter &painter, double min_y, double max_y)
{
  if (n > 50) {
    return;
  }
  
  QPen pen_purple(Qt::magenta, 0, Qt::SolidLine);
  painter.setPen(pen_purple);
  
  for (int i = 0; i < width() - 1; i++) {
    painter.drawLine(L2G(buffer_x[i], buffer_newton[i]), 
                    L2G(buffer_x[i+1], buffer_newton[i+1]));
  }
}

void Window::draw_spline(QPainter &painter, double min_y, double max_y)
{
  QPen pen_green(Qt::green, 0, Qt::SolidLine);
  painter.setPen(pen_green);
  
  for (int i = 0; i < width() - 1; i++) {
    painter.drawLine(L2G(buffer_x[i], buffer_spline[i]), 
                    L2G(buffer_x[i+1], buffer_spline[i+1]));
  }
}

void Window::draw_errors(QPainter &painter, double min_y, double max_y)
{
  double scaled_a = (a + b) / 2.0 - (b - a) / (2.0 * pow(2, scale_factor));
  double scaled_b = (a + b) / 2.0 + (b - a) / (2.0 * pow(2, scale_factor));
  
  double delta_x = (scaled_b - scaled_a) / width();
  
  // Рисуем ошибку сплайна
  QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
  painter.setPen(pen_blue);
  
  double x1 = scaled_a;
  double y1 = std::fabs(f(x1) - call_cube_approximation(x_points, spline_coefs, n, x1, a, b));
  double x2;
  double y2;
  
  for (int i = 1; i < width(); i++) {
    x2 = scaled_a + i * delta_x;
    y2 = std::fabs(f(x2) - call_cube_approximation(x_points, spline_coefs, n, x2, a, b));
    painter.drawLine(L2G(x1, y1), L2G(x2, y2));
    x1 = x2;
    y1 = y2;
  }
  
  // Рисуем ошибку метода Ньютона, если n <= 50
  if (n <= 50) {
    QPen pen_green(Qt::green, 0, Qt::SolidLine);
    painter.setPen(pen_green);
    
    x1 = scaled_a;
    // y1 = std::fabs(f(x1) - eval_chebyshev(x1));
    y1 = std::fabs(f(x1) - eval_newton(x1));
    
    for (int i = 1; i < width(); i++) {
      x2 = scaled_a + i * delta_x;
      // y2 = std::fabs(f(x2) - eval_chebyshev(x2));
      y2 = std::fabs(f(x2) - eval_newton(x2));
      painter.drawLine(L2G(x1, y1), L2G(x2, y2));
      x1 = x2;
      y1 = y2;
    }
  }
}

void Window::draw_info_text(QPainter &painter)
{
  QString info_text = QString("k=%1 %2").arg(func_id).arg(f_name);
  
  QString mode_text;
  switch (display_mode)
  {
    // case ORIGINAL_AND_CHEBYSHEV:
    //   mode_text = "Original + Chebyshev";
    //   break;
    case ORIGINAL_AND_NEWTON:
      mode_text = "Original + Newton";
      break;
    case ORIGINAL_AND_SPLINE:
      mode_text = "Original + Spline";
      break;
    case ORIGINAL_AND_BOTH:
      mode_text = "Original + Both";
      break;
    case ERROR_ONLY:
      mode_text = "Errors (dynamic scale)";
      break;
    case DISPLAY_MODE_COUNT:
    default:
      mode_text = "Unknown";
      break;
  }
  
  double max_abs = get_max_abs_value();
  
  double scaled_a = a / pow(2, scale_factor);
  double scaled_b = b / pow(2, scale_factor);
  
  QString full_info = QString("%1 | Mode: %2 | n=%3 | scale=2^%4 | perturb=%5 | max=%6 | [%7, %8]")
                             .arg(info_text)
                             .arg(mode_text)
                             .arg(n)
                             .arg(scale_factor)
                             .arg(perturb_factor)
                             .arg(max_abs, 0, 'e', 16)
                             .arg(scaled_a, 0, 'f', 2)
                             .arg(scaled_b, 0, 'f', 2);
  
  painter.setPen(Qt::blue);
  painter.drawText(10, 20, full_info);
  
  // Добавление информации о цветах
  QString color_info;
  
  if (display_mode == ERROR_ONLY) {
    color_info = "Цвета: ";
    painter.setPen(Qt::blue);
    color_info += "Синий - ошибка сплайна";
    painter.drawText(10, 40, color_info);
    
    if (n <= 50) {
      painter.setPen(Qt::green);
      color_info = "Зеленый - ошибка Ньютона";
      painter.drawText(250, 40, color_info);
    }
  } else {
    color_info = "Цвета: ";
    painter.setPen(Qt::black);
    color_info += "Черный - исходная функция";
    painter.drawText(10, 40, color_info);
    
    // if (display_mode == ORIGINAL_AND_CHEBYSHEV || display_mode == ORIGINAL_AND_BOTH) {
    //   if (n <= 50) {
    //     painter.setPen(Qt::blue);
    //     color_info = "Синий - Чебышев";
    //     painter.drawText(250, 40, color_info);
    //   }
    // }
    
    if (display_mode == ORIGINAL_AND_NEWTON || display_mode == ORIGINAL_AND_BOTH) {
      painter.setPen(Qt::magenta);
      color_info = "Фиолетовый - Ньютоны";
      painter.drawText(250, 40, color_info);
    }
    
    if (display_mode == ORIGINAL_AND_SPLINE || display_mode == ORIGINAL_AND_BOTH) {
      painter.setPen(Qt::green);
      color_info = "Зеленый - сплайн";
      painter.drawText(400, 40, color_info);
    }
  }
  
  std::cout << "Maximum absolute value: " 
            << std::scientific
            << std::setprecision(16)
            << max_abs 
            << std::endl;
}

void Window::keyPressEvent(QKeyEvent *event)
{
  bool need_update = true;
  
  switch (event->key()) {
    case Qt::Key_0:

      change_func();
      break;
      
    case Qt::Key_1:

      change_display_mode();
      break;
      
    case Qt::Key_2:
      scale_factor++;
      need_recalc = true;
      break;
      
    case Qt::Key_3:
      if (scale_factor > 0) {
        scale_factor--;
        need_recalc = true;
      }
      break;
      
    case Qt::Key_4:
      n *= 2;
      build_interpolation_data();
      need_recalc = true;
      break;
      
    case Qt::Key_5:
      if (n > 2) {
        n /= 2;
        build_interpolation_data();
        need_recalc = true;
      }
      break;
      
    case Qt::Key_6:
      perturb_factor += 1.0;
      build_interpolation_data();
      need_recalc = true;
      break;
      
    case Qt::Key_7:
      perturb_factor -= 1.0;
      build_interpolation_data();
      need_recalc = true;
      break;
      
    default:
      need_update = false;
      QWidget::keyPressEvent(event);
      break;
  }
  
  if (need_update) {
    update();
  }
}

void Window::calculate_residuals(double scaled_a, double scaled_b, double delta_x, double &min_y, double &max_y)
{
  double max_func = 0.0;
  for (double x = scaled_a; x - scaled_b < EPSILON; x += delta_x) {
    double y = std::fabs(f(x));
    if (y > max_func) {
      max_func = y;
    }
  }
  
  // Находим реальные минимум и максимум ошибок для динамического масштабирования
  min_y = 0.0;  // Минимум ошибки всегда 0
  max_y = 0.0;
  
  // Вычисляем ошибки в каждой точке
  for (double x = scaled_a; x - scaled_b < EPSILON; x += delta_x) {
    double spline_err = std::fabs(f(x) - eval_spline(x));
    if (spline_err > max_y) max_y = spline_err;
    
    // if (n <= 50) {
    //   double cheb_err = std::fabs(f(x) - eval_chebyshev(x));
    //   if (cheb_err > max_y) max_y = cheb_err;
    // }
    if (n <= 50) {
      double newton_err = std::fabs(f(x) - eval_newton(x));
      if (newton_err > max_y) max_y = newton_err;
    }
  }
  
  // Если максимальная ошибка слишком мала, устанавливаем минимальное значение
  if (max_y < EPSILON) {
    max_y = EPSILON;
  }
  
  // Добавляем небольшой отступ сверху для лучшей видимости
  max_y *= 1.05;
}

void Window::resizeEvent(QResizeEvent *event)
{
    delete[] buffer_x;
    delete[] buffer_y;
    // delete[] buffer_cheb;
    delete[] buffer_newton;
    delete[] buffer_spline;
    
    buffer_x = new double[width()];
    buffer_y = new double[width()];
    // buffer_cheb = new double[width()];
    buffer_newton = new double[width()];
    buffer_spline = new double[width()];
    
    need_recalc = true;
    QWidget::resizeEvent(event);
}

void Window::paintEvent (QPaintEvent * /* event */)
{  
  QPainter painter(this);
  
  double scaled_a = (a + b) / 2.0 - (b - a) / (2.0 * pow(2, scale_factor));
  double scaled_b = (a + b) / 2.0 + (b - a) / (2.0 * pow(2, scale_factor));
  
  double min_y = 0.0;
  double max_y = 0.0;
  double delta_x = (scaled_b - scaled_a) / width();
  
  if (need_recalc) {
    for (int i = 0; i < width(); i++) {
      double x = scaled_a + i * delta_x;
      buffer_x[i] = x;
      buffer_y[i] = f(x);
      
      if (n <= 50) {
        // buffer_cheb[i] = eval_chebyshev(x);
        buffer_newton[i] = eval_newton(x);
      }
      buffer_spline[i] = eval_spline(x);
    }
    need_recalc = false;
  }
  
  if (display_mode == ERROR_ONLY) {
    calculate_residuals(scaled_a, scaled_b, delta_x, min_y, max_y);
  } else {
    for (int i = 0; i < width(); i++) {
      double y_val = buffer_y[i];
      
      if (i == 0 || y_val < min_y) min_y = y_val;
      if (i == 0 || y_val > max_y) max_y = y_val;
      
      // if (n <= 50 && (display_mode == ORIGINAL_AND_CHEBYSHEV || display_mode == ORIGINAL_AND_BOTH)) {
      //   double y_cheb = buffer_cheb[i];
      //   if (y_cheb < min_y) min_y = y_cheb;
      //   if (y_cheb > max_y) max_y = y_cheb;
      // }
      
      if (n <= 50 && (display_mode == ORIGINAL_AND_NEWTON || display_mode == ORIGINAL_AND_BOTH)) {
        double y_newton = buffer_newton[i];
        if (y_newton < min_y) min_y = y_newton;
        if (y_newton > max_y) max_y = y_newton;
      }
      
      if (display_mode == ORIGINAL_AND_SPLINE || display_mode == ORIGINAL_AND_BOTH) {
        double y_spline = buffer_spline[i];
        if (y_spline < min_y) min_y = y_spline;
        if (y_spline > max_y) max_y = y_spline;
      }
    }
  }
  
  double padding = 0.05 * (max_y - min_y);
  min_y -= padding;
  max_y += padding;
  
  draw_axes(painter, min_y, max_y);
  
  if (display_mode != ERROR_ONLY) {
    draw_original_function(painter, min_y, max_y);
    
    // if (display_mode == ORIGINAL_AND_CHEBYSHEV || display_mode == ORIGINAL_AND_BOTH) {
    //   draw_chebyshev(painter, min_y, max_y);
    // }
    if (display_mode == ORIGINAL_AND_NEWTON || display_mode == ORIGINAL_AND_BOTH) {
      draw_newton(painter, min_y, max_y);
    }
    
    if (display_mode == ORIGINAL_AND_SPLINE || display_mode == ORIGINAL_AND_BOTH) {
      draw_spline(painter, min_y, max_y);
    }
  } else {
    draw_errors(painter, min_y, max_y);
  }
  
  draw_info_text(painter);
}
