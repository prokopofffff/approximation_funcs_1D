#include "interpolation.h"
#include <cmath>
#include <cstring>

// Chebyshev interpolation implementation
// void build_chebyshev_interpolation(int n, double *x, double *f, double *a)
// {
//     bool is_constant = true;
//     for(int i = 1; i < n; i++) {
//         if(std::abs(f[i] - f[0]) > 1e-10) {
//             is_constant = false;
//             break;
//         }
//     }
    
//     if(is_constant) {
//         a[0] = f[0];
//         for(int i = 1; i < n; i++) {
//             a[i] = 0.0;
//         }
//         return;
//     }
    
//     (void)x;
    
//     double *chebyshev_nodes = new double[n];
//     for (int i = 0; i < n; i++)
//     {
//         chebyshev_nodes[i] = cos((2.0 * i + 1.0) * M_PI / (2.0 * n));
//     }

//     for (int i = 0; i < n; i++)
//     {
//         a[i] = 0.0;
//         for (int j = 0; j < n; j++)
//         {
//             double T_ij = 1.0;
//             if (i == 1)
//                 T_ij = chebyshev_nodes[j];
//             else if (i > 1)
//             {
//                 double T_i_minus_2 = 1.0;
//                 double T_i_minus_1 = chebyshev_nodes[j];
//                 for (int k = 2; k <= i; k++)
//                 {
//                     T_ij = 2.0 * chebyshev_nodes[j] * T_i_minus_1 - T_i_minus_2;
//                     T_i_minus_2 = T_i_minus_1;
//                     T_i_minus_1 = T_ij;
//                 }
//             }
//             a[i] += f[j] * T_ij;
//         }
//         a[i] *= (2.0 / n);
//     }
//     a[0] /= 2.0;

//     delete[] chebyshev_nodes;
// }

// double eval_chebyshev_interpolation(double x_val, double a, double b, int n, double *x, double *coef)
// {
//     (void)x;
    
//     double t = (2.0 * x_val - a - b) / (b - a);
    
//     double b_n = 0.0;
//     double b_n_1 = 0.0;
//     double b_n_2;
    
//     for (int j = n - 1; j >= 0; j--)
//     {
//         b_n_2 = b_n_1;
//         b_n_1 = b_n;
//         b_n = 2.0 * t * b_n_1 - b_n_2 + coef[j];
//     }
    
//     return b_n - t * b_n_1;
// }

// Вычисление коэффициентов Ньютона (разделённых разностей)
void calculate_newton_coefficients(int n, double *x, double *f_x, double *c) {
    for (int i = 0; i < n; ++i)
        c[i] = f_x[i];
    for (int k = 1; k < n; ++k) {
        for (int i = n - 1; i >= k; --i) {
            c[i] = (c[i] - c[i-1]) / (x[i] - x[i-k]);
        }
    }
}

// Вычисление значения многочлена Ньютона в точке t
// n - число узлов, x - массив узлов, c - коэффициенты Ньютона
// Возвращает значение P(t)
double evaluate_newton_polynomial(double t, int n, double *x, double *c) {
    double result = c[n-1];
    for (int i = n-2; i >= 0; --i) {
        result = result * (t - x[i]) + c[i];
    }
    return result;
}

void build_cubic_spline(int n, double *x, double *f, double *a, double (*func)(double))
{
    // Проверка на константную функцию
    bool is_constant = true;
    for(int i = 1; i < n; i++) {
        if(std::abs(f[i] - f[0]) > 1e-10) {
            is_constant = false;
            break;
        }
    }
    
    if(is_constant) {
        // Для константной функции устанавливаем только коэффициенты a₀
        for(int i = 0; i < n-1; i++) {
            a[4*i] = f[0];      // a₀ = значение функции
            a[4*i+1] = 0.0;     // a₁ = 0 (первая производная)
            a[4*i+2] = 0.0;     // a₂ = 0 (вторая производная / 2)
            a[4*i+3] = 0.0;     // a₃ = 0 (третья производная / 6)
        }
        return;
    }
    
    double *h = new double[n - 1]; 
    double *alpha = new double[n - 1];
    double *l = new double[n];
    double *mu = new double[n - 1];
    double *z = new double[n];
    double *c = new double[n]; 
    
    for (int i = 0; i < n - 1; i++)
    {
        h[i] = x[i + 1] - x[i];
    }
    
    double h_small = 1e-16;
    double f_minus = func(x[0] - h_small);
    double f_center = f[0];
    double f_plus = func(x[0] + h_small);
    c[0] = (f_minus - 2 * f_center + f_plus) / (h_small * h_small);
    
    f_minus = func(x[n-1] - h_small);
    f_center = f[n-1];
    f_plus = func(x[n-1] + h_small);
    c[n-1] = (f_minus - 2 * f_center + f_plus) / (h_small * h_small);
    
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = c[0];
    
    for (int i = 1; i < n - 1; i++)
    {
        alpha[i] = 3.0 * ((f[i+1] - f[i]) / h[i] - (f[i] - f[i-1]) / h[i-1]);
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }
    
    l[n-1] = 1.0;
    z[n-1] = c[n-1];
    
    for (int j = n - 2; j >= 0; j--)
    {
        c[j] = z[j] - mu[j] * c[j+1];
    }
    
    for (int i = 0; i < n - 1; i++)
    {
        double b = (f[i+1] - f[i]) / h[i] - h[i] * (c[i+1] + 2.0 * c[i]) / 3.0;
        double d = (c[i+1] - c[i]) / (3.0 * h[i]);
        
        a[4*i] = f[i];
        a[4*i+1] = b;
        a[4*i+2] = c[i];
        a[4*i+3] = d;
    }
    
    // Clean up
    delete[] h;
    delete[] alpha;
    delete[] l;
    delete[] mu;
    delete[] z;
    delete[] c;
}

double eval_cubic_spline(double x_val, double a, double b, int n, double *x, double *coef)
{
    (void)a;
    (void)b;
    
    int segment = 0;
    while (segment < n - 2 && x_val > x[segment + 1])
    {
        segment++;
    }
    
    double dx = x_val - x[segment];
    double result = coef[4*segment] + coef[4*segment+1] * dx + 
                   coef[4*segment+2] * dx * dx + 
                   coef[4*segment+3] * dx * dx * dx;
    
    return result;
} 