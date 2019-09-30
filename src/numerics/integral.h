#pragma once

#include <functional>
#include <iostream>
#include <cmath>

using namespace std;

//! adaptive Archimedes quadrature in 1D, this method only computes surplusses ∫_a^b f(s)ds, i.e. f(a)=f(b)=0, fa=f(a), fb=f(b)
template<typename Fn, Fn *f, typename... Args>
double archimedes1d_adaptive_surplus(double a, double b, double fa, double fb, int level, Args&&... args)
{
  double x_middle = 0.5*(a + b);             // middle point of the integration interval
  double f_middle = f(x_middle, std::forward<Args>(args)...);            // evaluate integrand at center of interval
  double height = f_middle - 0.5*(fa + fb);  // height of triangle

  double height_ba = height * (b-a);
  
  const double epsilon = 1e-5;   // error tolerance

  if (level >= 3 || fabs(height_ba) < epsilon)    // recursion end criterion
  {
    return 0.5 * 4./3. * height_ba;
  }
  else
  {
    return 0.5 * height_ba
      + archimedes1d_adaptive_surplus<Fn,f,Args...>(a,        x_middle, fa,       f_middle, level+1, std::forward<Args>(args)...)
      + archimedes1d_adaptive_surplus<Fn,f,Args...>(x_middle, b,        f_middle, fb,       level+1, std::forward<Args>(args)...);
  }
};

// source: https://stackoverflow.com/a/25403872/10290071
/** Compute numerical integral of f over the domain [x0,xend]
 * @param x0 begin of integration interval
 * @param xend end of integration interval
 * @param fa   integrand evaluated at x0, sometimes this is known and can be hardcoded to avoid function calls
 */
template<typename Fn, Fn *f, typename... Args>
double Integral_(double x0, double xend, double fa, Args&&... args)
{
#if 0
  // normal simpson rule
  double h = xend - x0
  return h/6. * (
    f(x0, std::forward<Args>(args)...)
    + 4*f(x0 + 0.5*h, std::forward<Args>(args)...)
    + f(xend, std::forward<Args>(args)...)
  );
#endif

#if 0
  // composite simpson rule
  const int iter = 14;   // number of sampling points
  double erg;
  double l = (xend - x0);
  double h = l / iter;
  erg = f(x0, std::forward<Args>(args)...) - f(xend, std::forward<Args>(args)...);
  for (int i = 0; i < iter; i++) {
    erg +=   4 * f(x0 + i * h + h / 2, std::forward<Args>(args)...)
           + 2 * f(x0 + (i + 1) * h, std::forward<Args>(args)...);
  }
  erg = h / 6 * erg;
  return erg;
#endif

#if 1
  // adaptive archimedes quadrature
  const double a = x0;
  const double b = xend;
  const double fb = f(b, std::forward<Args>(args)...);

  const double trapezoid = 0.5 * (b-a) * (fa + fb);

  const double result = trapezoid + archimedes1d_adaptive_surplus<Fn,f,Args...>(a, b, fa, fb, 0, std::forward<Args>(args)...);
  return result;
#endif
}

#define Integral(FUNC, ...) Integral_<decltype(FUNC), FUNC>(__VA_ARGS__)
// usage:
// example: function to integrate is cos(a*x + b) with parameters a and b
//
// double function(double x, double a, double b) {
//   return cos(a*x + b);
// }
// ...
// Integral(function, x0, xend, function(x,a,b), a, b);
