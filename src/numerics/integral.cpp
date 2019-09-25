#include "numerics/integral.h"
/*
// diese Funktion wird nicht verwendet, stattdessen die Version in integral.h
double Integral(const function<double(double)> &f, double x0, double xend,
    int iter)
//Simpsonsumme zur Berechnung des Integrals der Funktion f(x) von x0 bis xend
//Anzahl der Intervalle: iter
    {
  double erg;
  double l = (xend - x0);
  double h = l / iter;
  erg = f(x0) - f(xend);
  for (int i = 0; i < iter; i++) {
    erg = erg + 4 * f(x0 + i * h + h / 2) + 2 * f(x0 +  (i + 1) * h);
  }
  erg = h / 6 * erg;
  return erg;

}
*/
