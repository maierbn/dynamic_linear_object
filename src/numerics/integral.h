#pragma once

#include <functional>

using namespace std;

//Simpsonsumme zur Berechnung des Integrals der Funktion f(x) von x0 bis xend
//Anzahl der Intervalle: iter
double Integral(const function<double(double)> &f, double x0, double xend, int iter);
/*
// source: https://stackoverflow.com/a/25403872/10290071
template<typename Fn, Fn fn, typename... Args>
typename std::result_of<Fn(Args...)>::type
Integral_(double x0, double xend, int iter, Args&&... args)
{
  double erg;
  double l = (xend - x0);
  double h = l / iter;
  erg = f(x0) - f(xend);
  for (int i = 0; i < iter; i++) {
    erg = erg + 4 * f(x0 + i * h + h / 2, std::forward<Args>(args)...) + 2 * f(x0 + (i + 1) * h, std::forward<Args>(args)...);
  }
  erg = h / 6 * erg;
  return erg;
}

#define Integral(FUNC, ...) Integral_<decltype(&FUNC), &FUNC>(__VA_ARGS__)
// Verwendung:
// Beispiel: zu integrierende Funktion cos(a*x + b) mit Parametern a und b
// Integral([](double t, double a, double b), x0, xend, iter, a, b);
*/
