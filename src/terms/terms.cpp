#include "terms/terms.h"

#include "terms/globals.h"
#include "utility/operators.h"
#include "numerics/integral.h"
#include "terms/problem_definition.h"

#include <cmath>

double N(int i, double h, double s) {
  //Hütchenfunktion bzw lineare Basisfunktion -> 1 für i*h, linear abfallend auf 0 zu (i-1)*h, (i+1)*h, ansonsten konstant 0
  double erg = 0;
  if (s <= (i - 1) * h || s >= (i + 1) * h) {
    erg = 0;
  } else {
    if (s < i * h) {
      erg = (s - (i - 1) * h) / h;
    } else {
      erg = ((i + 1) * h - s) / h;
    }
  }
  return erg;
}
double Theta(const vector<double> &thetaN, double h, double s) {
  //Theta gibt anhand Linearkombination der N+1 Koeffizienten mit den jeweiligen Basisfunktionen
  //den Wert von Theta an der Stelle s zurück
  double erg = 0;
  size_t i1 = s / h;
  size_t i2 = i1+1;
  if (i1 >= 0 && i1 <= thetaN.size() - 1) {
    double erg2 = N(i1, h, s);
    erg = erg + thetaN[i1] * erg2;
  }
  if (i2 >= 0 && i2 <= thetaN.size() - 1 && i1 != i2) {
    double erg3 = N(i2, h, s);
    erg = erg + thetaN[i2] * erg3;
  }
  return erg;

}

double CosIntIntegrand(double x, int i, double h, const vector<double> &thetaN) {
  //Berechnet den Integranden für CosInt
  return cos(Theta(thetaN,h,x))*N(i,h,x);
}

double CosInt(int i, double s, double h, const vector<double> &thetaN,
    int iter) {
  //Berechnet das Integral über den Cosinus von Theta multipliziert mit der i-ten Basisfunktion
  double x1 = max(i - 1.0, 0.0) * h;
  double x2 = min((i + 1.0) * h, s);

  return Integral(CosIntIntegrand, x1, x2, iter,
                  i, h, thetaN);

  /*
   return Integral([i,h,thetaN](double x) {
   return cos(Theta(thetaN,h,x))*N(i,h,x);
   }, 0, s, iter);*/
}

double SinIntIntegrand(double x, int i, double h, const vector<double> &thetaN) {
  //Berechnet den Integranden für SinInt
  return sin(Theta(thetaN,h,x))*N(i,h,x);
}

double SinInt(int i, double s, double h, const vector<double> &thetaN,
    int iter) {
  //Berechnet das Integral über den Sinus von Theta multipliziert mit der i-ten Basisfunktion

  double x1 = max(i - 1.0, 0.0) * h;
  double x2 = min((i + 1.0) * h, s);
  return Integral(SinIntIntegrand, x1, x2, iter,
                  i, h, thetaN);
  /*
   return Integral([i,h,thetaN](double x) {
   return sin(Theta(thetaN,h,x))*N(i,h,x);
   }, 0, s, iter);*/
}

double mIntegrand(double x, int i, int k, double h, int iter, const vector<double> &thetaN)
{
  return SinInt( i,x,h,thetaN,iter)*SinInt( k,x,h,thetaN,iter)+CosInt( i,x,h,thetaN,iter)*CosInt( k,x,h,thetaN,iter);
}

double m(int i, int k, double rho, double h, const vector<double> &thetaN,
    double l, int iter) {
  //Berechnet den Wert m_(i,k) aus Gleichung 14
  return rho
      * Integral(mIntegrand, 0, l, iter,
                 i,k,h,iter,thetaN);
}

double dSinIntIntegrand(double x, int i, int r, double h, const vector<double> &thetaN)
{
  return cos(Theta(thetaN,h,x))*N(i,h,x)*N(r,h,x);
}

double dSinInt(int i, int r, double s, double h, const vector<double> &thetaN,
    int iter) {
  //Berechnet die Ableitung der Funktion SInInt  nach thetar

  if (abs(r - i) <= 1 && s > (r - 1) * h && s > (i - 1) * h) {
    double x1 = max(min(r, i) - 1.0, 0.0) * h;
    double x2 = min((max(r, i) + 1.0) * h, s);
    return Integral(dSinIntIntegrand, x1, x2, iter,
                    i,r,h,thetaN);
  } else {
    return 0.0;
  }
  /*
   return Integral([i,r,h,thetaN](double x) {
   return cos(Theta(thetaN,h,x))*N(i,h,x)*N(r,h,x);
   }, 0, s, iter);*/
}

double dCosIntIntegrand(double x, int i, int r, double h, const vector<double> &thetaN)
{
  return -sin(Theta(thetaN,h,x))*N(i,h,x)*N(r,h,x);
}

double dCosInt(int i, int r, double s, double h, const vector<double> &thetaN,
    int iter) {
  //Berechnet die Ableitung der Funktion CosInt nach thetar

  if (abs(r - i) <= 1 && s > (r - 1) * h && s > (i - 1) * h) {
    double x1 = max(min(r, i) - 1.0, 0.0) * h;
    double x2 = min((max(r, i) + 1.0) * h, s);

    return Integral(dCosIntIntegrand, x1, x2, iter,
                    i,r,h,thetaN);

  } else {
    return 0.0;
  }/*
   return Integral([i,r,h,thetaN](double x) {
   return -sin(Theta(thetaN,h,x))*N(i,h,x)*N(r,h,x);
   }, 0, s, iter);*/
}

double AIntegrand(double x, int i, int r, int k, double h, vector<double> &thetaN, int iter)
{
  return dSinInt(i,r,x,h,thetaN,iter)*SinInt(k,x,h,thetaN,iter);
}

vector<vector<vector<double>>> A(double l, double h, int n,
    vector<double> &thetaN, int iter) {
  //Berechnet für alle i,k,r das Integral über dSinInt(i,r)*SinInt(k) zur Wiederverwendung
  vector<vector<vector<double>>> out;
  out.resize(n + 1);
  for (int i = 0; i <= n; i++) {
    out[i].resize(n + 1);
    for (int k = 0; k <= n; k++) {
      out[i][k].resize(n + 1);
    }
  }
  for (int i = 0; i <= n; i++) {
    for (int k = 0; k <= n; k++) {
      for (int r = 0; r <= n; r++) {
        if (r < i) {
          out[i][k][r] = out[r][k][i];
        } else {
          out[i][k][r] =
              Integral(AIntegrand, 0, l, iter,
                       i,r,k,h,thetaN,iter);
        }
      }
    }
  }
  return out;
}

double BIntegrand(double x, int i, int r, int k, double h, vector<double> &thetaN, int iter)
{
  return dCosInt(i,r,x,h,thetaN,iter)*CosInt(k,x,h,thetaN,iter);
}

vector<vector<vector<double>>> B(double l, double h, int n,
    vector<double> &thetaN, int iter) {
  //Berechnet für alle i,k,r das Integral über dCosInt(i,r)*CosInt(k) zur Wiederverwendung
  vector<vector<vector<double>>> out;
  out.resize(n + 1);
  for (int i = 0; i <= n; i++) {
    out[i].resize(n + 1);
    for (int k = 0; k <= n; k++) {
      out[i][k].resize(n + 1);
    }
  }
  for (int i = 0; i <= n; i++) {
    for (int k = 0; k <= n; k++) {
      for (int r = 0; r <= n; r++) {
        if (r < i) {
          out[i][k][r] = out[r][k][i];
        } else {
          out[i][k][r] =
              Integral(BIntegrand, 0, l, iter,
                       i,r,k,h,thetaN,iter);
        }
      }
    }
  }
  return out;
}

vector<vector<double>> K(int n, double RFlex, double h)
//Berechnet die Steifigkeitsmatrix des Kabels mit N Stützstellen mit Abstand h und der Steifigkeit RFlex
    {
  vector<vector<double>> out;
  out.resize(n + 1);
  for (int i = 0; i <= n; i++) {
    out[i].resize(n + 1);
  }

  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      if (i == j && i > 0 && i < n) {
        out[i][j] = 2 * RFlex / h;
      } else if (i == j) {
        out[i][j] = RFlex / h;

      } else if (j - i == 1 || i - j == 1) {
        out[i][j] = -RFlex / h;

      }
    }
  }
  return out;
}

double dm(int i, int k, int r, double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B) {
  //Berechnet die Ableitung der Funktion m(i,k) nach Theta_r
  return rho * (A[i][k][r] + A[k][i][r] + B[i][k][r] + B[k][i][r]);
}

double Z(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B) {
  //Berechnet den Wer Z aus Gleichung 16
  int imax = thetaN.size();
  double temp = 0;
  for (int i = 0; i < imax; i++) {
    temp += 0.5 * (dm(i, k, r, rho, A, B) - dm(r, i, k, rho, A, B))
        * thetaNp[i];

  }
  return temp;
}

double Y(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B) {
  //Berechnet den Wer Y aus Gleichung 16
  int imax = thetaN.size();
  double temp = 0;
  for (int i = 0; i < imax; i++) {
    temp -= 0.5 * dm(r, i, k, rho, A, B) * thetaNp[i];

  }
  return temp;

}
vector<vector<double>> mM(double rho, double h, vector<double> &thetaN,
    double l, int iter) {
  //Berechnet die Matrix M, mit Mij = m(i,j,...)
  vector<vector<double>> out;

  int imax = thetaN.size();
  out.resize(imax);
  for (int i = 0; i < imax; i++) {
    out[i].resize(imax);
  }

  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        out[i][j] = out[j][i];
      } else {
        out[i][j] = m(i, j, rho, h, thetaN, l, iter);
      }
    }
  }
  return out;
}

vector<vector<double>> mZ(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B) {
  //Berechnet die Matrix Z, mit Zij = Z(i,j,...)
  vector<vector<double>> out;
  int imax = thetaN.size();
  out.resize(imax);
  for (int i = 0; i < imax; i++) {
    out[i].resize(imax);
  }

  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        out[i][j] = -out[j][i];
      } else {
        if (i == j) {
          out[i][j] = 0;
        } else {
          out[i][j] = Z(i, j, thetaN, thetaNp, rho, A, B);
        }
      }
    }
  }
  return out;
}

vector<vector<double>> mY(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B) {
  //Berechnet die Matrix Y, mit Yij = Y(i,j,...)
  vector<vector<double>> out;
  int imax = thetaN.size();
  out.resize(imax);
  for (int i = 0; i < imax; i++) {
    out[i].resize(imax);
  }

  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        out[i][j] = out[j][i];
      } else {
        out[i][j] = Y(i, j, thetaN, thetaNp, rho, A, B);
      }
    }
  }
  return out;
}

double V(int r, int j, vector<double> &thetaN, vector<double> &thetaNp,
    int iter, double h, double mu, double t, double rho) {
  //Berechnet den Wert V aus Gleichung 33
  double out = 0;
  int imax = thetaN.size();
  double ki;
  double xpki;
  double ypki;
  vector<double> temp;
  for (int i = 0; i < imax; i++) {
    ki = i * h;
    xpki = 0;
    ypki = 0;
    for (int k = 0; k < imax; k++) {
      xpki = xpki - thetaNp[k] * SinInt(k, ki, h, thetaN, iter);
      ypki = ypki + thetaNp[k] * CosInt(k, ki, h, thetaN, iter);
    }

    temp = dversch(ki, t);
    xpki = xpki + temp[0];
    ypki = ypki + temp[1];
    if ((xpki * xpki + ypki * ypki) != 0)
      out = out
          - Reibf(ki, h, rho) * mu
              * (SinInt(j, ki, h, thetaN, iter)
                  * SinInt(r, ki, h, thetaN, iter)
                  + CosInt(j, ki, h, thetaN, iter)
                      * CosInt(r, ki, h, thetaN, iter))
              / sqrt(xpki * xpki + ypki * ypki);
  }

  return out;
}

vector<vector<double>> mV(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double l, double mu, double t, int iter, double rho) {
  vector<vector<double>> out;
  //Berechnet die Matrix V, mit Vij = V(i,j,...)
  int imax = thetaN.size();
  out.resize(imax);
  for (int i = 0; i < imax; i++) {
    out[i].resize(imax);
  }

  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        out[i][j] = out[j][i];
      } else {
        out[i][j] = V(i, j, thetaN, thetaNp, iter, h, mu, t, rho);
      }
    }
  }
  return out;
}

double u(int r, vector<double> &thetaN, vector<double> &thetaNp, int iter,
    double h, double mu, double t, double rho) {
  //Berechnet den Wert u aus Gleichung 34
  double out = 0;
  int imax = thetaN.size();
  double ki;
  double xpki;
  double ypki;
  vector<double> temp;
  for (int i = 0; i < imax; ++i) {
    ki = i * h;
    xpki = 0;
    ypki = 0;
    for (int k = 0; k < imax; ++k) {
      xpki = xpki - thetaNp[k] * SinInt(k, ki, h, thetaN, iter);
      ypki = ypki + thetaNp[k] * CosInt(k, ki, h, thetaN, iter);
    }

    temp = dversch(ki, t);
    xpki = xpki + temp[0];
    ypki = ypki + temp[1];
    if ((xpki * xpki + ypki * ypki) != 0)
      out = out
          - Reibf(ki, h, rho) * mu
              * (-temp[0] * SinInt(r, ki, h, thetaN, iter)
                  + temp[1] * CosInt(r, ki, h, thetaN, iter))
              / sqrt(xpki * xpki + ypki * ypki);

  }

  return out;
}

vector<double> uVek(vector<double> &thetaN, vector<double> &thetaNp, double h,
    double l, double mu, double t, int iter, double rho) {
  //Berechnet den Vektor u mit ui= u(i,...)
  vector<double> out;

  int imax = thetaN.size();
  out.resize(imax);

  for (int i = 0; i < imax; i++) {

    out[i] = u(i, thetaN, thetaNp, iter, h, mu, t, rho);

  }
  return out;
}

double verschVekIntegrand(double x, int i, double h, vector<double> &thetaN, int iter, double l, double t)
{
  return -SinInt(i,x,h,thetaN,iter)*ddversch(x,t)[0]+CosInt(i,x,h,thetaN,iter)*ddversch(x,t)[1];
}

vector<double> verschVek(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double l, double t, double rho, int iter) {
  //Berechnet den Vektor v aus Gleichung 25
  vector<double> out;
  int imax = thetaN.size();
  out.resize(imax);
  for (int i = 0; i < imax; ++i) {
    out[i] =
        Integral(verschVekIntegrand, 0, l, iter,
                 i,h,thetaN,iter,l,t);
  }
  return out;
}

vector<vector<double>> convertMMatrix(vector<vector<double>> A) {
  //Verändert die Massenmatrix so, dass sich der erste Winkel nur aus der gegebenen Funktion berechnet
  size_t n = A.size();
  vector<vector<double>> B;
  B.resize(n);
  for (size_t i = 0; i < n; ++i) {
    B[i].resize(n);
    for (size_t j = 0; j < n; ++j) {
      if (i == 0) {
        if (j == 0) {
          B[i][j] = 1;
        } else {
          B[i][j] = 0;
        }
      } else {
        B[i][j] = A[i][j];
      }

    }
  }
  return B;
}
