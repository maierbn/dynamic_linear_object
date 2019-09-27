#include "terms/terms.h"

#include "terms/globals.h"
#include "utility/operators.h"
#include "numerics/integral.h"
#include "terms/problem_definition.h"

#include <cmath>
#include <iostream>
#include <fstream>

double N(int i, double h, double hInv, double s) {
  //Hütchenfunktion bzw lineare Basisfunktion -> 1 für i*h, linear abfallend auf 0 zu (i-1)*h, (i+1)*h, ansonsten konstant 0

  double shInv = s * hInv;
  int is = shInv;
  is -= i;

  double erg = 0;     // falls s < (i-1)*h   or   s > (i+1)*h
  switch (is)
  {
  case -1:            // falls (i-1)*h < s < i*h:     erg = (s - (i - 1) * h) / h
    //erg = (s - (i - 1) * h) * hInv;
    erg = shInv - i + 1;
    break;

  case 0:             // falls i*h <= s < (i+1)*h:    erg = ((i + 1) * h - s) / h
    //erg = ((i + 1) * h - s) * hInv;
    erg = i + 1 - shInv;
    break;

  default:
    break;
  }
  return erg;
}
double Theta(const vector<double> &thetaN, double h, double hInv, double s) {
  //Theta gibt anhand Linearkombination der N+1 Koeffizienten mit den jeweiligen Basisfunktionen
  //den Wert von Theta an der Stelle s zurück
  double shInv = s * hInv;
  int i1 = shInv;
  int i2 = i1+1;

  double alpha = (shInv - i1);
  //return alpha * thetaN[i2] + (1-alpha) * thetaN[i1];
  return alpha * (thetaN[i2]-thetaN[i1]) + thetaN[i1];

#if 0
  double alpha = (s*hInv - i1);
  double erg1 = alpha * thetaN[i2] + (1-alpha) * thetaN[i1];
  // Annahme, dass s immer im richtigen intervall liegt und thetaN groß genug ist.
  double erg2 = thetaN[i1] * N(i1, h, hInv, s) + thetaN[i2] * N(i2, h, hInv, s);

  if (fabs(erg1 - erg2) > 1e-16)
    cout << alpha << ": " << erg1 << "," << erg2 << endl;

  return erg2;
  // Code ohne diese Annahme:
  double erg = 0;
  if (i1 >= 0 && i1 <= thetaN.size() - 1) {
    double erg2 = N(i1, h, hInv, s);
    erg = erg + thetaN[i1] * erg2;
  }
  if (i2 >= 0 && i2 <= thetaN.size() - 1 && i1 != i2) {
    double erg3 = N(i2, h, hInv, s);
    erg = erg + thetaN[i2] * erg3;
  }

  return erg;
#endif
}

double CosIntIntegrand(double x, int i, double h, double hInv, const vector<double> &thetaN) {
  //Berechnet den Integranden für CosInt
  return cos(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x);
}

double CosInt(int i, double s, double h, double hInv, const vector<double> &thetaN,
    int iter) {
  //Berechnet das Integral über den Cosinus von Theta multipliziert mit der i-ten Basisfunktion
  double x1 = max(i - (double)1.0, (double)0.0) * h;
  double x2 = min((i + (double)1.0) * h, s);

  return Integral(CosIntIntegrand, x1, x2, iter,
                  CosIntIntegrand(x1, i, h, hInv, thetaN),
                  i, h, hInv, thetaN);

  /*
   return Integral([i,h,thetaN](double x) {
   return cos(Theta(thetaN,h,x))*N(i,h,hInv,x);
   }, 0, s, iter);*/
}

double SinIntIntegrand(double x, int i, double h, double hInv, const vector<double> &thetaN) {
  //Berechnet den Integranden für SinInt
  return sin(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x);
}

double SinInt(int i, double s, double h, double hInv, const vector<double> &thetaN,
    int iter) {
  //Berechnet das Integral über den Sinus von Theta multipliziert mit der i-ten Basisfunktion

  double x1 = max(i - (double)1.0, (double)0.0) * h;
  double x2 = min((i + (double)1.0) * h, s);

  return Integral(SinIntIntegrand, x1, x2, iter,
                  SinIntIntegrand(x1, i, h, hInv, thetaN),
                  i, h, hInv, thetaN);
  /*
   return Integral([i,h,thetaN](double x) {
   return sin(Theta(thetaN,h,x))*N(i,h,hInv,x);
   }, 0, s, iter);*/
}

double mIntegrand(double x, int i, int k, double h, double hInv, int iter, const vector<double> &thetaN)
{
  return SinInt( i,x,h,hInv,thetaN,iter)*SinInt( k,x,h,hInv,thetaN,iter)+CosInt( i,x,h,hInv,thetaN,iter)*CosInt( k,x,h,hInv,thetaN,iter);
}

double m(int i, int k, double rho, double h, double hInv, const vector<double> &thetaN,
    double l, int iter) {

  //Berechnet den Wert m_(i,k) aus Gleichung 14
  return rho
      * Integral(mIntegrand, 0, l, iter, 0,
                 i,k,h,hInv,iter,thetaN);
}

double dSinIntIntegrand(double x, int i, int r, double h, double hInv, const vector<double> &thetaN)
{
  return cos(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x)*N(r,h,hInv,x);
}

double dSinInt(int i, int r, double s, double h, double hInv, const vector<double> &thetaN,
    int iter) {
  //Berechnet die Ableitung der Funktion SInInt  nach thetar

  if (abs(r - i) <= 1 && s > (r - 1) * h && s > (i - 1) * h) {
    double x1 = max(min(r, i) - (double)1.0, (double)0.0) * h;
    double x2 = min((max(r, i) + (double)1.0) * h, s);

    return Integral(dSinIntIntegrand, x1, x2, iter,
                    dSinIntIntegrand(x1, i,r,h,hInv,thetaN),
                    i,r,h,hInv,thetaN);
  } else {
    return 0.0;
  }
  /*
   return Integral([i,r,h,thetaN](double x) {
   return cos(Theta(thetaN,h,x))*N(i,h,hInv,x)*N(r,h,x);
   }, 0, s, iter);*/
}

double dCosIntIntegrand(double x, int i, int r, double h, double hInv, const vector<double> &thetaN)
{
  return -sin(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x)*N(r,h,hInv,x);
}

double dCosInt(int i, int r, double s, double h, double hInv, const vector<double> &thetaN,
    int iter) {
  //Berechnet die Ableitung der Funktion CosInt nach thetar

  if (abs(r - i) <= 1 && s > (r - 1) * h && s > (i - 1) * h) {
    double x1 = max(min(r, i) - (double)1.0, (double)0.0) * h;
    double x2 = min((max(r, i) + (double)1.0) * h, s);

    return Integral(dCosIntIntegrand, x1, x2, iter,
                    dCosIntIntegrand(x1, i,r,h,hInv,thetaN),
                    i,r,h,hInv,thetaN);

  } else {
    return 0.0;
  }/*
   return Integral([i,r,h,thetaN](double x) {
   return -sin(Theta(thetaN,h,x))*N(i,h,hInv,x)*N(r,h,x);
   }, 0, s, iter);*/
}

double AIntegrand(double x, int i, int r, int k, double h, double hInv, vector<double> &thetaN, int iter)
{
  return dSinInt(i,r,x,h,hInv,thetaN,iter)*SinInt(k,x,h,hInv,thetaN,iter);
}

void computeMatrixA(double l, double h, double hInv, int n,
    vector<double> &thetaN, int iter) {
  //Berechnet für alle i,k,r das Integral über dSinInt(i,r)*SinInt(k) zur Wiederverwendung

  if (Am.empty())
  {
    Am.resize(n + 1);
    for (int i = 0; i <= n; i++) {
      Am[i].resize(n + 1);
      for (int k = 0; k <= n; k++) {
        Am[i][k].resize(n + 1);
      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(dynamic) shared(Am)
  for (int i = 0; i <= n; i++) {
    for (int k = 0; k <= n; k++) {
      for (int r = 0; r <= n; r++) {
        if (r >= i) {

          Am[i][k][r] =
              Integral(AIntegrand, 0, l, iter,
                       0,
                       i,r,k,h,hInv,thetaN,iter);

          Am[r][k][i] = Am[i][k][r];
        }
      }
    }
  }
/*
  #pragma omp parallel for collapse(3) schedule(dynamic) shared(Am)
  for (int i = 0; i <= n; i++) {
    for (int k = 0; k <= n; k++) {
      for (int r = 0; r <= n; r++) {
        if (r < i) {
          Am[i][k][r] = Am[r][k][i];
        }
      }
    }
  }*/
}

double BIntegrand(double x, int i, int r, int k, double h, double hInv, vector<double> &thetaN, int iter)
{
  return dCosInt(i,r,x,h,hInv,thetaN,iter)*CosInt(k,x,h,hInv,thetaN,iter);
}

void computeMatrixB(double l, double h, double hInv, int n,
    vector<double> &thetaN, int iter) {
  //Berechnet für alle i,k,r das Integral über dCosInt(i,r)*CosInt(k) zur Wiederverwendung

  if (Bm.empty())
  {
    Bm.resize(n + 1);
    for (int i = 0; i <= n; i++) {
      Bm[i].resize(n + 1);
      for (int k = 0; k <= n; k++) {
        Bm[i][k].resize(n + 1);
      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(dynamic) shared(Bm)
  for (int i = 0; i <= n; i++) {
    for (int k = 0; k <= n; k++) {
      for (int r = 0; r <= n; r++) {
        if (r >= i) {

          Bm[i][k][r] =
              Integral(BIntegrand, 0, l, iter,
                       0,
                       i,r,k,h,hInv,thetaN,iter);
          Bm[r][k][i] = Bm[i][k][r];
        }
      }
    }
  }
/*
  #pragma omp parallel for collapse(3) schedule(dynamic) shared(Bm)
  for (int i = 0; i <= n; i++) {
    for (int k = 0; k <= n; k++) {
      for (int r = 0; r <= n; r++) {
        if (r < i) {
          Bm[i][k][r] = Bm[r][k][i];
        }
      }
    }
  }*/
}

void computeMatrixK(int n, double RFlex, double h, double hInv) {
  //Berechnet die Steifigkeitsmatrix des Kabels mit N Stützstellen mit Abstand h und der Steifigkeit RFlex

  if (Km.empty())
  {
    Km.resize(n + 1);
    for (int i = 0; i <= n; i++) {
      Km[i].resize(n + 1);
    }
  }

  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      if (i == j && i > 0 && i < n) {
        Km[i][j] = 2 * RFlex * hInv;
      } else if (i == j) {
        Km[i][j] = RFlex * hInv;

      } else if (j - i == 1 || i - j == 1) {
        Km[i][j] = -RFlex * hInv;

      }
    }
  }
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
void computeMatrixM(double rho, double h, double hInv, vector<double> &thetaN,
    double l, int iter) {
  //Berechnet die Matrix M, mit Mij = m(i,j,...)

  int imax = thetaN.size();
  if (Mm.empty())
  {
    Mm.resize(imax);
    for (int i = 0; i < imax; i++) {
      Mm[i].resize(imax);
    }
  }

  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Mm)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j >= i) {
        Mm[i][j] = m(i, j, rho, h, hInv, thetaN, l, iter);
        Mm[j][i] = Mm[i][j];
      }
    }
  }
  /*
  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Mm)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        Mm[i][j] = Mm[j][i];
      }
    }
  }*/
}

void computeMatrixZ(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B) {
  //Berechnet die Matrix Z, mit Zij = Z(i,j,...)

  int imax = thetaN.size();
  if (Zm.empty())
  {
    Zm.resize(imax);
    for (int i = 0; i < imax; i++) {
      Zm[i].resize(imax);
    }
  }

  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Zm)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j >= i) {
        if (i == j) {
          Zm[i][j] = 0;
        } else {
          Zm[i][j] = Z(i, j, thetaN, thetaNp, rho, A, B);
          Zm[j][i] = -Zm[i][j];
        }
      }
    }
  }
/*
  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Zm)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        Zm[i][j] = -Zm[j][i];
      }
    }
  }*/
}

void computeMatrixY(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B) {
  //Berechnet die Matrix Y, mit Yij = Y(i,j,...)

  int imax = thetaN.size();
  if (Ym.empty())
  {
    Ym.resize(imax);
    for (int i = 0; i < imax; i++) {
      Ym[i].resize(imax);
    }
  }

  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Ym)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j >= i) {
        Ym[i][j] = Y(i, j, thetaN, thetaNp, rho, A, B);
        Ym[j][i] = Ym[i][j];
      }
    }
  }

  /*
  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Ym)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        Ym[i][j] = Ym[j][i];
      }
    }
  }
  */
}

double V(int r, int j, vector<double> &thetaN, vector<double> &thetaNp,
    int iter, double h, double hInv, double mu, double t, double rho) {
  //Berechnet den Wert V aus Gleichung 33
  double out = 0;
  int imax = thetaN.size();
  double ki;
  double xpki;
  double ypki;
  array<double,2> prescribed_displacement;
  for (int i = 0; i < imax; i++) {
    ki = i * h;
    xpki = 0;
    ypki = 0;
    for (int k = 0; k < imax; k++) {
      xpki = xpki - thetaNp[k] * SinInt(k, ki, h, hInv, thetaN, iter);
      ypki = ypki + thetaNp[k] * CosInt(k, ki, h, hInv, thetaN, iter);
    }

    prescribed_displacement = dversch(ki, t);
    xpki = xpki + prescribed_displacement[0];
    ypki = ypki + prescribed_displacement[1];
    if ((xpki * xpki + ypki * ypki) != 0)
      out = out
          - Reibf(ki, h, rho) * mu
              * (SinInt(j, ki, h, hInv, thetaN, iter)
                  * SinInt(r, ki, h, hInv, thetaN, iter)
                  + CosInt(j, ki, h, hInv, thetaN, iter)
                      * CosInt(r, ki, h, hInv, thetaN, iter))
              / sqrt(xpki * xpki + ypki * ypki);
  }

  return out;
}

void computeMatrixV(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double l, double mu, double t, int iter, double rho) {
  //Berechnet die Matrix V, mit Vij = V(i,j,...)

  int imax = thetaN.size();
  if (Vm.empty())
  {
    Vm.resize(imax);
    for (int i = 0; i < imax; i++) {
      Vm[i].resize(imax);
    }
  }

  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Vm)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j >= i) {
        Vm[i][j] = V(i, j, thetaN, thetaNp, iter, h, hInv, mu, t, rho);
        Vm[j][i] = Vm[i][j];
      }
    }
  }
/*
  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Vm)
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < imax; j++) {
      if (j < i) {
        Vm[i][j] = Vm[j][i];
      }
    }
  }*/
}

double u(int r, vector<double> &thetaN, vector<double> &thetaNp, int iter,
    double h, double hInv, double mu, double t, double rho) {
  //Berechnet den Wert u aus Gleichung 34
  double out = 0;
  int imax = thetaN.size();
  double ki;
  double xpki;
  double ypki;
  array<double,2> d_displacement;
  for (int i = 0; i < imax; ++i) {
    ki = i * h;
    xpki = 0;
    ypki = 0;
    for (int k = 0; k < imax; ++k) {
      xpki = xpki - thetaNp[k] * SinInt(k, ki, h, hInv, thetaN, iter);
      ypki = ypki + thetaNp[k] * CosInt(k, ki, h, hInv, thetaN, iter);
    }

    d_displacement = dversch(ki, t);
    xpki = xpki + d_displacement[0];
    ypki = ypki + d_displacement[1];
    if ((xpki * xpki + ypki * ypki) != 0)
      out = out
          - Reibf(ki, h, rho) * mu
              * (-d_displacement[0] * SinInt(r, ki, h, hInv, thetaN, iter)
                  + d_displacement[1] * CosInt(r, ki, h, hInv, thetaN, iter))
              / sqrt(xpki * xpki + ypki * ypki);

  }

  return out;
}

vector<double> uVek(vector<double> &thetaN, vector<double> &thetaNp, double h, double hInv,
    double l, double mu, double t, int iter, double rho) {
  //Berechnet den Vektor u mit ui= u(i,...)
  vector<double> out;

  int imax = thetaN.size();
  out.resize(imax);

  for (int i = 0; i < imax; i++) {

    out[i] = u(i, thetaN, thetaNp, iter, h, hInv, mu, t, rho);

  }
  return out;
}

double verschVekIntegrand(double x, int i, double h, double hInv, vector<double> &thetaN, int iter, double l, double t)
{
  return -SinInt(i,x,h,hInv, thetaN,iter)*ddversch(x,t)[0]+CosInt(i,x,h,hInv, thetaN,iter)*ddversch(x,t)[1];
}

vector<double> verschVek(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double l, double t, double rho, int iter) {
  //Berechnet den Vektor v aus Gleichung 25
  vector<double> out;
  int imax = thetaN.size();
  out.resize(imax);
  for (int i = 0; i < imax; ++i) {

    out[i] =
        Integral(verschVekIntegrand, 0, l, iter, 0,
                 i,h,hInv,thetaN,iter,l,t);
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

void convertMMatrix() {
  //Verändert die Massenmatrix so, dass sich der erste Winkel nur aus der gegebenen Funktion berechnet
  size_t n = Mm.size();

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i == 0) {
        if (j == 0) {
          Mm[i][j] = 1;
        } else {
          Mm[i][j] = 0;
        }
      }
    }
  }
}
