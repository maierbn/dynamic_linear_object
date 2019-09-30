#include "terms/terms.h"

#include "terms/globals.h"
#include "utility/operators.h"
#include "numerics/integral.h"
#include "terms/problem_definition.h"

#include <cmath>
#include <iostream>
#include <fstream>

// linear hat function, Eq. (3)
double N(int i, double h, double hInv, double s)
{
  // linear nodal basis function -> 1 for i*h, linear decreasing to 0 at (i-1)*h, (i+1)*h, else constant 0

  double shInv = s * hInv;
  int is = shInv;
  is -= i;

  double result = 0;      // if s < (i-1)*h   or   s > (i+1)*h
  switch (is)
  {
  case -1:                // if (i-1)*h < s < i*h:     result = (s - (i - 1) * h) / h
    //result = (s - (i - 1) * h) * hInv;
    result = shInv - i + 1;
    break;

  case 0:                 // if i*h <= s < (i+1)*h:    result = ((i + 1) * h - s) / h
    //result = ((i + 1) * h - s) * hInv;
    result = i + 1 - shInv;
    break;

  default:
    break;
  }
  return result;
}

// discretized angle
double Theta(const vector<double> &thetaN, double h, double hInv, double s)
{
  // discretized angle theta
  double shInv = s * hInv;
  int i1 = shInv;
  int i2 = i1+1;

  double alpha = (shInv - i1);

  return alpha * (thetaN[i2]-thetaN[i1]) + thetaN[i1];
}

double CosIntIntegrand(double x, int i, double h, double hInv, const vector<double> &thetaN) {
  //Berechnet den Integranden für CosInt
  return cos(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x);
}

// compute integral over cosine of theta, multiplied with ith basis function, Eq. (7.1)
double CosInt(int i, double s, double h, double hInv, const vector<double> &thetaN)
{
  double x1 = max(i - (double)1.0, (double)0.0) * h;
  double x2 = min((i + (double)1.0) * h, s);

  return Integral(CosIntIntegrand, x1, x2,
                  CosIntIntegrand(x1, i, h, hInv, thetaN),
                  i, h, hInv, thetaN);
}

double SinIntIntegrand(double x, int i, double h, double hInv, const vector<double> &thetaN) {
  //Berechnet den Integranden für SinInt
  return sin(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x);
}

// compute integral over sine of theta, multiplied with ith basis function, Eq. (7.2)
double SinInt(int i, double s, double h, double hInv, const vector<double> &thetaN)
{
  double x1 = max(i - (double)1.0, (double)0.0) * h;
  double x2 = min((i + (double)1.0) * h, s);

  return Integral(SinIntIntegrand, x1, x2,
                  SinIntIntegrand(x1, i, h, hInv, thetaN),
                  i, h, hInv, thetaN);
}

double mIntegrand(double x, int i, int k, double h, double hInv, const vector<double> &thetaN)
{
  return SinInt( i,x,h,hInv,thetaN)*SinInt( k,x,h,hInv,thetaN)+CosInt( i,x,h,hInv,thetaN)*CosInt( k,x,h,hInv,thetaN);
}

// entry of mass matrix, m_(i,k), Eq. (6.1)
double m(int i, int k, double rho, double h, double hInv, const vector<double> &thetaN, double l)
{
  return rho
      * Integral(mIntegrand, 0, l,  0,
                 i,k,h,hInv,thetaN);
}

double dSinIntIntegrand(double x, int i, int r, double h, double hInv, const vector<double> &thetaN)
{
  return cos(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x)*N(r,h,hInv,x);
}

// dSinInt/dTheta
double dSinInt(int i, int r, double s, double h, double hInv, const vector<double> &thetaN)
{
  if (abs(r - i) <= 1 && s > (r - 1) * h && s > (i - 1) * h)
  {
    double x1 = max(min(r, i) - (double)1.0, (double)0.0) * h;
    double x2 = min((max(r, i) + (double)1.0) * h, s);

    return Integral(dSinIntIntegrand, x1, x2,
                    dSinIntIntegrand(x1, i,r,h,hInv,thetaN),
                    i,r,h,hInv,thetaN);
  }
  else
  {
    return 0.0;
  }
}

double dCosIntIntegrand(double x, int i, int r, double h, double hInv, const vector<double> &thetaN)
{
  return -sin(Theta(thetaN,h,hInv,x))*N(i,h,hInv,x)*N(r,h,hInv,x);
}

// dCosInt/dTheta
double dCosInt(int i, int r, double s, double h, double hInv, const vector<double> &thetaN)
{
  if (abs(r - i) <= 1 && s > (r - 1) * h && s > (i - 1) * h)
  {
    double x1 = max(min(r, i) - (double)1.0, (double)0.0) * h;
    double x2 = min((max(r, i) + (double)1.0) * h, s);

    return Integral(dCosIntIntegrand, x1, x2,
                    dCosIntIntegrand(x1, i,r,h,hInv,thetaN),
                    i,r,h,hInv,thetaN);
  }
  else
  {
    return 0.0;
  }
}

double AIntegrand(double x, int i, int r, int k, double h, double hInv, vector<double> &thetaN)
{
  return dSinInt(i,r,x,h,hInv,thetaN)*SinInt(k,x,h,hInv,thetaN);
}

// compute the integral over dSinInt(i,r)*SinInt(k) for reuse, for all i,k,r
void computeMatrixA(double l, double h, double hInv, int n,
    vector<double> &thetaN)
{
  // initialize memory for A matrix, if not yet done
  if (Am.empty())
  {
    Am.resize(n + 1);
    for (int i = 0; i <= n; i++)
    {
      Am[i].resize(n + 1);
      for (int k = 0; k <= n; k++)
      {
        Am[i][k].resize(n + 1);
      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(dynamic) shared(Am)
  for (int i = 0; i <= n; i++)
  {
    for (int k = 0; k <= n; k++)
    {
      for (int r = 0; r <= n; r++)
      {
        if (r >= i)
        {
          Am[i][k][r] =
              Integral(AIntegrand, 0, l,
                       0,
                       i,r,k,h,hInv,thetaN);

          Am[r][k][i] = Am[i][k][r];
        }
      }
    }
  }
}

double BIntegrand(double x, int i, int r, int k, double h, double hInv, vector<double> &thetaN)
{
  return dCosInt(i,r,x,h,hInv,thetaN)*CosInt(k,x,h,hInv,thetaN);
}

// compute the integral over dCosInt(i,r)*CosInt(k) for reuse, for all i,k,r
void computeMatrixB(double l, double h, double hInv, int n,
    vector<double> &thetaN)
{
  // initialize memory for B matrix, if not yet done
  if (Bm.empty())
  {
    Bm.resize(n + 1);
    for (int i = 0; i <= n; i++)
    {
      Bm[i].resize(n + 1);
      for (int k = 0; k <= n; k++)
      {
        Bm[i][k].resize(n + 1);
      }
    }
  }

  #pragma omp parallel for collapse(3) schedule(dynamic) shared(Bm)
  for (int i = 0; i <= n; i++)
  {
    for (int k = 0; k <= n; k++)
    {
      for (int r = 0; r <= n; r++)
      {
        if (r >= i)
        {
          Bm[i][k][r] =
              Integral(BIntegrand, 0, l,
                       0,
                       i,r,k,h,hInv,thetaN);
          Bm[r][k][i] = Bm[i][k][r];
        }
      }
    }
  }
}

// compute the stiffness matrix K
void computeMatrixK(int n, double RFlex, double h, double hInv)
{
  // initialize memory for K matrix, if not yet done
  if (Km.empty())
  {
    Km.resize(n + 1);
    for (int i = 0; i <= n; i++)
    {
      Km[i].resize(n + 1);
    }
  }

  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j <= n; j++)
    {
      if (i == j && i > 0 && i < n)
      {
        Km[i][j] = 2 * RFlex * hInv;
      }
      else if (i == j)
      {
        Km[i][j] = RFlex * hInv;
      }
      else if (j - i == 1 || i - j == 1)
      {
        Km[i][j] = -RFlex * hInv;
      }
    }
  }
}

// compute the derivative of m(i,k) w.r.t Theta_r
double dm(int i, int k, int r, double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B)
{
  return rho * (A[i][k][r] + A[k][i][r] + B[i][k][r] + B[k][i][r]);
}

// compute Z_r,k, Eq. (6.2)
double Z(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B)
{
  int imax = thetaN.size();
  double temp = 0;
  for (int i = 0; i < imax; i++)
  {
    temp += 0.5 * (dm(i, k, r, rho, A, B) - dm(r, i, k, rho, A, B))
        * thetaNp[i];
  }

  return temp;
}

// compute Y_r,k, Eq. (6.3)
double Y(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B)
{
  int imax = thetaN.size();
  double temp = 0;
  for (int i = 0; i < imax; i++)
  {
    temp -= 0.5 * dm(r, i, k, rho, A, B) * thetaNp[i];
  }

  return temp;
}
// compute mass matrix M with M_ij = m(i,j,...)
void computeMatrixM(double rho, double h, double hInv, vector<double> &thetaN, double l)
{
  //Berechnet die Matrix M, mit Mij = m(i,j,...)

  int imax = thetaN.size();
  // initialize memory for M matrix, if not yet done
  if (Mm.empty())
  {
    Mm.resize(imax);
    for (int i = 0; i < imax; i++)
    {
      Mm[i].resize(imax);
    }
  }

  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Mm)
  for (int i = 0; i < imax; i++)
  {
    for (int j = 0; j < imax; j++)
    {
      if (j >= i)
      {
        Mm[i][j] = m(i, j, rho, h, hInv, thetaN, l);
        Mm[j][i] = Mm[i][j];
      }
    }
  }
}

// compute matrix Z, with Z_ij = Z(i,j,...)
void computeMatrixZ(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B)
{

  int imax = thetaN.size();
  // initialize memory for Z matrix, if not yet done
  if (Zm.empty())
  {
    Zm.resize(imax);
    for (int i = 0; i < imax; i++)
    {
      Zm[i].resize(imax);
    }
  }

  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Zm)
  for (int i = 0; i < imax; i++)
  {
    for (int j = 0; j < imax; j++)
    {
      if (j >= i)
      {
        if (i == j)
        {
          Zm[i][j] = 0;
        }
        else
        {
          Zm[i][j] = Z(i, j, thetaN, thetaNp, rho, A, B);
          Zm[j][i] = -Zm[i][j];
        }
      }
    }
  }
}

// compute matrix Y, with Y_ij = Y(i,j,...)
void computeMatrixY(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B)
{
  int imax = thetaN.size();
  // initialize memory for Y matrix, if not yet done
  if (Ym.empty())
  {
    Ym.resize(imax);
    for (int i = 0; i < imax; i++)
    {
      Ym[i].resize(imax);
    }
  }

  #pragma omp parallel for collapse(2) schedule(dynamic) shared(Ym)
  for (int i = 0; i < imax; i++)
  {
    for (int j = 0; j < imax; j++)
    {
      if (j >= i)
      {
        Ym[i][j] = Y(i, j, thetaN, thetaNp, rho, A, B);
        Ym[j][i] = Ym[i][j];
      }
    }
  }
}

// compute entry V_rj, Eq. (11)
double V(int r, int j, vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double mu, double t, double rho)
{
  double out = 0;
  int imax = thetaN.size();
  double ki;
  double xpki;
  double ypki;
  array<double,2> prescribed_displacement;

  for (int i = 0; i < imax; i++)
  {
    ki = i * h;
    xpki = 0;
    ypki = 0;

    for (int k = 0; k < imax; k++)
    {
      xpki = xpki - thetaNp[k] * SinInt(k, ki, h, hInv, thetaN);
      ypki = ypki + thetaNp[k] * CosInt(k, ki, h, hInv, thetaN);
    }

    prescribed_displacement = ddispl(ki, t);
    xpki = xpki + prescribed_displacement[0];
    ypki = ypki + prescribed_displacement[1];

    if ((xpki * xpki + ypki * ypki) != 0)
    {
      out = out
          - frictionForce(ki, h, rho) * mu
              * (SinInt(j, ki, h, hInv, thetaN)
                  * SinInt(r, ki, h, hInv, thetaN)
                  + CosInt(j, ki, h, hInv, thetaN)
                      * CosInt(r, ki, h, hInv, thetaN))
              / sqrt(xpki * xpki + ypki * ypki);
    }
  }

  return out;
}

// compute matrix V, with V_ij = V(i,j,...)
void computeMatrixV(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double l, double mu, double t, double rho)
{

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
        Vm[i][j] = V(i, j, thetaN, thetaNp,  h, hInv, mu, t, rho);
        Vm[j][i] = Vm[i][j];
      }
    }
  }
}

// compute entry u_r, Eq. (12)
double u(int r, vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double mu, double t, double rho)
{
  double out = 0;
  int imax = thetaN.size();
  double ki;
  double xpki;
  double ypki;
  array<double,2> d_displacement;

  for (int i = 0; i < imax; ++i)
  {
    ki = i * h;
    xpki = 0;
    ypki = 0;

    for (int k = 0; k < imax; ++k)
    {
      xpki = xpki - thetaNp[k] * SinInt(k, ki, h, hInv, thetaN);
      ypki = ypki + thetaNp[k] * CosInt(k, ki, h, hInv, thetaN);
    }

    d_displacement = ddispl(ki, t);
    xpki = xpki + d_displacement[0];
    ypki = ypki + d_displacement[1];

    if ((xpki * xpki + ypki * ypki) != 0)
    {
      out = out
          - frictionForce(ki, h, rho) * mu
              * (-d_displacement[0] * SinInt(r, ki, h, hInv, thetaN)
                  + d_displacement[1] * CosInt(r, ki, h, hInv, thetaN))
              / sqrt(xpki * xpki + ypki * ypki);
    }
  }

  return out;
}

// compute the vector u with u_i= u(i,...)
vector<double> uVector(vector<double> &thetaN, vector<double> &thetaNp, double h, double hInv,
    double l, double mu, double t, double rho)
{
  vector<double> out;
  int imax = thetaN.size();
  out.resize(imax);

  for (int i = 0; i < imax; i++)
  {
    out[i] = u(i, thetaN, thetaNp,  h, hInv, mu, t, rho);
  }
  return out;
}

// integrand for vector w, Eq. (9)
double wVectorIntegrand(double x, int i, double h, double hInv, vector<double> &thetaN, double l, double t)
{
  return -SinInt(i,x,h,hInv, thetaN)*dddispl(x,t)[0]+CosInt(i,x,h,hInv, thetaN)*dddispl(x,t)[1];
}

// compute vector w, Eq. (9)
vector<double> wVector(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double l, double t, double rho)
{
  vector<double> out;
  int imax = thetaN.size();
  out.resize(imax);
  for (int i = 0; i < imax; ++i)
  {
    out[i] =
        Integral(wVectorIntegrand, 0, l,  0,
                 i,h,hInv,thetaN,l,t);
  }
  return out;
}

// adjust mass matrix M, such that the first angle, \theta_0 is given by the prescribed angle function in problem_definition.cpp
void convertMMatrix()
{
  size_t n = Mm.size();

  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j)
    {
      if (i == 0)
      {
        if (j == 0)
        {
          Mm[i][j] = 1;
        }
        else
        {
          Mm[i][j] = 0;
        }
      }
    }
  }
}
