#include "numerics/linear_solve.h"

#include "alglib/solvers.h"
#include "alglib/ap.h"
#include "alglib/linalg.h"

#include "terms/globals.h"

using namespace alglib;

vector<double> solve(vector<vector<double>> &A, vector<double> &b) {
  //berechnet den Lösungsvektor x des LGS Ax=b
  int n = b.size();
  vector<double> x(n);
  Adata.resize(n * n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Adata[i * n + j] = A[i][j];
    }
  }
  real_2d_array h;
  real_1d_array xl;
  real_1d_array rhs;
  h.attach_to_ptr(n, n, Adata.data());
  rhs.attach_to_ptr(n, b.data());
  vector<double> xlsg(n);

  xl.attach_to_ptr(n, xlsg.data());
  ae_int_t info = 0;
  densesolverreport rep;

  rmatrixsolve(h, n, rhs, info, rep, xl);
  for (int i = 0; i < n; ++i) {

    x[i] = xl[i];
  }
  return x;
}
