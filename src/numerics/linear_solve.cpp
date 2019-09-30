#include "numerics/linear_solve.h"

#include "alglib/solvers.h"
#include "alglib/ap.h"
#include "alglib/linalg.h"

#include "terms/globals.h"

using namespace alglib;

vector<double> solve(vector<vector<double>> &A, vector<double> &b)
{
  // solve linear system Ax = b for x

  int n = b.size();
  vector<double> x(n);

  // populate A
  Adata.resize(n * n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      Adata[i * n + j] = A[i][j];
    }
  }

  // define helper vectors
  real_2d_array h;
  real_1d_array xl;
  real_1d_array rhs;
  h.attach_to_ptr(n, n, Adata.data());
  rhs.attach_to_ptr(n, b.data());
  vector<double> xlsg(n);

  xl.attach_to_ptr(n, xlsg.data());
  ae_int_t info = 0;
  densesolverreport rep;

  // solve system
  rmatrixsolve(h, n, rhs, info, rep, xl);

  // assign result
  for (int i = 0; i < n; ++i)
  {
    x[i] = xl[i];
  }
  return x;
}
