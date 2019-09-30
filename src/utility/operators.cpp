#include "utility/operators.h"

vector<double> operator+(vector<double> a, vector<double> b) {
  vector<double> c(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    c[i] = a[i] + b[i];
  }
  return c;

}
vector<double> operator-(vector<double> a, vector<double> b) {
  vector<double> c(a.size());
  for (size_t i = 0; i < a.size(); ++i)
  {
    c[i] = a[i] - b[i];
  }
  return c;

}
vector<double> multipl(vector<vector<double>> &A, vector<double> &x) {
  vector<double> c(A.size());
  for (size_t i = 0; i < A.size(); i++)
  {
    for (size_t j = 0; j < x.size(); j++)
    {
      c[i] = c[i] + A[i][j] * x[j];
    }
  }
  return c;

}
vector<double> multipl(double a, vector<double> &b) {
  vector<double> c(b.size());
  for (size_t i = 0; i < b.size(); ++i)
  {
    c[i] = a * b[i];
  }
  return c;
}
