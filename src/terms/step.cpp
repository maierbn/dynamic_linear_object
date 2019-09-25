#include "terms/step.h"

#include "terms/terms.h"
#include "terms/problem_definition.h"
#include "utility/operators.h"
#include "numerics/linear_solve.h"
#include "terms/globals.h"

vector<double> step(int n, double Rflex, double L, double rho, double mu,
    bool ver, bool reib, bool winkelkontrolle, vector<double> &thetaN,
    vector<double> &thetaNp, int iter, double t, string path) {
  //berechnet einen Schritt einer Iteration zur Lösung der ODE
  vector<double> temp(n + 1);
  vector<double> erg(n + 1);
  double h = L / n;
  double hInv = 1.0 / h;
  if (winkelkontrolle) {
    thetaN[0] = theta0(t);
    thetaNp[0] = theta0p(t);
  }
  Km = K(n, Rflex, h, hInv);
  Am = A(L, h, hInv, n, thetaN, iter);
  Bm = B(L, h, hInv, n, thetaN, iter);
  Zm = mZ(thetaN, thetaNp, rho, Am, Bm);
  Ym = mY(thetaN, thetaNp, rho, Am, Bm);

  temp = temp - (multipl(Km, thetaN)) + (multipl(Zm, thetaNp))
      + (multipl(Ym, thetaNp));
  if (reib) {
    Vm = mV(thetaN, thetaNp, h, hInv, L, mu, t, iter, rho);

    temp = temp + multipl(Vm, thetaNp);
    if (ver) {
      uv = uVek(thetaN, thetaNp, h, hInv, L, mu, t, iter, rho);

      temp = temp + uv;
    }
  }
  if (ver) {
    verschv = verschVek(thetaN, thetaNp, h, hInv, L, t, rho, iter);
    temp = temp - verschv;
  }
  vector<vector<double>> Mm = mM(rho, h, hInv, thetaN, L, iter);
  if (winkelkontrolle) {
    Mm = convertMMatrix(Mm);
    temp[0] = theta0pp(t);
  }
  erg = solve(Mm, temp);
  return erg;

}
