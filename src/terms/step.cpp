#include "terms/step.h"

#include "terms/terms.h"
#include "terms/problem_definition.h"
#include "utility/operators.h"
#include "numerics/linear_solve.h"
#include "terms/globals.h"

// compute a single timestep of the timestepping iteration to solve the ODE
vector<double> step(int n, double Rflex, double L, double rho, double mu,
    bool enablePrescribedDisplacement, bool enableFriction, bool enablePrescribedAngle, vector<double> &thetaN,
    vector<double> &thetaNp, double t, string path)
{
  vector<double> temp(n + 1);
  vector<double> result;

  double h = L / n;
  double hInv = 1.0 / h;

  // set the prescribed value of the angle
  if (enablePrescribedAngle)
  {
    thetaN[0] = theta0(t);
    thetaNp[0] = theta0p(t);
  }

  // compute all coefficient matrices
  computeMatrixK(n, Rflex, h, hInv);
  computeMatrixA(L, h, hInv, n, thetaN);
  computeMatrixB(L, h, hInv, n, thetaN);
  computeMatrixZ(thetaN, thetaNp, rho, Am, Bm);
  computeMatrixY(thetaN, thetaNp, rho, Am, Bm);

  // compute Z ω_n + Y ω_n - K Θ_n, for right hand side of Eq. (14.2)
  temp = (multipl(Zm, thetaNp)) + (multipl(Ym, thetaNp)) -(multipl(Km, thetaN));

  // account for friction force
  if (enableFriction)
  {

    // add + V ω_n
    computeMatrixV(thetaN, thetaNp, h, hInv, L, mu, t, rho);
    temp = temp + multipl(Vm, thetaNp);

    if (enablePrescribedDisplacement)
    {
      // add + u
      uV = uVector(thetaN, thetaNp, h, hInv, L, mu, t, rho);
      temp = temp + uV;
    }
  }

  if (enablePrescribedDisplacement)
  {
    // add - w
    wV = wVector(thetaN, thetaNp, h, hInv, L, t, rho);
    temp = temp - wV;
  }

  // compute M matrix for left hand side of Eq. (14.2)
  computeMatrixM(rho, h, hInv, thetaN, L);

  // account for prescribed angle
  if (enablePrescribedAngle)
  {
    convertMMatrix();
    temp[0] = theta0pp(t);
  }

  // solve the linear system with matrix M, Eq. (14.2)
  result = solve(Mm, temp);

  return result;

}
