#include "numerics/timestepping.h"

#include "utility/utility.h"
#include "utility/operators.h"
#include "terms/step.h"
#include "terms/globals.h"
#include "terms/problem_definition.h"

#include <iostream>

// forward euler to solve ODE
void euler(int n, double Rflex, double L, double rho, double mu, bool enablePrescribedDisplacement,
    bool enableFriction, bool enablePrescribedAngle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, string path,
    double dt, double tend)
{
  const double t0 = 0;   // start time

  int nIterations = (tend - t0) / dt;   // number of  iterations
  dt = (tend - t0) / nIterations;
  double t = t0;
  vector<double> thetaNtemp = thetaNstart;
  vector<double> thetaNptemp = thetaNpstart;
  vector<double> omega(n + 1);
  save(thetaNtemp, t, path);

  cout << 100 * (t - t0) / (tend - t0) << "%" << endl;

  for (int var = 0; var < nIterations; ++var)
  {
    omega = step(n, Rflex, L, rho, mu, enablePrescribedDisplacement, enableFriction, enablePrescribedAngle,
        thetaNtemp, thetaNptemp, t, path);
    t = t + dt;

    thetaNtemp = thetaNtemp + multipl(dt, thetaNptemp);
    thetaNptemp = thetaNptemp + multipl(dt, omega);

    if (enablePrescribedAngle)
    {
      thetaNtemp[0] = theta0(t);
      thetaNptemp[0] = theta0p(t);
    }
    save(thetaNtemp, t, path);

    cout << 100 * (t - t0) / (tend - t0) << "%" << endl;
  }

  double h = L / n;
  save(h, path);
}

void rungeKutta4(int n, double Rflex, double L, double rho, double mu, bool enablePrescribedDisplacement,
    bool enableFriction, bool enablePrescribedAngle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, string path,
    double dt, double tend)
{
  // Runge-Kutta-4 method to solve ODE
  const double t0 = 0;   // start time

  int nIterations = (tend - t0) / dt;   // number of iterations
  dt = (tend - t0) / nIterations;
  int saveFreq = max(1.0,0.01/dt);    // write to file every 0.01s
  double t = t0;

  vector<double> thetaNtemp = thetaNstart;
  vector<double> thetaNptemp = thetaNpstart;
  vector<double> temp(n + 1);
  vector<double> tempp(n + 1);

  for (int var = 0; var < nIterations; ++var)
  {
    if (var % 100 == 0)
    {
      cout << 100 * (t - t0) / (tend - t0) << "%" << endl;
    }
    if (var % saveFreq == 0)
      save(thetaNtemp, t, path);

    k1 = thetaNptemp;
    k1p = step(n, Rflex, L, rho, mu, enablePrescribedDisplacement, enableFriction, enablePrescribedAngle, thetaNtemp,
        thetaNptemp, t, path);

    temp = thetaNtemp + multipl(dt / 2, k1);
    tempp = thetaNptemp + multipl(dt / 2, k1p);
    k2 = tempp;
    k2p = step(n, Rflex, L, rho, mu, enablePrescribedDisplacement, enableFriction, enablePrescribedAngle, temp,
        tempp, t + dt / 2, path);

    temp = thetaNtemp + multipl(dt / 2, k2);
    tempp = thetaNptemp + multipl(dt / 2, k2p);
    k3 = tempp;
    k3p = step(n, Rflex, L, rho, mu, enablePrescribedDisplacement, enableFriction, enablePrescribedAngle, temp,
        tempp, t + dt / 2, path);

    temp = thetaNtemp + multipl(dt, k3);
    tempp = thetaNptemp + multipl(dt, k3p);
    k4 = tempp;
    k4p = step(n, Rflex, L, rho, mu, enablePrescribedDisplacement, enableFriction, enablePrescribedAngle, temp,
        tempp, t + dt, path);

    t = t + dt;
    thetaNtemp = thetaNtemp + multipl(dt / 6, k1) + multipl(dt / 3, k2)
        + multipl(dt / 3, k3) + multipl(dt / 6, k4);
    thetaNptemp = thetaNptemp + multipl(dt / 6, k1p) + multipl(dt / 3, k2p)
        + multipl(dt / 3, k3p) + multipl(dt / 6, k4p);

    if (enablePrescribedAngle)
    {
      thetaNtemp[0] = theta0(t);
      thetaNptemp[0] = theta0p(t);
    }
  }

  save(thetaNtemp, t, path);
  cout << 100 * (t - t0) / (tend - t0) << "%" << endl;

  double h = L / n;
  save(h, path);
}
