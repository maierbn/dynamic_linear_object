#include "numerics/timestepping.h"

#include "utility/utility.h"
#include "utility/operators.h"
#include "terms/step.h"
#include "terms/globals.h"
#include "terms/problem_definition.h"

#include <iostream>

void euler(int n, double Rflex, double L, double rho, double mu, bool ver,
    bool reib, bool winkelkontrolle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, int iter, double t0, string path,
    string path2, double dt, double tend) {
  //Explizites Euler Verfahren zur Lösung der ODE

  int wdh = (tend - t0) / dt;   // Anzahl Iterationen
  dt = (tend - t0) / wdh;
  double t = t0;
  vector<double> thetaNtemp = thetaNstart;
  vector<double> thetaNptemp = thetaNpstart;
  vector<double> omega(n + 1);
  save(thetaNtemp, t, path);
  cout << 100 * (t - t0) / (tend - t0) << "%" << endl;
  for (int var = 0; var < wdh; ++var) {
    omega = step(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle,
        thetaNtemp, thetaNptemp, iter, t, path);
    t = t + dt;

    thetaNtemp = thetaNtemp + multipl(dt, thetaNptemp);
    thetaNptemp = thetaNptemp + multipl(dt, omega);
    if (winkelkontrolle) {
      thetaNtemp[0] = theta0(t);
      thetaNptemp[0] = theta0p(t);
    }
    save(thetaNtemp, t, path);
    cout << 100 * (t - t0) / (tend - t0) << "%" << endl;

  }
  double h = L / n;
  save(h, path);

}

void rungeKutta4(int n, double Rflex, double L, double rho, double mu, bool ver,
    bool reib, bool winkelkontrolle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, int iter, double t0, string path,
    string path2, double dt, double tend) {
  //Runge-Kutta-4 Verfahren zur Lösung der ODE

  int wdh = (tend - t0) / dt;   // Anzahl Iterationen
  dt = (tend - t0) / wdh;
  int saveFreq = max(1.0,0.01/dt);    // Schreibe Datei alle 0.01 Simulations-Sekunden
  double t = t0;
  vector<double> thetaNtemp = thetaNstart;
  vector<double> thetaNptemp = thetaNpstart;
  vector<double> temp(n + 1);
  vector<double> tempp(n + 1);
  for (int var = 0; var < wdh; ++var) {

    if (var % 100 == 0)
    {
      cout << 100 * (t - t0) / (tend - t0) << "%" << endl;
    }
    if (var % saveFreq == 0)
      save(thetaNtemp, t, path);

    k1 = thetaNptemp;
    k1p = step(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaNtemp,
        thetaNptemp, iter, t, path);

    temp = thetaNtemp + multipl(dt / 2, k1);
    tempp = thetaNptemp + multipl(dt / 2, k1p);
    k2 = tempp;
    k2p = step(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, temp,
        tempp, iter, t + dt / 2, path);

    temp = thetaNtemp + multipl(dt / 2, k2);
    tempp = thetaNptemp + multipl(dt / 2, k2p);
    k3 = tempp;
    k3p = step(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, temp,
        tempp, iter, t + dt / 2, path);

    temp = thetaNtemp + multipl(dt, k3);
    tempp = thetaNptemp + multipl(dt, k3p);
    k4 = tempp;
    k4p = step(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, temp,
        tempp, iter, t + dt, path);

    t = t + dt;
    thetaNtemp = thetaNtemp + multipl(dt / 6, k1) + multipl(dt / 3, k2)
        + multipl(dt / 3, k3) + multipl(dt / 6, k4);
    thetaNptemp = thetaNptemp + multipl(dt / 6, k1p) + multipl(dt / 3, k2p)
        + multipl(dt / 3, k3p) + multipl(dt / 6, k4p);
    if (winkelkontrolle) {
      thetaNtemp[0] = theta0(t);
      thetaNptemp[0] = theta0p(t);
    }
  }
  save(thetaNtemp, t, path);
  cout << 100 * (t - t0) / (tend - t0) << "%" << endl;

  double h = L / n;
  save(h, path);

}
