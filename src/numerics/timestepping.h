#pragma once

#include <vector>
#include <iostream>

using namespace std;

//Explizites Euler Verfahren zur Lösung der ODE
void euler(int n, double Rflex, double L, double rho, double mu, bool ver,
    bool reib, bool winkelkontrolle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, int iter, double t0, string path,
    string path2, double dt, double tend);

//Runge-Kutta-4 Verfahren zur Lösung der ODE
void rungeKutta4(int n, double Rflex, double L, double rho, double mu, bool ver,
    bool reib, bool winkelkontrolle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, int iter, double t0, string path,
    string path2, double dt, double tend);
