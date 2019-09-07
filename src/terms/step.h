#pragma once

#include <vector>
#include <iostream>

using namespace std;

//berechnet einen Schritt einer Iteration zur Lösung der ODE
vector<double> step(int n, double Rflex, double L, double rho, double mu,
    bool ver, bool reib, bool winkelkontrolle, vector<double> &thetaN,
    vector<double> &thetaNp, int iter, double t, string path);
