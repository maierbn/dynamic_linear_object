#pragma once

#include <vector>
#include <iostream>

using namespace std;

// forward euler to solve ODE
void euler(int n, double Rflex, double L, double rho, double mu, bool enablePrescribedDisplacement,
    bool enableFriction, bool enablePrescribedAngle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, string path,
    double dt, double tend);

// Runge-Kutta-4 method to solve ODE
void rungeKutta4(int n, double Rflex, double L, double rho, double mu, bool enablePrescribedDisplacement,
    bool enableFriction, bool enablePrescribedAngle, vector<double> thetaNstart,
    const vector<double> thetaNpstart, string path,
    double dt, double tend);
