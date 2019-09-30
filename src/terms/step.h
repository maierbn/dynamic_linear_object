#pragma once

#include <vector>
#include <iostream>

using namespace std;

// compute a single timestep of the timestepping iteration to solve the ODE
vector<double> step(int n, double Rflex, double L, double rho, double mu,
    bool enablePrescribedDisplacement, bool enableFriction, bool enablePrescribedAngle, vector<double> &thetaN,
    vector<double> &thetaNp, double t, string path);
