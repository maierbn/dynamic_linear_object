#pragma once

#include <vector>
#include <iostream>
#include <chrono>

using namespace std;

// save vector x and time point t in a csv file at path
void save(vector<double> &x, double t, string path);

// save vector x in a csv file at path
void save(double &x, string path);

// save parameter values in a csv file at path
void save(int n, double Rflex, double L, double rho, double mu,
          bool enablePrescribedDisplacement, bool enableFriction, bool enablePrescribedAngle, string path);

// save duration computed from current time and t1 in csv file at path
void finalize(std::chrono::time_point<std::chrono::steady_clock> t1, string path);

//Load look-up table containing prescribed displacements and angles for the gripper position
//File format is csv based with separator ",", each line contains: t,x,y,phi,x',y',phi',x'',y'',phi'' (i.e. also 1st and 2nd derivatives)
void loadLUT(std::string path);
