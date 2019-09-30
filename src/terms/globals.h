#pragma once

#include <vector>

using namespace std;

// declaration of global variables that are used during computation
extern double tend;                       //< end time of simulation
extern vector<vector<vector<double>>> Am; //< A matrix, precomputed values to use for Y and Z
extern vector<vector<vector<double>>> Bm; //< B matrix, precomputed values to use for Y and Z
extern vector<vector<double>> Zm;         //< Z matrix
extern vector<vector<double>> Km;         //< K matrix
extern vector<vector<double>> Vm;         //< V matrix
extern vector<vector<double>> Ym;         //< Y matrix
extern vector<vector<double>> Mm;         //< M matrix
extern vector<double> uV;                 //< u vector
extern vector<double> wV;                 //< w vector
extern vector<double> Adata;              //< buffer used for linear system solve

// variables needed for the Runge-Kutta-4 scheme
extern vector<double> k1;
extern vector<double> k2;
extern vector<double> k3;
extern vector<double> k4;
extern vector<double> k1p;
extern vector<double> k2p;
extern vector<double> k3p;
extern vector<double> k4p;
