#pragma once

#include <vector>

using namespace std;

// Addition zweier Vektoren
vector<double> operator+(vector<double> a, vector<double> b);

// Substraktion zweier Vektoren
vector<double> operator-(vector<double> a, vector<double> b);

// Matrix-Vektor-Multiplikation
vector<double> multipl(vector<vector<double>> &A, vector<double> &x);

// Skalierung eines Vektors
vector<double> multipl(double a, vector<double> &b);
