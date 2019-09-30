#pragma once

#include <vector>

using namespace std;

// vector addiotn
vector<double> operator+(vector<double> a, vector<double> b);

// vector substraction
vector<double> operator-(vector<double> a, vector<double> b);

// matrix vector multiplication
vector<double> multipl(vector<vector<double>> &A, vector<double> &x);

// scaling of a vector
vector<double> multipl(double a, vector<double> &b);
