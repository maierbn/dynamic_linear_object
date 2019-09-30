#pragma once

#include <vector>

using namespace std;

// solve linear system Ax = b for x
// this uses the alglib numeric library
vector<double> solve(vector<vector<double>> &A, vector<double> &b);
