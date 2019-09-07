#pragma once

#include <vector>

using namespace std;

//berechnet den Lösungsvektor x des LGS Ax=b
vector<double> solve(vector<vector<double>> &A, vector<double> &b);
