#pragma once

#include <vector>

using namespace std;

// Definition aller globalen Variablen, die intern für die Berechnung
// als Zwischenspeicher verwendet werden

extern vector<vector<vector<double>>> Am;
extern vector<vector<vector<double>>> Bm;
extern vector<vector<double>> Zm;
extern vector<vector<double>> Km;
extern vector<vector<double>> Vm;
extern vector<vector<double>> Ym;
extern vector<vector<double>> Mm;
extern vector<double> uv;
extern vector<double> verschv;
extern vector<double> Adata;

//Für Runge Kutta 4
extern vector<double> k1;
extern vector<double> k2;
extern vector<double> k3;
extern vector<double> k4;
extern vector<double> k1p;
extern vector<double> k2p;
extern vector<double> k3p;
extern vector<double> k4p;
