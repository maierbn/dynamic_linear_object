#pragma once

#include <vector>
#include <iostream>
#include <chrono>

using namespace std;

//Plottet eine MxN Matrix
void plotteMatrix(vector<vector<double>> &Matrix);

//Speichert den Vektor x und den Zeitpunkt t in einer csv Datei am Ort path
void save(vector<double> &x, double t, string path);

//Speichert den Wert x in einer csv Datei am Ort path
void save(double &x, string path);

//Speichert die Werte in der csv Datei am Ort path
void save(int n, double Rflex, double L, double rho, double mu,
          int iter, bool ver, bool reib, bool winkelkontrolle, string path);

//Speichert die Zeitdifferent aus t1 und der aktuellen Zeit in einer csv Datei am Ort path
void finish(std::chrono::time_point<std::chrono::steady_clock> t1, string path);
