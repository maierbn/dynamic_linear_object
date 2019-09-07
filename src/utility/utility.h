#pragma once

#include <vector>
#include <iostream>

using namespace std;

//Plottet eine MxN Matrix
void plotteMatrix(vector<vector<double>> &Matrix);

//Speichert den Vektor x in einer csv Datei am Ort path
void save(vector<double> &x, string path);

//Speichert den Wert x in einer csv Datei am Ort path
void save(double &x, string path);

//Speichert die Zeitdifferent aus t1 und der aktuellen Zeit in einer csv Datei am Ort path
void finish(time_t t1, string path);
