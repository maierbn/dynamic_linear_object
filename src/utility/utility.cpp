#include "utility/utility.h"

#include <iostream>
#include <fstream>

void plotteMatrix(vector<vector<double>> &Matrix) {
//Plottet eine MxN Matrix
  int imax = Matrix.size();
  int jmax = Matrix[0].size();
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < jmax; j++) {
      cout << Matrix[i][j] << " ";
    }
    cout << endl;
  }

}

void save(vector<double> &x, string path) {
  //Speichert den Vektor x in einer csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei konnte nicht geöffnet werden" << endl;
    return;
  }
  for (size_t i = 0; i < x.size(); ++i) {
    if (i != 0) {
      file << ",";
    }
    file << x[i];
  }
  file << "\n";

}

void save(double &x, string path) {
  //Speichert den Wert x in einer csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei konnte nicht geöffnet werden" << endl;
    return;
  }

  file << x << "\n";

}

void finish(time_t t1, string path) {
  //Speichert die Zeitdifferent aus t1 und der aktuellen Zeit in einer csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei konnte nicht geöffnet werden" << endl;
    return;
  }
  time_t t2 = time(0) - t1;
  file << t2;
  cout << t2 << "s" << endl;

}
