#include "utility/utility.h"

#include <iostream>
#include <fstream>

#include "terms/problem_definition.h"

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

void save(vector<double> &x, double t, string path) {
  //Speichert den Vektor x in einer csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei \"" << path << "\" konnte nicht geöffnet werden" << endl;
    return;
  }
  vector<double> verschiebung = versch(0,t);
  file << t << "," << verschiebung[0] << "," << verschiebung[1] << "," << theta0(t);
  for (size_t i = 0; i < x.size(); ++i) {
    file << "," << x[i];
  }
  file << "\n";

}

void save(double &x, string path) {
  //Speichert den Wert x in einer csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei \"" << path << "\" konnte nicht geöffnet werden" << endl;
    return;
  }

  file << x << "\n";
  file.close();
}

void save(int n, double Rflex, double L, double rho, double mu,
          int iter, bool ver, bool reib, bool winkelkontrolle, string path)
{
  //Speichert die Werte in der csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei \"" << path << "\" konnte nicht geöffnet werden" << endl;
    return;
  }

  file << n << "," << Rflex << "," << L << "," << rho << "," << mu << ","
    << iter << "," << ver << "," << reib << "," << winkelkontrolle << ","
    << t00 << "," << t01 << "," << t02 << "," << vx1 << "," << vy1 << ","
    << vx2 << "," << vy2 << "\n";

  file.close();
}

void finish(time_t t1, string path) {
  //Speichert die Zeitdifferent aus t1 und der aktuellen Zeit in einer csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei \"" << path << "\" konnte nicht geöffnet werden" << endl;
    return;
  }
  time_t t2 = time(0) - t1;
  file << t2;
  file.close();

  cout << "Dauer: " << t2 << " s" << endl;
}
