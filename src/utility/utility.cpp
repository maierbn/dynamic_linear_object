#include "utility/utility.h"

#include <iostream>
#include <fstream>
#include <tuple>

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
  array<double,2> verschiebung = versch(0,t);
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

void finish(std::chrono::time_point<std::chrono::steady_clock> t1, string path) {
  //Speichert die Zeitdifferent aus t1 und der aktuellen Zeit in einer csv Datei am Ort path
  ofstream file;
  file.open(path.c_str(), ios::out | ios::app);
  if (!file.is_open()) {
    cout << "Datei \"" << path << "\" konnte nicht geöffnet werden" << endl;
    return;
  }
  auto now = std::chrono::steady_clock::now();
  auto t2 = std::chrono::duration_cast<chrono::milliseconds>(now - t1).count();
  file << t2/1000.0;
  file.close();

  cout << "Dauer: " << t2 << " ms" << endl;
}

void loadLUT(std::string path)
{
  //File format is csv based with separator ",", each line contains: t,x,y,theta,x',y',theta',x'',y'',theta'' (i.e. also 1st and 2nd derivatives)
  ifstream file(path);
  if (!file.is_open())
  {
    cout << "Could not load file \"" << path << "\".";
  }
  for (int lineNo = 0; !file.eof(); lineNo++)
  {
    std::string line;
    getline(file,line);

    // load time point
    double t = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    // load values
    double x = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    double y = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    double angle = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    // load 1st derivatives
    double dx = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    double dy = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    double dangle = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    // load 2nd derivatives
    double ddx = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    double ddy = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    double ddangle = atof(line.c_str());
    line.erase(0,line.find(",")+1);
    line.erase(0,line.find_first_not_of(" \t")+1);

    // store to look-up table
    if (lineNo == 0)
    {
      LUTtimestepWidth = t;
    }
    /*
    extern std::vector<std::tuple<double,double,double>> displacementsX;
    extern std::vector<std::tuple<double,double,double>> displacementsY;
    extern std::vector<std::tuple<double,double,double>> thetas;*/
    displacementsX.push_back(std::make_tuple(x,dx,ddx));
    displacementsY.push_back(std::make_tuple(y,dy,ddy));
    thetas.push_back(std::make_tuple(angle,dangle,ddangle));
  }

  file.close();
}
