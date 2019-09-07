#pragma once

#include <vector>

using namespace std;

//Hütchenfunktion bzw lineare Basisfunktion -> 1 für i*h, linear abfallend auf 0 zu (i-1)*h, (i+1)*h, ansonsten konstant 0
double N(int i, double h, double s);

//Theta gibt anhand Lineakombination der N+1 Koeffizienten mit den jeweiligen Basisfunktionen
//den Wert von Theta an der Stelle s zurück
double Theta(const vector<double> &thetaN, double h, double s);

//Berechnet das Integral über den Cosinus von Theta multipliziert mit der i-ten Basisfunktion
double CosInt(int i, double s, double h, const vector<double> &thetaN,
    int iter);

//Berechnet das Integral über den Sinus von Theta multipliziert mit der i-ten Basisfunktion
double SinInt(int i, double s, double h, const vector<double> &thetaN,
    int iter);

//Berechnet den Wert m_(i,k) aus Gleichung 14
double m(int i, int k, double rho, double h, const vector<double> &thetaN,
    double l, int iter);

//Berechnet die Ableitung der Funktion SInInt  nach thetar
double dSinInt(int i, int r, double s, double h, const vector<double> &thetaN,
    int iter);

//Berechnet die Ableitung der Funktion CosInt nach thetar
double dCosInt(int i, int r, double s, double h, const vector<double> &thetaN,
    int iter);

//Berechnet für alle i,k,r das Integral über dSinInt(i,r)*SinInt(k) zur Wiederverwendung
vector<vector<vector<double>>> A(double l, double h, int n,
    vector<double> &thetaN, int iter);

//Berechnet für alle i,k,r das Integral über dCosInt(i,r)*CosInt(k) zur Wiederverwendung
vector<vector<vector<double>>> B(double l, double h, int n,
    vector<double> &thetaN, int iter);

//Berechnet die Steifigkeitsmatrix des Kabels mit N Stützstellen mit Abstand h und der Steifigkeit RFlex
vector<vector<double>> K(int n, double RFlex, double h);

//Berechnet die Ableitung der Funktion m(i,k) nach Theta_r
double dm(int i, int k, int r, double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

//Berechnet den Wer Z aus Gleichung 16
double Z(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

//Berechnet den Wer Y aus Gleichung 16
double Y(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

//Berechnet die Matrix M, mit Mij = m(i,j,...)
vector<vector<double>> mM(double rho, double h, vector<double> &thetaN,
    double l, int iter);

//Berechnet die Matrix Z, mit Zij = Z(i,j,...)
vector<vector<double>> mZ(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

//Berechnet die Matrix Y, mit Yij = Y(i,j,...)
vector<vector<double>> mY(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

//Berechnet den Wert V aus Gleichung 33
double V(int r, int j, vector<double> &thetaN, vector<double> &thetaNp,
    int iter, double h, double mu, double t, double rho);

//Berechnet die Matrix V, mit Vij = V(i,j,...)
vector<vector<double>> mV(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double l, double mu, double t, int iter, double rho);

//Berechnet den Wert u aus Gleichung 34
double u(int r, vector<double> &thetaN, vector<double> &thetaNp, int iter,
    double h, double mu, double t, double rho);

//Berechnet den Vektor u mit ui= u(i,...)
vector<double> uVek(vector<double> &thetaN, vector<double> &thetaNp, double h,
    double l, double mu, double t, int iter, double rho);

//Berechnet den Vektor v aus Gleichung 25
vector<double> verschVek(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double l, double t, double rho, int iter);

//Verändert die Massenmatrix so, dass sich der erste Winkel nur aus der gegebenen Funktion berechnet
vector<vector<double>> convertMMatrix(vector<vector<double>> A);
