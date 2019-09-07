#pragma once

#include <vector>

using namespace std;

// Die folgenden Randbedingungen werden gesetzt:

// aufgeprägte Verschiebung:
// verschiebung(t) = v2 * t^2 + v1 * t,    mit v2=(vx2,vy2), v1=(vx1,vy1)

// aufgeprägter Winkel:
// theta0(t) = t02 * t^2 + t01 * t + t00

// Reibkoeffizient:
// Reibf(s, h, rho)     Gibt die Reibung an der Stelle s in [0,L] an

//Für zusätzliche Bewegung und Kontrolle über den Winkel
extern double vx2;
extern double vx1;
extern double vy2;
extern double vy1;
extern double t02;
extern double t01;
extern double t00;

//Gibt die aktuelle aufgebrachte Verschiebung zum Zeitpunkt t an
vector<double> versch(double s, double t);

//Gibt die aktuelle 1. Ableitung der aufgebrachten Verschiebung zum Zeitpunkt t an
vector<double> dversch(double s, double t);

//Gibt die aktuelle 2. Ableitung der aufgebrachten Verschiebung zum Zeitpunkt t an
vector<double> ddversch(double s, double t);

//Gibt die Reibung an der Stelle s in [0,L] an
double Reibf(double s, double h, double rho);

//Gibt den aktuellen aufgebrachten Winkel im Ursprung zum Zeitpunkt t an
double theta0(double t);

//Gibt die 1. Ableitung des aktuellen aufgebrachten Winkels im Ursprung zum Zeitpunkt t an
double theta0p(double t);

//Gibt die 2. Ableitung des aktuellen aufgebrachten Winkels im Ursprung zum Zeitpunkt t an
double theta0pp(double t);
