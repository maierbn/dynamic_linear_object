#include "terms/problem_definition.h"

//Für zusätzliche Bewegung und Kontrolle über den Winkel
double vx2;
double vx1;
double vy2;
double vy1;
double t02;
double t01;
double t00;

//Start des Programms

vector<double> versch(double s, double t) {
  //Gibt die aktuelle aufgebrachte Verschiebung zum Zeitpunkt t an
  vector<double> out;
  out.resize(2);
  out[0] = vx2 * t * t + vx1 * t;
  out[1] = vy2 * t * t + vy1 * t;
  return out;
}
vector<double> dversch(double s, double t) {
  //Gibt die aktuelle 1. Ableitung der aufgebrachten Verschiebung zum Zeitpunkt t an

  vector<double> out;
  out.resize(2);
  out[0] = vx2 * 2 * t + vx1;
  out[1] = vy2 * 2 * t + vy1;
  return out;
}
vector<double> ddversch(double s, double t) {
  //Gibt die aktuelle 2. Ableitung der aufgebrachten Verschiebung zum Zeitpunkt t an
  vector<double> out;
  out.resize(2);
  out[0] = vx2 * 2;
  out[1] = vy2 * 2;
  return out;
}
double Reibf(double s, double h, double rho) {
  //Gibt die Reibung an der Stelle s in [0,L] an //TODO Reibf
  /*
  if(s<0.5){
    return 0.0;
  }
  return 30-(s-0.5)*40;*/
  return 10.0;
  //return 9.81 * rho * h;
}
double theta0(double t) {
  //Gibt den aktuellen aufgebrachten Winkel im Ursprung zum Zeitpunkt t an
  return t02 * t * t + t01 * t + t00;
}
double theta0p(double t) {
  //Gibt die 1. Ableitung des aktuellen aufgebrachten Winkels im Ursprung zum Zeitpunkt t an
  return t02 * 2 * t + t01;
}
double theta0pp(double t) {
  //Gibt die 2. Ableitung des aktuellen aufgebrachten Winkels im Ursprung zum Zeitpunkt t an
  return t02 * 2;
}
