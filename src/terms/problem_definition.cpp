#include "terms/problem_definition.h"

#include <tuple>
#include <cmath>

#include "terms/globals.h"

//Für zusätzliche Bewegung und Kontrolle über den Winkel
double vx2;
double vx1;
double vy2;
double vy1;
double t02;
double t01;
double t00;
bool useLUTdisplacement;
bool useLUTtheta;
double LUTtimestepWidth;
bool useSmoothFunction;
std::vector<std::tuple<double,double,double>> displacementsX;
std::vector<std::tuple<double,double,double>> displacementsY;
std::vector<std::tuple<double,double,double>> thetas;

using std::pow;

//Start des Programms

array<double,2> versch(double s, double t) {
  //Gibt die aktuelle aufgebrachte Verschiebung zum Zeitpunkt t an
  array<double,2> out;

  if (useLUTdisplacement)
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // if there is no more data, repeat the last recorded value
    if (index+1 >= displacementsX.size())
    {
      index = displacementsX.size()-2;
      alpha = 1;
    }

    out[0] = (1-alpha) * std::get<0>(displacementsX[index]) + alpha * std::get<0>(displacementsX[index+1]);
    out[1] = (1-alpha) * std::get<0>(displacementsY[index]) + alpha * std::get<0>(displacementsY[index+1]);
  }
  else if (useSmoothFunction)
  {
    double tt = t/tend;
    double s = pow(tt,4)*35. + pow(tt,5)*-84. + pow(tt,6)*70. + pow(tt,7)*-20.;

    out[0] = s*vx1*tend;
    out[1] = s*vy1*tend;
  }
  else
  {
    out[0] = vx2 * t * t + vx1 * t;
    out[1] = vy2 * t * t + vy1 * t;
  }
  return out;
}

array<double,2> dversch(double s, double t) {
  //Gibt die aktuelle 1. Ableitung der aufgebrachten Verschiebung zum Zeitpunkt t an
  array<double,2> out;

  if (useLUTdisplacement)
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // if there is no more data, repeat the last recorded value
    if (index+1 >= displacementsX.size())
    {
      index = displacementsX.size()-2;
      alpha = 1;
    }

    out[0] = (1-alpha) * std::get<1>(displacementsX[index]) + alpha * std::get<1>(displacementsX[index+1]);
    out[1] = (1-alpha) * std::get<1>(displacementsY[index]) + alpha * std::get<1>(displacementsY[index+1]);
  }
  else if (useSmoothFunction)
  {
    double tt = t/tend;
    double ds = 1./tend * (pow(tt,3)*4*35. + pow(tt,4)*5*-84. + pow(tt,5)*6*70. + pow(tt,6)*7*-20.);

    out[0] = ds*vx1*tend;
    out[1] = ds*vy1*tend;
  }
  else
  {
    out[0] = vx2 * 2 * t + vx1;
    out[1] = vy2 * 2 * t + vy1;
  }
  return out;
}

array<double,2> ddversch(double s, double t) {
  //Gibt die aktuelle 2. Ableitung der aufgebrachten Verschiebung zum Zeitpunkt t an
  array<double,2> out;

  if (useLUTdisplacement)
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // if there is no more data, repeat the last recorded value
    if (index+1 >= displacementsX.size())
    {
      index = displacementsX.size()-2;
      alpha = 1;
    }

    out[0] = (1-alpha) * std::get<2>(displacementsX[index]) + alpha * std::get<2>(displacementsX[index+1]);
    out[1] = (1-alpha) * std::get<2>(displacementsY[index]) + alpha * std::get<2>(displacementsY[index+1]);
  }
  else if (useSmoothFunction)
  {
    double tt = t/tend;
    double dds = 1./tend*1./tend * (pow(tt,2)*3*4*35. + pow(tt,3)*4*5*-84. + pow(tt,4)*5*6*70. + pow(tt,5)*6*7*-20.);

    out[0] = dds*vx1*tend;
    out[1] = dds*vy1*tend;
  }
  else
  {
    out[0] = vx2 * 2;
    out[1] = vy2 * 2;
  }
  return out;
}
double Reibf(double s, double h, double rho) {
  //Gibt die Reibung an der Stelle s in [0,L] an
  /*
  if(s<0.5){
    return 0.0;
  }
  return 30-(s-0.5)*40;*/
  return 9.81 * rho * h;
}
double theta0(double t) {
  //Gibt den aktuellen aufgebrachten Winkel im Ursprung zum Zeitpunkt t an

  if (useLUTtheta)
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // if there is no more data, repeat the last recorded value
    if (index+1 >= thetas.size())
    {
      index = thetas.size()-2;
      alpha = 1;
    }

    return (1-alpha) * std::get<0>(thetas[index]) + alpha * std::get<0>(thetas[index+1]);
  }
  else if (useSmoothFunction)
  {
    double tt = t/tend;
    double s = pow(tt,4)*35. + pow(tt,5)*-84. + pow(tt,6)*70. + pow(tt,7)*-20.;

    return s*t01*tend + t00;
  }
  else
  {
    return t02 * t * t + t01 * t + t00;
  }
}

double theta0p(double t) {
  //Gibt die 1. Ableitung des aktuellen aufgebrachten Winkels im Ursprung zum Zeitpunkt t an

  if (useLUTtheta)
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // if there is no more data, repeat the last recorded value
    if (index+1 >= thetas.size())
    {
      index = thetas.size()-2;
      alpha = 1;
    }

    return (1-alpha) * std::get<1>(thetas[index]) + alpha * std::get<1>(thetas[index+1]);
  }
  else if (useSmoothFunction)
  {
    double tt = t/tend;
    double ds = 1./tend * (pow(tt,3)*4*35. + pow(tt,4)*5*-84. + pow(tt,5)*6*70. + pow(tt,6)*7*-20.);

    return ds*t01*tend;
  }
  else
  {
    return t02 * 2 * t + t01;
  }
}

double theta0pp(double t) {
  //Gibt die 2. Ableitung des aktuellen aufgebrachten Winkels im Ursprung zum Zeitpunkt t an

  if (useLUTtheta)
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // if there is no more data, repeat the last recorded value
    if (index+1 >= thetas.size())
    {
      index = thetas.size()-2;
      alpha = 1;
    }

    return (1-alpha) * std::get<2>(thetas[index]) + alpha * std::get<2>(thetas[index+1]);
  }
  else if (useSmoothFunction)
  {
    double tt = t/tend;
    double dds = 1./tend*1./tend * (pow(tt,2)*3*4*35. + pow(tt,3)*4*5*-84. + pow(tt,4)*5*6*70. + pow(tt,5)*6*7*-20.);

    return dds*t01*tend;
  }
  else
  {

    return t02 * 2;
  }
}
