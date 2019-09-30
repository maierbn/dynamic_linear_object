#include "terms/problem_definition.h"

#include <tuple>
#include <cmath>

#include "terms/globals.h"

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
bool enableSmoothVelocityProfile;
std::vector<std::tuple<double,double,double>> displacementsX;
std::vector<std::tuple<double,double,double>> displacementsY;
std::vector<std::tuple<double,double,double>> thetas;

using std::pow;

// computes the current displacement
array<double,2> displ(double s, double t)
{
  array<double,2> out;

  if (useLUTdisplacement)    // use recorded data
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // when there is no more data, repeat the last recorded value
    if (index+1 >= displacementsX.size())
    {
      index = displacementsX.size()-2;
      alpha = 1;
    }

    out[0] = (1-alpha) * std::get<0>(displacementsX[index]) + alpha * std::get<0>(displacementsX[index+1]);
    out[1] = (1-alpha) * std::get<0>(displacementsY[index]) + alpha * std::get<0>(displacementsY[index+1]);
  }
  else if (enableSmoothVelocityProfile)    // superimpose the smooth velocity profile
  {
    double tt = t/tend;
    double s = pow(tt,4)*35. + pow(tt,5)*-84. + pow(tt,6)*70. + pow(tt,7)*-20.;

    out[0] = s*vx1*tend;
    out[1] = s*vy1*tend;
  }
  else   // normal quadratic function
  {
    out[0] = vx2 * t * t + vx1 * t;
    out[1] = vy2 * t * t + vy1 * t;
  }
  return out;
}

// computes the derivative of the current displacement (=prescribed velocity)
array<double,2> ddispl(double s, double t)
{
  array<double,2> out;

  if (useLUTdisplacement)    // use recorded data
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // when there is no more data, repeat the last recorded value
    if (index+1 >= displacementsX.size())
    {
      index = displacementsX.size()-2;
      alpha = 1;
    }

    out[0] = (1-alpha) * std::get<1>(displacementsX[index]) + alpha * std::get<1>(displacementsX[index+1]);
    out[1] = (1-alpha) * std::get<1>(displacementsY[index]) + alpha * std::get<1>(displacementsY[index+1]);
  }
  else if (enableSmoothVelocityProfile)    // superimpose the smooth velocity profile
  {
    double tt = t/tend;
    double ds = 1./tend * (pow(tt,3)*4*35. + pow(tt,4)*5*-84. + pow(tt,5)*6*70. + pow(tt,6)*7*-20.);

    out[0] = ds*vx1*tend;
    out[1] = ds*vy1*tend;
  }
  else   // normal quadratic function
  {
    out[0] = vx2 * 2 * t + vx1;
    out[1] = vy2 * 2 * t + vy1;
  }
  return out;
}

// computes the 2nd derivative of the current displacement (=prescribed acceleration)
array<double,2> dddispl(double s, double t)
{
  array<double,2> out;

  if (useLUTdisplacement)    // use recorded data
  {
    double timeStepNo = t/LUTtimestepWidth;
    int index = int(timeStepNo);
    double alpha = timeStepNo - index;

    // when there is no more data, repeat the last recorded value
    if (index+1 >= displacementsX.size())
    {
      index = displacementsX.size()-2;
      alpha = 1;
    }

    out[0] = (1-alpha) * std::get<2>(displacementsX[index]) + alpha * std::get<2>(displacementsX[index+1]);
    out[1] = (1-alpha) * std::get<2>(displacementsY[index]) + alpha * std::get<2>(displacementsY[index+1]);
  }
  else if (enableSmoothVelocityProfile)    // superimpose the smooth velocity profile
  {
    double tt = t/tend;
    double dds = 1./tend*1./tend * (pow(tt,2)*3*4*35. + pow(tt,3)*4*5*-84. + pow(tt,4)*5*6*70. + pow(tt,5)*6*7*-20.);

    out[0] = dds*vx1*tend;
    out[1] = dds*vy1*tend;
  }
  else   // normal quadratic function
  {
    out[0] = vx2 * 2;
    out[1] = vy2 * 2;
  }
  return out;
}

double frictionForce(double s, double h, double rho)
{
  //compute the friction force at coordinate s
  return 9.81 * rho * h;
}

// current prescribed angle
double theta0(double t)
{
  if (useLUTtheta)    // use recorded data
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
  else if (enableSmoothVelocityProfile)    // superimpose the smooth velocity profile
  {
    double tt = t/tend;
    double s = pow(tt,4)*35. + pow(tt,5)*-84. + pow(tt,6)*70. + pow(tt,7)*-20.;

    return s*t01*tend + t00;
  }
  else   // normal quadratic function
  {
    return t02 * t * t + t01 * t + t00;
  }
}

// current 1st derivative of prescribed angle
double theta0p(double t)
{
  if (useLUTtheta)    // use recorded data
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
  else if (enableSmoothVelocityProfile)    // superimpose the smooth velocity profile
  {
    double tt = t/tend;
    double ds = 1./tend * (pow(tt,3)*4*35. + pow(tt,4)*5*-84. + pow(tt,5)*6*70. + pow(tt,6)*7*-20.);

    return ds*t01*tend;
  }
  else   // normal quadratic function
  {
    return t02 * 2 * t + t01;
  }
}

// current 2nd derivative of prescribed angle
double theta0pp(double t)
{
  if (useLUTtheta)    // use recorded data
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
  else if (enableSmoothVelocityProfile)    // superimpose the smooth velocity profile
  {
    double tt = t/tend;
    double dds = 1./tend*1./tend * (pow(tt,2)*3*4*35. + pow(tt,3)*4*5*-84. + pow(tt,4)*5*6*70. + pow(tt,5)*6*7*-20.);

    return dds*t01*tend;
  }
  else   // normal quadratic function
  {
    return t02 * 2;
  }
}
