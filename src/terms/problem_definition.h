#pragma once

#include <vector>
#include <array>

using namespace std;

// The following prescribed values are set:

// displacement:
// verschiebung(t) = v2 * t^2 + v1 * t,    with v2=(vx2,vy2), v1=(vx1,vy1)

// angle
// theta0(t) = t02 * t^2 + t01 * t + t00

// friction
// frictionForce(s, h, rho)     computes friction force at s in [0,L]

// global variables that specify the prescribed functions
extern double vx2;
extern double vx1;
extern double vy2;
extern double vy1;
extern bool enableSmoothVelocityProfile;    ///< if the smooth velocity profile with with u(0)=u'(0)=u''(0)=u'''(0) = u(tend)=u'(tend)=u''(tend)=u'''(tend) = 0 should be used

extern double t02;
extern double t01;
extern double t00;

// data needed for using look-up-tables for the prescribed displacements and angles
extern bool useLUTdisplacement;   //< if prescribed displacements from a look up table should be used
extern bool useLUTtheta;          //< if prescribed theta values from a look up table should be used
extern double LUTtimestepWidth;
extern std::vector<std::tuple<double,double,double>> displacementsX;
extern std::vector<std::tuple<double,double,double>> displacementsY;
extern std::vector<std::tuple<double,double,double>> thetas;

// computes the current displacement
array<double,2> displ(double s, double t);

// computes the derivative of the current displacement (=prescribed velocity)
array<double,2> ddispl(double s, double t);

// computes the 2nd derivative of the current displacement (=prescribed acceleration)
array<double,2> dddispl(double s, double t);

// compute the friction force
double frictionForce(double s, double h, double rho);

// current prescribed angle
double theta0(double t);

// current 1st derivative of prescribed angle
double theta0p(double t);

// current 2nd derivative of prescribed angle
double theta0pp(double t);
