#pragma once

#include <vector>

using namespace std;

// linear hat function, Eq. (3)
double N(int i, double h, double hInv, double s);

// discretized angle
double Theta(const vector<double> &thetaN, double h, double hInv, double s);

// compute integral over cosine of theta, multiplied with ith basis function, Eq. (7.1)
double CosInt(int i, double s, double h, double hInv, const vector<double> &thetaN);

// compute integral over sine of theta, multiplied with ith basis function, Eq. (7.2)
double SinInt(int i, double s, double h, double hInv, const vector<double> &thetaN);

// entry of mass matrix, m_(i,k), Eq. (6.1)
double m(int i, int k, double rho, double h, double hInv, const vector<double> &thetaN,
    double l);

// dSinInt/dTheta
double dSinInt(int i, int r, double s, double h, double hInv, const vector<double> &thetaN);

// dCosInt/dTheta
double dCosInt(int i, int r, double s, double h, double hInv, const vector<double> &thetaN);

// compute the integral over dSinInt(i,r)*SinInt(k) for reuse, for all i,k,r
void computeMatrixA(double l, double h, double hInv, int n,
    vector<double> &thetaN);

// compute the integral over dCosInt(i,r)*CosInt(k) for reuse, for all i,k,r
void computeMatrixB(double l, double h, double hInv, int n,
    vector<double> &thetaN);

// compute the stiffness matrix K
void computeMatrixK(int n, double RFlex, double h, double hInv);

// compute the derivative of m(i,k) w.r.t Theta_r
double dm(int i, int k, int r, double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

// compute Z_r,k, Eq. (6.2)
double Z(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

// compute Y_r,k, Eq. (6.3)
double Y(int r, int k, vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

// compute mass matrix M with M_ij = m(i,j,...)
void computeMatrixM(double rho, double h, double hInv, vector<double> &thetaN,
    double l);

// compute matrix Z, with Z_ij = Z(i,j,...)
void computeMatrixZ(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

// compute matrix Y, with Y_ij = Y(i,j,...)
void computeMatrixY(vector<double> &thetaN, vector<double> &thetaNp,
    double rho, vector<vector<vector<double>>> &A,
    vector<vector<vector<double>>> &B);

// compute entry V_rj, Eq. (11)
double V(int r, int j, vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double mu, double t, double rho);

// compute matrix V, with V_ij = V(i,j,...)
void computeMatrixV(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double l, double mu, double t, double rho);

// compute entry u_r, Eq. (12)
double u(int r, vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double mu, double t, double rho);

// compute the vector u with u_i= u(i,...)
vector<double> uVector(vector<double> &thetaN, vector<double> &thetaNp, double h,
    double hInv, double l, double mu, double t, double rho);

// compute vector w, Eq. (9)
vector<double> wVector(vector<double> &thetaN, vector<double> &thetaNp,
    double h, double hInv, double l, double t, double rho);

// adjust mass matrix M, such that the first angle, \theta_0 is given by the prescribed angle function in problem_definition.cpp
void convertMMatrix();
