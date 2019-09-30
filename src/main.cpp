/*
 * main.cpp
 *
 *  Created on: 17.01.2019, revised 30.09.2019
 *      Author: Marius Stach, Benjamin Maier
 */
#include <vector>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <sstream>

#include "terms/problem_definition.h"
#include "utility/utility.h"
#include "numerics/timestepping.h"
#include "numerics/integral.h"
#include "terms/globals.h"

using namespace std;

/** All parameters that are needed to specify a simulation.
 */
struct SimulationConfiguration
{
  bool enablePrescribedDisplacement = false;  ///< if there are prescribed displacements of the start point P0
  bool enableFriction = true;                 ///< if friction is enabled, if set to false, value of mu is not considered
  bool enablePrescribedAngle = true;          ///< if there is a prescribed angle at the start point P0, \hat\theta
  std::string outputFilename;                 ///< file name of the output file

  // parameters of the prescribed displacement:
  // displacement(t) = v2 * t^2 + v1 * t,    where v2 = (vx2,vy2), v1 = (vx1,vy1)
  double vx1;
  double vy1;
  double vx2;
  double vy2;
  bool enableSmoothVelocityProfile = false;   ///< if the prescribed velocity follows a smooth velocity profile with u(0)=u'(0)=u''(0)=u'''(0) = u(tend)=u'(tend)=u''(tend)=u'''(tend) = 0

  // parameters of the prescribed angle:
  // theta0(t) = t02 * t^2 + t01 * t + t00
  double t02;
  double t01;
  double t00;
  double initialAngle;                        ///< the initial angle at which the object is positioned

  double rho;      ///< linear density, mass per cross section area
  double L;        ///< length of the linear object
  double Rflex;    ///< flexural rigidity
  double mu;       ///< constant friction coefficient

  int n;           ///< number of Finite Elements
  double dt;       ///< time step width
  double tend;     ///< end time of the simulation

  void runSimulation()
  {
    // start duration measurement
    std::chrono::time_point<std::chrono::steady_clock> t1 = std::chrono::steady_clock::now();

    // print information
    double h = L/n;
    cout << "output file: " << outputFilename << endl
      << "n: " << n << ", dt: " << dt << ", tend: " << tend << ", L: " << L << " m, Rflex: " << Rflex << " N*m^2, rho: "
      << rho << " kg/m" << ", h: " << h << ", rho*h: " << rho*h << " kg, ReibF: " << frictionForce(L/2., h, rho) << endl;

    // delete existing output files
    remove(outputFilename.c_str());

    // save options
    save(n, Rflex, L, rho, mu, enablePrescribedDisplacement, enableFriction, enablePrescribedAngle, outputFilename);

    // initialize global values for problem definition
    ::vx1 = vx1;
    ::vy1 = vy1;
    ::vx2 = vx2;
    ::vy2 = vy2;
    ::t00 = t00;
    ::t01 = t01;
    ::t02 = t02;
    ::enableSmoothVelocityProfile = enableSmoothVelocityProfile;
    ::tend = tend;

    // initialize variables
    vector<double> thetaN(n + 1);     ///< vector of degrees of freedom, \Theta_n
    vector<double> thetaNp(n + 1);    ///< first derivative of degrees of freedom, \omega_n = \Theta_n'
    for (int i = 0; i < n+1; i++)
    {
      thetaN[i] = initialAngle;
      thetaNp[i] = 0.0;
    }

    // run time stepping scheme
    rungeKutta4(n, Rflex, L, rho, mu, enablePrescribedDisplacement, enableFriction, enablePrescribedAngle, thetaN,
      thetaNp, outputFilename, dt, tend);

    // finalize output file and measure duration
    finalize(t1, outputFilename);
  }
};

int main()
{
  // define, which scenario should be run
  std::string scenario = "experiment4";       // set this to experiment1 - experiment4 to reproduce the simulations of the paper

  std::cout << "Scenario " << scenario << " ";

  SimulationConfiguration s;

  if (scenario == "experiment1")
  {
    // Experiment #1: Demonstrate effects of friction

    // problem definition
    s.enablePrescribedDisplacement = true;
    s.enableFriction = true;     // false
    s.enablePrescribedAngle = false;
    s.outputFilename = "experiment1.csv";

    // prescribed displacement:
    // displacement(t) = v2 * t^2 + v1 * t,    where v2 = (vx2,vy2), v1 = (vx1,vy1)
    s.vx2 = 0;
    s.vx1 = 0;
    s.vy2 = 1;
    s.vy1 = 0;

    // prescribed angle:
    // theta0(t) = t02 * t^2 + t01 * t + t00
    s.t02 = 0;
    s.t01 = 0;
    s.t00 = 0;

    // material parameters
    s.Rflex = 1.0;             // flexural rigidity [Nm^2]
    s.L = 1;                   // length [m]
    s.rho = 1.0;               // linear density (density * cross section area)  [kg/m^3]*[m^2] = [kg/m]
    s.mu = 1.0;                // sliding friction coefficient [-]

    // numeric parameters
    s.n = 3;                   // number of Finite Elements [-]
    s.dt = 1e-2;               // time step width [s]
    s.tend = 2;                // end time [s]

    s.runSimulation();
  }
  else if (scenario == "experiment2")
  {
    // Experiment #2: Demonstrate effects of flexural rigidity

    // problem definition
    s.enablePrescribedDisplacement = true;
    s.enableFriction = true;
    s.enablePrescribedAngle = true;
    s.outputFilename = "experiment2.csv";

    // prescribed displacement:
    // displacement(t) = v2 * t^2 + v1 * t,    where v2 = (vx2,vy2), v1 = (vx1,vy1)
    s.vx2 = 6;
    s.vx1 = 3;
    s.vy2 = 0;
    s.vy1 = 0;

    // prescribed angle:
    // theta0(t) = t02 * t^2 + t01 * t + t00
    s.t02 = 0;
    s.t01 = M_PI/2.;
    s.t00 = 0;

    s.initialAngle = 0;

    // material parameters
    s.Rflex = 1.0; // 1.0 or 5.0      flexural rigidity [Nm^2]
    s.L = 1;                   // length [m]
    s.rho = 1.0;               // linear density (density * cross section area)  [kg/m^3]*[m^2] = [kg/m]
    s.mu = 0.3;                // sliding friction coefficient [-]

    // numeric parameters
    s.n = 3;                   // number of Finite Elements [-]
    s.dt = 1e-2;               // time step width [s]
    s.tend = 1.0;              // end time [s]

    s.runSimulation();
  }
  else if (scenario == "experiment3")
  {
    // Experiment #3: Rotation of white cable

    // problem definition
    s.enablePrescribedDisplacement = false;
    s.enableFriction = true;
    s.enablePrescribedAngle = true;
    s.outputFilename = "experiment3.csv";

    // prescribed displacement:
    // displacement(t) = v2 * t^2 + v1 * t,    where v2 = (vx2,vy2), v1 = (vx1,vy1)
    s.vx1 = 0;
    s.vy1 = 0;
    s.vx2 = 0;
    s.vy2 = 0;

    // prescribed angle:
    // theta0(t) = t02 * t^2 + t01 * t + t00
    s.t02 = 0;
    s.t01 = M_PI / 10.0;
    s.t00 = M_PI;
    s.initialAngle = M_PI;
    s.enableSmoothVelocityProfile = false;

    // material parameters
    double length = 44.5e-2;  // [m]
    double weight = 8e-3;    // [kg]

    s.rho = weight / length;   // linear density (density * cross section area)  [kg/m^3]*[m^2] = [kg/m]
    s.L = 40e-2;               // length [m]
    s.Rflex = 9e-4;            // flexural rigidity [Nm^2]
    s.mu = 0.3;                // sliding friction coefficient [-]

    // numeric parameters
    s.n = 3;                   // number of Finite Elements [-]
    s.tend = 10;               // end time [s]
    s.dt = 5e-3;               // time step width [s]

    s.runSimulation();
  }
  else if (scenario == "experiment4")
  {
    // Experiment #4: Movement and rotation of flexible tube

    // problem definition
    s.enablePrescribedDisplacement = true;
    s.enableFriction = true;
    s.enablePrescribedAngle = true;
    s.outputFilename = "experiment4.csv";

    // prescribed displacement:
    // displacement(t) = v2 * t^2 + v1 * t,    where v2 = (vx2,vy2), v1 = (vx1,vy1)
    s.vx1 = 0.4 / 4;
    s.vy1 = -0.1 / 4;
    s.vx2 = 0;
    s.vy2 = 0;

    // prescribed angle:
    // theta0(t) = t02 * t^2 + t01 * t + t00
    s.t02 = 0;
    s.t01 = M_PI_2*1.5 / 4.0;
    s.t00 = M_PI;
    s.initialAngle = M_PI;
    s.enableSmoothVelocityProfile = true;

    // material parameters
    double length = 80e-2;  // [m]
    double weight = 120e-3;    // [kg]

    s.rho = weight / length;    // linear density (density * cross section area)  [kg/m^3]*[m^2] = [kg/m]
    s.L = length;               // length [m]
    s.Rflex = 5e-2;             // flexural rigidity [Nm^2]
    s.mu = 0.6;                 // sliding friction coefficient [-]

    // numeric parameters
    s.n = 3;                   // number of Finite Elements [-]
    s.dt = 1e-2;               // time step width [s]
    s.tend = 4;                // end time [s]

    s.runSimulation();
  }

  // call python visualization script with output file name
  cout << "Now, for visualization call:       ./plot.py " << s.outputFilename << endl;

  std::stringstream command;
  command << "./plot.py " << s.outputFilename;
  int ret = system(command.str().c_str()); ret++;

  return EXIT_SUCCESS;
}
