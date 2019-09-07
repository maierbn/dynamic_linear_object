/*
 * main.cpp
 *
 *  Created on: 17.01.2019
 *      Author: mariu
 */
#include <vector>
#include <iostream>
#include <ctime>

#include "terms/problem_definition.h"
#include "utility/utility.h"
#include "numerics/timestepping.h"

using namespace std;

int main() { //TODO Main

  //Main part
  //Eingabe durch den Benutzer

  int n = 8;
  double Rflex = 1.0;
  double L = 1;
  double rho =1.0;
  double mu = 0.3;

  double t0 = 0;
  double tend = 0.5/30.;  // 0.5
  bool ver;
  bool reib;
  bool winkelkontrolle;
  int wdh = 600/30.;
  int iter = 80;
  string path1;
  string path2;

  vector<double> thetaStart(2 * n + 2);
  vector<double> thetaN(n + 1);
  vector<double> thetaNp(n + 1);
  for (int i = 0; i <= n; i++) {
    thetaN[i] = thetaStart[i];
    thetaNp[i] = thetaStart[n + 1 + i];
  }

  time_t t1;

  /*
//1. Fall
  ver = true;
  reib = true;
  winkelkontrolle = false;
  path1 = "Fall1t5.csv";
  path2 = "Fall1t5t.csv";

  vx2 = 0;
  vx1 = 0;
  vy2 = 0.2;
  vy1 = 0;
  t02 = 0;
  t01 = 0;
  t00 = 0;

  t1 = time(0);

  remove(path1.c_str());
  remove(path2.c_str());
  rungeKutta4(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN,
      thetaNp, iter, t0, path1, path2, wdh, tend);
  finish(t1, path1);

*/
  /*
  //1,5. Fall lineare Bewegung
  ver = true;
  reib = true;
  winkelkontrolle = false;
  path1 = "Fall1lint5rho1000R05.csv";
  path2 = "Fall1lint5rho1000R05t.csv";

  vx2 = 0;
  vx1 = 0;
  vy2 = 0;
  vy1 = 1;
  t02 = 0;
  t01 = 0;
  t00 = 0;
  t1 = time(0);

  remove(path1.c_str());
  remove(path2.c_str());
  rungeKutta4(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN,
      thetaNp, iter, t0, path1, path2, wdh, tend);
  finish(t1, path1);


*/
  /*
   //2. Fall
   ver = true;
   reib = false;
   winkelkontrolle = false;
   path1 = "Fall2R10.csv";
   path2 = "Fall2R10t.csv";

   vx2 = 0;
   vx1 = 0;
   vy2 = 18;
   vy1 = 0;
   t02 = 0;
   t01 = 0;
   t00 = 0;
   t1 = time(0);

   remove(path1.c_str());
   remove(path2.c_str());
   rungeKutta4(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN,
   thetaNp, iter, t0, path1, path2, wdh, tend);
   finish(t1, path1);
   */


   //3. Fall
   ver = true;
   reib = true;
   winkelkontrolle = true;
   path1 = "Fall3Reibung.csv";
   path2 = "Fall3Reibungt.csv";

   vx2 = 0;
   vx1 = 6;
   vy2 = 16;
   vy1 = 0;
   t02 = 0;
   t01 = -4.7;
   t00 = 0;
   t1 = time(0);

   remove(path1.c_str());
   remove(path2.c_str());
   rungeKutta4(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN,
   thetaNp, iter, t0, path1, path2, wdh, tend);
   finish(t1, path1);

  /*
   //4. Fall
   ver = false;
   reib = false;
   winkelkontrolle = true;
   path1 = "Fall4t075.csv";
   path2 = "Fall4t075t.csv";
   vx2 = 0;
   vx1 = 0;
   vy2 = 0;
   vy1 = 0;
   t02 = 6;
   t01 = 3;
   t00 = 0;
   t1 = time(0);

   remove(path1.c_str());
   remove(path2.c_str());
   rungeKutta4(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN,
   thetaNp, iter, t0, path1, path2, wdh, tend);
   finish(t1, path1);
   */
  /*
   //5. Fall
   ver = true;
   reib = true;
   winkelkontrolle = true;
   path1 = "Fall5R5n8.csv";
   path2 = "Fall5R5n8t.csv";

   vx2 = 6;
   vx1 = 3;
   vy2 = 0;
   vy1 = 0;
   t02 = 0;
   t01 = 1.6;
   t00 = 0;
   t1 = time(0);

   remove(path1.c_str());
   remove(path2.c_str());
   rungeKutta4(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN,
   thetaNp, iter, t0, path1, path2, wdh, tend);
   finish(t1, path1);

   */
  /*
   //Euler
   string path1 = "ErgebnisseEuler_horizontal_nach_oben_v9_n8_reib.csv";
   string path2 = "ErgebnisseEuler_horizontal_nach_oben_v9_n8_reibt.csv";
   remove(path1.c_str());
   remove(path2.c_str());
   euler(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN, thetaNp,
   iter, t0, path1, path2, wdh, tend);
   finish(t1, path1);*/
  /*

   //RungeKutta
   t1 = time(0);
   string path1 = "ErgebnisseRK4_horizontal_nach_oben_v9_reib.csv";
   string path2 = "ErgebnisseRK4_horizontal_nach_oben_v9_reibt.csv";
   remove(path1.c_str());
   remove(path2.c_str());

   wdh = 600;
   rungeKutta4(n, Rflex, L, rho, mu, ver, reib, winkelkontrolle, thetaN,
   thetaNp, iter, t0, path1, path2, wdh, tend);
   finish(t1, path1);
   */

  /*
   //Test von Alglib
   real_2d_array h;
   vector<double> a(9);
   a[0] = 4;
   a[1] = 5;
   a[2] = -2;
   a[3] = 7;
   a[4] = -1;
   a[5] = 2;
   a[6] = 3;
   a[7] = 1;
   a[8] = 4;
   h.attach_to_ptr(3, 3, a.data());
   vector<double> rhs(3);
   rhs[0] = -14;
   rhs[1] = 42;
   rhs[2] = 28;
   real_1d_array r;
   r.attach_to_ptr(3, rhs.data());
   real_1d_array x;
   vector<double> xlsg(3);
   x.attach_to_ptr(3, xlsg.data());
   int info = 0;
   densesolverreport rep;
   rmatrixsolve(h, 3, r, info, rep, x);
   cout << info << " " << rep.r1 << " " << rep.rinf << endl;
   for (int i = 0; i < 3; ++i) {
   cout << x(i) << endl;
   }
   save(x, path);
   vector<double> lhs(3);
   lhs[0] = -14;
   lhs[1] = 42;
   lhs[2] = 28;
   real_1d_array p;
   vector<vector<double>> A(3);
   A = K(2, 1, 0.5);
   plotteMatrix(A);
   vector<double> c = multipl(A, lhs);
   p.attach_to_ptr(3, c.data());
   save(p, path);

   */

  /*
   // Test von step
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0] = 1;
   thetaN[1] = 2;
   thetaN[2] = 3;
   vector<double> thetaNp;
   thetaNp.resize(3);
   thetaNp[0] = 1;
   thetaNp[1] = 2;
   thetaNp[2] = 3;
   int n = 2;
   int iter = 100;
   double rho = 1;
   double l = 1;
   double h = l / n;
   double mu = 0.3;
   double t = 0;
   double Rflex = 1;
   double tau = 0;
   bool ver = true;
   bool reib = true;
   bool winkelkontrolle = false;
   string path = "";
   time_t t1 = time(0);
   vector<double> test;
   for (int w = 0; w < 1; ++w) {

   test = step(n, Rflex, l, rho, mu, ver, reib,
   winkelkontrolle, thetaN, thetaNp, iter, t, path);
   }
   cout << time(0) - t1 << endl;
   for (int var = 0; var < (n + 1); ++var) {
   cout << test[var] << endl;
   }*/

  /*


   //TEst von uvek,verschvek
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0] = 0;
   thetaN[1] = 5;
   thetaN[2] = 8;
   vector<double> thetaNp;
   thetaNp.resize(3);
   thetaNp[0] = 0;
   thetaNp[1] = 5;
   thetaNp[2] = 1;
   int n = 2;
   int iter = 100;
   double rho = 1;
   double l = 1;
   double h = l / n;
   double mu = 0.3;
   double t = 0.3;
   vector<double> uvek;
   vector<double> verschvek;
   uvek = uVek(thetaN, thetaNp, h, l, mu, t, iter);
   verschvek=verschVek(thetaN, thetaNp,
   h,  l,  t,  rho,  iter);
   for (int var = 0; var < (n+1); ++var) {
   cout<<uvek[var]<<" "<<verschvek[var]<<endl;
   }
   */
  /*
   //Test von mZ bzw anderen Matrizen
   vector<vector<vector<double>>> AMatrix;
   vector<vector<vector<double>>> BMatrix;
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0] = 0;
   thetaN[1] = 5;
   thetaN[2] = 8;
   vector<double> thetaNp;
   thetaNp.resize(3);
   thetaNp[0] = 0;
   thetaNp[1] = 5;
   thetaNp[2] = 1;
   int n = 4;
   int iter = 100;
   //double rho = 1;
   double l = 1;
   double h = l / n;
   double mu = 0.3;
   double t = 0.3;
   double rho = 1;
   double Rflex=1;
   AMatrix = A(1.0, 0.5, n, thetaN, iter);
   BMatrix = B(1.0, 0.5, n, thetaN, iter);
   vector<vector<double>> Test;
   Test = mZ(thetaN, thetaNp, rho,AMatrix,BMatrix);
   plotteMatrix(Test);
   Test= convertMMatrix(Test);
   plotteMatrix(Test);
   */

  /*
   //Test von B
   vector<vector<vector<double>>> BMatrix;
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0] = 0;
   thetaN[1] = 5;
   thetaN[2] = 8;
   int n = 2;
   BMatrix = A(1.0, 0.5, n, thetaN, 200);
   for (int k = 0; k <= n; k++) {
   for (int i = 0; i <= n; i++) {
   for (int j = 0; j <= n; j++) {
   cout << BMatrix[i][j][k] << " ";
   }
   cout << endl;
   }
   cout << endl << endl;
   }
   */
  /*
   //Test von dSinIn
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0] = 0;
   thetaN[1] = 5;
   thetaN[2] = 8;
   cout << SinInt( 1, 0.5, 0.5, thetaN, 200) << endl;
   */

  /*
   //Test von m
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0] = 0;
   thetaN[1] = 5;
   thetaN[2] = 8;
   cout << m(1, 1, 1, 0.5, thetaN, 1, 100) << endl;
   */
  /*
   //Test von CosInt
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0]=0;
   thetaN[1]=5;
   thetaN[2]=8;
   cout<<CosInt(1,0.5,0.5,thetaN,1000)<<endl;
   */

  /*
   //Test von K
   vector<vector<double>> KMatrix;
   int N = 3;
   double RFlex = 1;
   double h = 1.0 / N;
   KMatrix=K(N, RFlex, h);
   plotteMatrix(KMatrix);
   */

  /*
   //Test von Theta
   double erg;
   vector<double> vec;
   vector<double> thetaN;
   thetaN.resize(3);
   thetaN[0]=0;
   thetaN[1]=10;
   thetaN[2]=8;
   vec.resize(101);
   for(int i=0;i<101;i++){
   vec[i]=Theta(thetaN,0.5,i/100.0);
   cout<<(cos(vec[i])*N(1,0.5,i/100.0))<<endl;
   }
   */

  /*
   //Test von floor
   for (int i=0;i<=10;i++){
   cout<<i<<floor(i*0.1);}
   cout<<endl;
   */

  /*
   //Test von N
   double erg2;
   N(0,0.5,0.25,erg2);
   cout<<erg2<<endl;
   */

  /*
   //Test der Integralfunktion
   double (*f)(double)=funk;
   double erg;
   Integral(f,0,3.1415926535,10,erg);
   std::cout << erg<<endl ;
   */

  return 0;
}