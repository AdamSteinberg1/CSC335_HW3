#include "Polynomial.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;


//constructs a clamped cubic spline from points and fpo and fpn
//fpo is the derivative of f at x_0
//fpn is the derivative of f at x_n
//each cubic we make is defined by four constants: a_i, b_i, c_i, d_i
vector<Polynomial> clampedCubicSpline(vector<pair<double, double>> points, double fpo, double fpn)
{
  //how many cubics we're going to make
  int n = points.size()-1;

  //h contains the difference between consecutive x-values
  vec h(n, fill::zeros);
  for(int i = 0; i < n; i++)
    h(i) = points[i+1].first - points[i].first;

  //alpha is a vector representing the constants on the right side of our system of linear equations
  vec alpha(n+1, fill::zeros);
  alpha(0) = 3*(points[1].second - points[0].second)/h(0) - 3*fpo;
  alpha(n) = 3*fpn - 3*(points[n].second - points[n-1].second)/h(n-1);
  for(int i = 1; i < n; i++)
    alpha(i) = 3/h[i]*(points[i+1].second - points[i].second) - 3/h[i-1]*(points[i].second - points[i-1].second);

  //A is a matrix representing the coefficients in our system of linear equations
  mat A(n+1, n+1, fill::zeros);
  A(0,0) = 2*h(0);
  A(0,1) = h(0);
  A(n,n-1) = h(n-1);
  A(n,n) = 2*h(n-1);
  for(int i = 1; i < n; i++)
  {
    A(i,i-1) = h[i-1];
    A(i,i) = 2*(h[i-1]+h[i]);
    A(i,i+1) = h[i];
  }

  //solve A*c = alpha for c
  //c contains all our c_i constants
  vec c = solve(A, alpha);

  vector<Polynomial> result;
  for(int i = 0; i < n; i++)
  {
    //calculate all the constants that define the ith cubic
    double a_i = points[i].second;
    double b_i = (points[i+1].second - points[i].second)/h[i] - h[i]*(c(i+1)+2*c(i))/3;
    double c_i = c(i);
    double d_i = (c(i+1)-c(i))/(3*h[i]);
    //the x-values of the input points
    double x_i = points[i].first;

    //construct the ith cubic p
    Polynomial diff = {-x_i, 1}; // (x - x_i)
    Polynomial p = a_i + b_i*diff + c_i*diff.power(2) + d_i*diff.power(3);
    result.push_back(p);
  }
  return result;
}


//returns a string representing a polynomial in gnuplot's format
string gnuPrint(Polynomial p)
{
  stringstream result;
  result.precision(numeric_limits<long double>::max_digits10);
  int n = p.getDegree();
  for (int i =0; i <= n; i++)
  {
    if(p[i] == 0)
      continue;

    result << p[i];
    if (i != 0)
      result << "*x";

    if(i > 1)
      result << "**" << i;

    if(i < n)
      result << " + ";
  }
  return result.str();
}

//graphs the three
void graph(vector<Polynomial> spline1, vector<pair<double, double>> points1, vector<Polynomial> spline2, vector<pair<double, double>> points2, vector<Polynomial> spline3, vector<pair<double, double>> points3)
{
  ofstream script("tmp.plt");

  script  << "set terminal pngcairo" << endl
          << "set output 'nobleBeast.png'" << endl
          << "set xrange [0:30]" << endl
          << "set yrange [0:30]" << endl
          << "set samples 200" << endl
          << "p(x) = ";

  for(int i = 0; i < spline1.size(); i++)
  {
    script << "x<" << points1[i+1].first << " ? " << gnuPrint(spline1[i]) << " : ";
  }
  for(int i = 0; i < spline2.size(); i++)
  {
    script << "x<" << points2[i+1].first << " ? " << gnuPrint(spline2[i]) << " : ";
  }
  for(int i = 0; i < spline3.size()-1; i++)
  {
    script << "x<" << points3[i+1].first << " ? " << gnuPrint(spline3[i]) << " : ";
  }
  script << gnuPrint(spline3.back()) << endl;

  script  << "plot p(x) notitle" <<  endl;
  script.close();

  system("gnuplot tmp.plt");
  system("rm tmp.plt");
}

int main()
{
  vector<pair<double, double>> curve1 = {
    {1, 3},
    {2, 3.7},
    {5, 3.9},
    {6, 4.2},
    {7, 5.7},
    {8, 6.6},
    {10, 7.1},
    {13, 6.7},
    {17, 4.5},
  };

  vector<pair<double, double>> curve2 = {
    {17, 4.5},
    {20, 7.0},
    {23, 6.1},
    {24, 5.6},
    {25, 5.8},
    {27, 5.2},
    {27.7, 4.1}
  };

  vector<pair<double, double>> curve3 = {
    {27.7, 4.1},
    {28, 4.3},
    {29, 4.1},
    {30, 3.0}
  };

  //construct the clamped cubic splines
  vector<Polynomial> p1 = clampedCubicSpline(curve1, 1.0, -2.0/3.0);
  vector<Polynomial> p2 = clampedCubicSpline(curve2, 3.0, -4.0);
  vector<Polynomial> p3 = clampedCubicSpline(curve3, 1.0/3.0, -1.5);

  int i = 1;
  for (auto & p : {p1,p2,p3})
  {
    cout << "Clamped Cubic Spline for Curve " << i << ":" << endl;
    for(int j = 0; j < p.size(); j++)
    {
      cout << "S" << j << "(x) = " << p[j] << endl;
    }
    cout << endl;
    i++;
  }

  graph(p1, curve1, p2, curve2, p3, curve3);

  return 0;
}
