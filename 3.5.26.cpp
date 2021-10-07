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


//constructs a natural cubic spline from points
//each cubic we make is defined by four constants: a_i, b_i, c_i, d_i
vector<Polynomial> cubicSpline(vector<pair<double, double>> points)
{
  //how many cubics we're going to make
  int n = points.size()-1;

  //h contains the difference between consecutive x-values
  vec h(n, fill::zeros);
  for(int i = 0; i < n; i++)
    h(i) = points[i+1].first - points[i].first;

  //alpha is a vector representing the constants on the right side of our system of linear equations
  vec alpha(n+1, fill::zeros);
  for(int i = 1; i < n; i++)
    alpha(i) = 3/h[i]*(points[i+1].second - points[i].second) - 3/h[i-1]*(points[i].second - points[i-1].second);

  //A is a matrix representing the coefficients in our system of linear equations
  mat A(n+1, n+1, fill::zeros);
  A(0,0) = 1.0;
  A(n,n) = 1.0;
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

int main()
{
  vector<pair<double, double>> data = {
    {0, 0},
    {0.25, 23.04},
    {0.5, 47.37},
    {1.0, 97.45},
    {1.25, 123.66},
  };

  //part a
  //construct the natual cubic spline
  vector<Polynomial> p = cubicSpline(data);

  cout << "Natural Cubic Spline:" << endl;
  for(int i = 0; i < p.size(); i++)
    cout << "S" << i << "(x) = " << p[i] << endl;
  cout << endl;

  //part b
  //predict time when distance is 0.75
  //use 3rd cubic polynomial because 0.5 < 0.75 < 1.0
  cout << "time at 0.75 miles = " << p[2].evaluate(0.75) << " seconds" << endl;

  //part c
  //find speed at start and end
  double startspeed = 3600/p.front().derivative().evaluate(0);
  double endspeed = 3600/p.back().derivative().evaluate(1.25);

  cout << "California Chrome's starting speed = " << startspeed << " mph" << endl;
  cout << "California Chrome's speed at finish line = " << endspeed << " mph" << endl;

  return 0;
}
