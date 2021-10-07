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

//returns a string representing a piecwise function where p is the list of
//polynomials that are stitched together at the x values
string gnuPrint(vector<Polynomial> p, vector<double> xValues)
{
  stringstream result;
  for(int i = 0; i < p.size()-1; i++)
  {
    result << "x<" << xValues[i+1] << " ? " << gnuPrint(p[i]) << " : ";
  }
  result << gnuPrint(p.back());

  return result.str();
}

void graphSpline(vector<Polynomial> spline, vector<pair<double, double>> points, string fileName)
{
  vector<double> xValues;
  for(auto point : points)
    xValues.push_back(point.first);

  ofstream script("tmp.plt");

  script  << "set terminal pngcairo" << endl
          << "set output '"<< fileName <<".png'" << endl
          << "p(x) = " << gnuPrint(spline, xValues) << endl
          << "plot p(x) title '"<< fileName <<"', '-' notitle" << endl;

  for (auto point : points)
    script << point.first << " " << point.second << endl;

  script.close();

  system("gnuplot tmp.plt");
  system("rm tmp.plt");
}

int main()
{
  vector<pair<double, double>> sample1 = {
    {0, 6.67},
    {6, 17.33},
    {10, 42.67},
    {13, 37.33},
    {17, 30.10},
    {20, 29.31},
    {28, 28.74},
  };

  vector<pair<double, double>> sample2 = {
    {0, 6.67},
    {6, 16.11},
    {10, 18.89},
    {13, 15.00},
    {17, 10.56},
    {20, 9.44},
    {28, 8.89},
  };

  //construct the natual cubic spline
  vector<Polynomial> p1 = cubicSpline(sample1);
  vector<Polynomial> p2 = cubicSpline(sample2);

  cout << "Natural Cubic Spline for Sample 1:" << endl;
  for(int i = 0; i < p1.size(); i++)
  {
    cout << "S" << i << "(x) = " << p1[i] << endl;
  }
  cout << endl;

  cout << "Natural Cubic Spline for Sample 2:" << endl;
  for(int i = 0; i < p2.size(); i++)
  {
    cout << "S" << i << "(x) = " << p2[i] << endl;
  }
  cout << endl;

  graphSpline(p1, sample1, "sample1");
  graphSpline(p2, sample2, "sample2");

  //from the graph of sample 1 we know the maximum is in
  //the 3rd cubic polynomial near x=10
  double maxLocation1 = p1[2].derivative().root(10);
  cout << "Maximum of sample 1 spline is (" << maxLocation1 <<", "<<p1[2].evaluate(maxLocation1) <<")"<< endl;

  //from the graph of sample 1 we know the maximum is in
  //the 2nd cubic polynomial near x=8
  double maxLocation2 = p2[1].derivative().root(8);
  cout << "Maximum of sample 2 spline is (" << maxLocation2 <<", "<<p2[1].evaluate(maxLocation1) <<")"<< endl;

  return 0;
}
