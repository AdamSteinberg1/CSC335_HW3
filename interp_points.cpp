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

//constructs a Lagrange Interpolating Polynomial from points
Polynomial lagrange(vector<pair<double, double>> points)
{
  Polynomial result;
  for(int k = 0; k < points.size(); k++)
  {
    //we're going to accumulate a numerator and denominator and then divide them
    Polynomial numerator = {1};
    double denominator = 1;
    for(int i = 0; i < points.size(); i++)
    {
      if(i == k)
        continue;

      double xk = points[k].first;
      double xi = points[i].first;
      numerator = numerator * Polynomial({-xi, 1});
      denominator *= (xk-xi);
    }

    Polynomial l = numerator/denominator;
    result = result + points[k].second * l;
  }
  return result;
}

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

//reads in the data and stores it in a vector of pairs
vector<pair<double, double>> readData(string fileName)
{
  vector<pair<double, double>> data;
  ifstream file(fileName.c_str());
  double x, y;
  while(file >> x >> y)
  {
    data.push_back(make_pair(x,y));
  }
  return data;
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


//writes the file "tmp.plt" which is a gnuplot script to graph:
//  the data in "interp_points.dat"
//  the Lagrange interpolating polynomial p1
//  the Natural Cubic Spline p2, where the boundaries between the cubics are at xValues
void writeGnuPlotScript(Polynomial p1, vector<Polynomial> p2, vector<double> xValues)
{
  ofstream script("tmp.plt");

  script  << "set terminal pngcairo" << endl
          << "set output 'graph.png'" << endl
          << "set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2" << endl
          << "set style line 2 linecolor rgb '#00a000' linetype 1 linewidth 2" << endl
          << "set xrange [0:105]" << endl
          << "set yrange [-10:100]" << endl
          << "p1(x) = " << gnuPrint(p1) << endl
          << "p2(x) = " << gnuPrint(p2, xValues) << endl
          << "plot 'interp_points.dat' pt 7 ps 1, "
          << "p1(x) title 'Lagrange Interpolating Polynomial' with lines linestyle 1, "
          << "p2(x) title 'Natural Cubic Spline' with lines linestyle 2" << endl;
}

int main()
{
  //read in the data from the problem statement
  vector<pair<double, double>> data = readData("interp_points.dat");

  //construct the Lagrange interpolating polynomial
  Polynomial p1 = lagrange(data);

  //make a list of just the x values from the data
  //this is for later when we need the x-values to stitch together the cubic spline
  vector<double> xValues;
  for(auto point : data)
    xValues.push_back(point.first);

  //construct the natual cubic spline
  vector<Polynomial> p2 = cubicSpline(data);

  //write a script for gnuplot
  writeGnuPlotScript(p1, p2, xValues);
  //run the script we just wrote
  system("gnuplot tmp.plt");
  //clean up by deleting the script
  system("rm tmp.plt");

  return 0;
}
