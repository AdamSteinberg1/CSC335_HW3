#include "Polynomial.h"
#include <vector>
#include <iostream>
using namespace std;

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

  //part a
  Polynomial p1 = lagrange(sample1);
  Polynomial p2 = lagrange(sample2);

  cout << "p1 =  "<<p1<<endl;
  cout << "p2 =  "<<p2<<endl;

  //part b
  Polynomial dp1 = p1.derivative();
  Polynomial dp2 = p2.derivative();
  cout << "derivative of p1 = " << dp1 << endl;
  cout << "derivative of p2 = " << dp2 << endl;

  //initial guesses are picked because they are the largest point's in the tabulated data
  double criticalValue1 = dp1.root(10);
  double criticalValue2 = dp2.root(10);

  cout << "p1 has a maximum of " << p1.evaluate(criticalValue1) << " at " << criticalValue1 << endl;
  cout << "p2 has a maximum of " << p2.evaluate(criticalValue1) << " at " << criticalValue2 << endl;

  return 0;
}
