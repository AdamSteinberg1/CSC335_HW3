#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <ostream>

class Polynomial
{
  private:
    std::vector<long double> coefficients;

    //returns the result of distributing ax^n to this polynomial
    Polynomial distribute(long double a, int n) const
    {
      std::vector<long double> resultCoefficients;
      for(int i = 0; i < n; i++)
      {
        resultCoefficients.push_back(0.0);
      }

      for(int i = 0; i < coefficients.size(); i++)
      {
        resultCoefficients.push_back(a*coefficients[i]);
      }
      return Polynomial(resultCoefficients);
    }

  public:
    //default constructor creates the polynomial 0x^0
    Polynomial()
    {
      coefficients = {0.0};
    }

    //constructs a polynomial where args are the coefficients listed from lowest to highest degree
    Polynomial(std::initializer_list<long double> args)
    {
      coefficients = args;
    }

    //constructs a polynomial where args are the coefficients listed from lowest to highest degree
    Polynomial(std::vector<long double> args)
    {
      coefficients = args;
    }

    std::vector<long double> getCoefficients() const
    {
      return coefficients;
    }

    long double operator[](int i) const
    {
      return coefficients[i];
    }

    int getDegree() const
    {
      return coefficients.size();
    }

    //returns the polynomial evaluated at a specific x value
    long double evaluate(long double x)
    {
      long double result = 0;
      for(int i = 0; i < getDegree(); i++)
      {
        result += coefficients[i]*pow(x,i);
      }
      return result;
    }

    friend Polynomial operator*(Polynomial a, Polynomial b);
};

Polynomial operator+(Polynomial a, Polynomial b)
{
  int resultDegree = std::max(a.getDegree(), b.getDegree());
  std::vector<long double> resultCoefficients(resultDegree, 0);

  for(int i = 0; i < a.getDegree(); i++)
  {
    resultCoefficients[i] += a[i];
  }

  for(int i = 0; i < b.getDegree(); i++)
  {
    resultCoefficients[i] += b[i];
  }
  return Polynomial(resultCoefficients);
}

Polynomial operator*(long double scalar, Polynomial p)
{
  std::vector<long double> resultCoefficients;
  int n = p.getDegree();
  for (int i = 0; i < n; i++)
  {
    resultCoefficients.push_back(scalar * p[i]);
  }
  return Polynomial(resultCoefficients);
}

Polynomial operator/(Polynomial p, long double scalar)
{
  std::vector<long double> resultCoefficients;
  int n = p.getDegree();
  for (int i = 0; i < n; i++)
  {
    resultCoefficients.push_back(p[i]/scalar);
  }
  return Polynomial(resultCoefficients);
}

Polynomial operator-(Polynomial a, Polynomial b)
{
  return a+(-1*b);
}

Polynomial operator*(Polynomial a, Polynomial b)
{
  Polynomial result;
  for(int i = 0; i < a.getDegree(); i++)
  {
      if(a[i]!=0)
      {
        result = result + b.distribute(a[i], i);
      }
  }

  return result;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& p)
{
  int n = p.getDegree();
  for (int i =0; i < n; i++)
  {
    if(p[i] == 0)
      continue;

    os << p[i];
    if (i != 0)
      os << "x";

    if(i > 1)
      os << "^" << i;

    if(i < n-1)
      os << " + ";
  }
  return os;
}

Polynomial lagrange(std::vector<std::pair<long double, long double>> points)
{
  Polynomial result;
  for(int k = 0; k < points.size(); k++)
  {
    Polynomial numerator = {1};
    long double denominator = 1;
    for(int i = 0; i < points.size(); i++)
    {
      if(i == k)
        continue;

      auto xk = points[k].first;
      auto xi = points[i].first;
      numerator = numerator * Polynomial({-xi, 1});
      denominator *= (xk-xi);
    }

    Polynomial l = numerator/denominator;
    result = result + points[k].second * l;
  }
  return result;
}
