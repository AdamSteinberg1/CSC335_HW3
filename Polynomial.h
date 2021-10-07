//class to represent a polynomial
//Written by Adam Steinberg
#pragma once
#include <vector>
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
      //ensure that the coefficients vector isn't bigger than it needs to be
      while(coefficients.back() == 0)
      {
        coefficients.pop_back();
      }
    }

    //constructs a polynomial where args are the coefficients listed from lowest to highest degree
    Polynomial(std::vector<long double> args)
    {
      coefficients = args;
      //ensure that the coefficients vector isn't bigger than it needs to be
      while(coefficients.back() == 0)
      {
        coefficients.pop_back();
      }
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
      return coefficients.size() - 1;
    }

    //returns the polynomial evaluated at a specific x value
    long double evaluate(long double x)
    {
      long double result = 0;
      for(int i = 0; i <= getDegree(); i++)
      {
        result += coefficients[i]*pow(x,i);
      }
      return result;
    }

    Polynomial derivative()
    {
      //derivative of a constant is 0
      if(getDegree() == 0)
        return Polynomial();

      //shift all the coefficients down by a degree
      //the x^0 term is removed because derivative of a constant is 0
      std::vector<long double> resultCoefficients(coefficients.begin()+1, coefficients.end());
      for(int i = 0; i < resultCoefficients.size(); i++)
      {
        resultCoefficients[i] *= i+1;
      }
      return Polynomial(resultCoefficients);
    }

    //finds a root of the polynomial using Newton's Method with an initial approximation of p0
    long double root(long double p0)
    {
      const int MAX_ITERATIONS = 10000;
      const long double TOLERANCE = 0.00000000001;

      Polynomial fPrime = derivative();

      for (int i = 1; i <= MAX_ITERATIONS; i++)
      {
        long double p = p0 - evaluate(p0)/fPrime.evaluate(p0);
        if(fabs(p - p0) < TOLERANCE)
        {
          return p;
        }
        p0 = p;
      }
      return p0;
    }

    Polynomial power(int n)
    {
      if (n == 0)
        return Polynomial({1});
      if (n == 1)
        return (*this);
      if (n % 2 == 0)
      {
          Polynomial m = this->power(n / 2);
          return m * m;
      }
      else
        return (*this) * (this->power(n - 1));
    }

    friend Polynomial operator*(Polynomial a, Polynomial b);
};

Polynomial operator+(Polynomial a, Polynomial b)
{
  int resultDegree = std::max(a.getDegree(), b.getDegree());
  std::vector<long double> resultCoefficients(resultDegree+1, 0);

  for(int i = 0; i <= a.getDegree(); i++)
  {
    resultCoefficients[i] += a[i];
  }

  for(int i = 0; i <= b.getDegree(); i++)
  {
    resultCoefficients[i] += b[i];
  }
  return Polynomial(resultCoefficients);
}

Polynomial operator*(long double scalar, Polynomial p)
{
  std::vector<long double> resultCoefficients;
  int n = p.getDegree();
  for (int i = 0; i <= n; i++)
  {
    resultCoefficients.push_back(scalar * p[i]);
  }
  return Polynomial(resultCoefficients);
}

Polynomial operator+(long double scalar, Polynomial p)
{
  std::vector<long double> resultCoefficients = p.getCoefficients();
  resultCoefficients[0] += scalar;
  return Polynomial(resultCoefficients);
}

Polynomial operator/(Polynomial p, long double scalar)
{
  std::vector<long double> resultCoefficients;
  int n = p.getDegree();
  for (int i = 0; i <= n; i++)
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
  for(int i = 0; i <= a.getDegree(); i++)
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
  for (int i =0; i <= n; i++)
  {
    if(p[i] == 0)
      continue;

    os << p[i];
    if (i != 0)
      os << "x";

    if(i > 1)
      os << "^" << i;

    if(i < n)
      os << " + ";
  }
  return os;
}
