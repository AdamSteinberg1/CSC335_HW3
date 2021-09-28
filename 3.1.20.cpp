#include "Polynomial.h"
#include <vector>
#include <iostream>
using namespace std;

int main()
{
  vector<pair<long double, long double>> sample1 = {
    {0, 6.67},
    {6, 17.33},
    {10, 42.67},
    {13, 37.33},
    {17, 30.10},
    {20, 29.31},
    {28, 28.74},
  };

  vector<pair<long double, long double>> sample2 = {
    {0, 6.67},
    {6, 16.11},
    {10, 18.89},
    {13, 15.00},
    {17, 10.56},
    {20, 9.44},
    {28, 8.89},
  };

  Polynomial p1 = lagrange(sample1);
  Polynomial p2 = lagrange(sample2);

  cout << "p1 =  "<<p1<<endl;
  cout << "p2 =  "<<p2<<endl;
  return 0;
}
