#include "Polynomial.h"
#include <vector>
#include <iostream>
using namespace std;

int main()
{
  //Hermite polynomial calculated by hand
  vector<double> coefficients = {0, 75, 0, 0.222222222, -0.031111111, -0.006444444, 0.002263889, -0.000913194, 0.000130527, -0.0000202236};
  Polynomial distance = coefficients[0]
                      + coefficients[1]*Polynomial({0, 1})
                      + coefficients[2]*Polynomial({0,0,1})
                      + coefficients[3]*Polynomial({0,0,1})*Polynomial({-3, 1})
                      + coefficients[4]*Polynomial({0,0,1})*Polynomial({-3, 1}).power(2)
                      + coefficients[5]*Polynomial({0,0,1})*Polynomial({-3, 1}).power(2)*Polynomial({-5, 1})
                      + coefficients[6]*Polynomial({0,0,1})*Polynomial({-3, 1}).power(2)*Polynomial({-5, 1}).power(2)
                      + coefficients[7]*Polynomial({0,0,1})*Polynomial({-3, 1}).power(2)*Polynomial({-5, 1}).power(2)*Polynomial({-8, 1})
                      + coefficients[8]*Polynomial({0,0,1})*Polynomial({-3, 1}).power(2)*Polynomial({-5, 1}).power(2)*Polynomial({-8, 1}).power(2)
                      + coefficients[9]*Polynomial({0,0,1})*Polynomial({-3, 1}).power(2)*Polynomial({-5, 1}).power(2)*Polynomial({-8, 1}).power(2)*Polynomial({-13, 1});

  Polynomial speed = distance.derivative();
  Polynomial b = -80.666666667 + speed; //55 mph = 80.666 ft/s

  //find root of b using Newton's method
  cout << "55 mph reached at t = " << b.root(5) << endl;
  Polynomial acceleration = speed.derivative();
  //speed has a critical value when acceleration=0
  double timeOfMaxSpeed = acceleration.root(12);
  cout << "Max speed of " << speed.evaluate(timeOfMaxSpeed) << " ft/s reached at t = " << timeOfMaxSpeed << endl;

  return 0;
}
