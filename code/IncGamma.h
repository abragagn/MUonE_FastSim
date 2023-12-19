#ifndef INCGAMMA_h
#define INCGAMMA_h

#include <TF1.h>

#include <cmath>
#include <iostream>
using namespace std;

class IncGamma 
{
public:

//costruttore
IncGamma();
//distruttore
~ IncGamma(){}

double MyGamma(double a, double x);
TF1* Set_a(double a, TF1* gamma);
//TF1* Set_lim(double lim, TF1* gamma);

};
#endif