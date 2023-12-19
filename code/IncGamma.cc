#include <cmath>

#include <iostream>
#include <TCanvas.h>
#include <Math/SpecFuncMathCore.h>
#include "IncGamma.h"

using namespace ROOT;
using namespace Math;


IncGamma::IncGamma(){}

/*TF1* IncGamma::MyGamma()
{
TF1 *MyIncompleteGamma= new TF1("Inc","ROOT::Math::inc_gamma([0],x)");
    
        return MyIncompleteGamma;
};*/

double IncGamma::MyGamma(double a, double x)
{
double MyIncompleteGamma= ROOT::Math::inc_gamma(a,x);
    
        return MyIncompleteGamma;
};

TF1* IncGamma::Set_a(double a, TF1* gamma)
{
  gamma->SetParameter(0,a);
return gamma;
};

/*TF1* IncGamma::Set_lim(double lim)
{
  gamma->SetParameter(1,lim);
return gamma;
};*/
