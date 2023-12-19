#ifndef ECAL_h
#define ECAL_h

#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>


#include <cmath>
#include <map>
#include <TMatrixD.h>
#include <TGraph.h>
#include "MuEana.h"

#include <iostream>

using namespace std;

class ECAL 
//: public TH2
{
public:

//costruttore
ECAL(double nbinsx, 
    double xlow, 
    double xup, 
    double nbinsy, 
    double ylow, 
    double yup);
//distruttore
~ ECAL(){}


//double radlen;

void CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup);
TH2F* GiveEcalGrid();
double GiveCentralCell(double coox,double cooy);
void SetEnergy(double energy);
int* GiveArray3x3(int n);
double AddHitCoo(double r,double phi,double xi,double yi,double w);
void AddHitCooDepth(double r, double phi,double xi, double yi, double w, double depth, double deX0depthoffset_pth);
double* Draw_ECAL(int i);
double* EnergyContent();
void Print_();
inline void setSpotEnergy(double e) { spotEnergy = e; }

private:
const double nbinX;
const double nbinY;
const double Xlow;
const double Xup;
const double Ylow;
const double Yup;
double spotEnergy;
double energy_IN;
typedef map<int, int>  n_cell;
n_cell number;
n_cell Rev_number;

TH1F* Energy_dist;
TH1F* Energy_dist1;
TH1F* Energy_dist3x3;
TH2F* EcalGrid;

int *Array9;
};
#endif