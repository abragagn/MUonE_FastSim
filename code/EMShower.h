#ifndef EMSHOWER_h
#define EMSHOWER_h

#include <TH2.h>
#include <cmath>
#include "EMECALShowerParametrization.h"
#include "ECAL.h"
#include "IncGamma.h"


//#include "RandomEngineAndDistribution.h"
#include "GammaFunctionGenerator.h"
#include "RadialInterval.h"
#include <iostream>
#include "Math/Vector3D.h"
#include <TRandom3.h>
#include <TMath.h>
#include <TGraph.h>


#include <vector>

class GammaFunctionGenerator;

class EMShower
{
  typedef ROOT::Math::XYZVector XYZPoint;
  typedef pair<XYZPoint, double> Spot;
  typedef pair<unsigned int, double> Step;
  typedef vector<Step> Steps;
  typedef Steps::const_iterator step_iterator;

public:
  EMShower(GammaFunctionGenerator* gamma,
           EMECALShowerParametrization* const myParam,
           ECAL* const myGrid,
           bool bFixedLength,
           int nPart,
           double X0depth,
           vector<double> energy_in,
           vector<double> coo);

  virtual ~EMShower() { ; }

  /// Computes the steps before the real compute
  void prepareSteps();

  /// Compute the shower longitudinal and lateral development
  void compute();

  /// get the depth of the centre of gravity of the shower(s)
  //  inline double getMeanDepth() const {return globalMeanDepth;};

  /// get the depth of the maximum of the shower
  inline double getMaximumOfShower() const { return globalMaximum; }

  /// set the grid address
  void setGrid(ECAL* const myGrid) { theGrid = myGrid; }


private:
  // The longitudinal development ersatzt.
  double gam(double x, double a) const;

  // Energy deposited in the layer t-dt-> t, in units of E0 (initial energy)
  double deposit(double t, double a, double b, double dt);

  // Energy deposited between 0 and t, in units of E0 (initial energy)
  double deposit(double a, double b, double t);

  // Set the intervals for the radial development
  void setIntervals(unsigned icomp, RadialInterval& rad);

  // integer gamma function generator
  GammaFunctionGenerator* myGammaGenerator;

  // The parametrization
  EMECALShowerParametrization* const theParam;

  // The Calorimeter properties
  const ECALProperties* theECAL;

  // The crystal grid
  ECAL* theGrid;
  //TH2F* EcalGrid;
  bool bFixedLength_;
  int nPart;
  // The incident particle(s)
  //std::vector<const RawParticle*>* const thePart;

  double X0depth;

  std::vector<double> energy_in;
  std::vector<double> coo;

  // The basic quantities for the shower development.

  std::vector<double> theNumberOfSpots;
  std::vector<double> Etot;
  //std::vector<std::vector<double> > Etot_step;
  std::vector<double> Etot_step;

  std::vector<double> E;
  std::vector<double> photos;
  std::vector<double> T;
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> Ti;
  std::vector<double> TSpot;
  std::vector<double> aSpot;
  std::vector<double> bSpot;

  // F.B : Use the maximum of the shower rather the center of gravity
  //  std::vector<double> meanDepth;
  //  double globalMeanDepth;
  std::vector<double> maximumOfShower;
  std::vector<std::vector<double>> depositedEnergy;
  std::vector<double> meanDepth;
  double innerDepth, outerDepth;
  double Xi;
  double Yi;
  double globalMaximum;

  double totalEnergy;

  double numberPart;

  double numberSpot;

  // The steps for the longitudinal development
  Steps steps;
  unsigned nSteps;
  bool stepsCalculated;

  // Histos
  TH1F* Rad1;
  TH1F* Rad2;
  TH1F* Rad3;
  TH1F* Rad4;
  TH1F* RadTot;
  TH1F* en_1cell;
  TH1F* en_3x3cell;

  TH1F* Longit;

  //  Histos* myHistos;
  IncGamma myIncompleteGamma;

  // Random engine
  //const RandomEngineAndDistribution* random;

};

#endif