#ifndef EMECALShowerParametrization_H
#define EMECALShowerParametrization_H

#include "ECALProprieties.h"

/** 
 * Electromagnetic Shower parametrization utilities according to 
 * G. Grindhammer and S. Peters, hep-ex/0001020, Appendix A
 *
 * \author Patrick Janot
 * \date: 25-Jan-2004
 */
#include <vector>
#include <cmath>

class EMECALShowerParametrization {
public:
  EMECALShowerParametrization(const ECALProperties* ecal,
                              const std::vector<double> coreIntervals,
                              const std::vector<double> tailIntervals,
                              double RCFact, //= 1.,
                              double RTFact)// = 1.)
      : theECAL(ecal),
        theCore(coreIntervals),
        theTail(tailIntervals),
        theRcfactor(RCFact),
        theRtfactor(RTFact) {}

  virtual ~EMECALShowerParametrization() {}

  //====== Longitudinal profiles =======

  // -------- Average -------

  inline double meanT(double lny) const {
      return meanTHom(lny);
  }

  inline double meanAlpha(double lny) const {
      return meanAlphaHom(lny);
  }

  // Average Homogeneous

  inline double meanTHom(double lny) const { return lny - 0.858; }
  //return lny-1.23; }

  inline double meanAlphaHom(double lny) const {return 0.21 + (0.492 + 2.38 /theECAL->theZeff()) * lny; } 
    /* / theECAL->theZeff()theECAL->theZeff() */
    
  // ---- Fluctuated longitudinal profiles ----

  inline double meanLnT(double lny) const {
      return meanLnTHom(lny);
}

  inline double sigmaLnT(double lny) const {
      return sigmaLnTHom(lny);
  }

  inline double meanLnAlpha(double lny) const {
      return meanLnAlphaHom(lny);
  }

  inline double sigmaLnAlpha(double lny) const {
      return sigmaLnAlphaHom(lny);
  }

  inline double correlationAlphaT(double lny) const {
      return correlationAlphaTHom(lny);
  }

  // Fluctuated longitudinal profiles homogeneous

  inline double meanLnTHom(double lny) const { return std::log(lny - 0.812); }

  inline double sigmaLnTHom(double lny) const { return 1. / (-1.4 + 1.26 * lny); }

  inline double meanLnAlphaHom(double lny) const { return std::log(0.81 + (0.458 + 2.26 / theECAL->theZeff()) * lny); }
    /* / theECAL->theZeff()theECAL->theZeff() */

  inline double sigmaLnAlphaHom(double lny) const { return 1. / (-0.58 + 0.86 * lny); }

  inline double correlationAlphaTHom(double lny) const { return 0.705 - 0.023 * lny; }

  //====== Radial profiles =======

  // ---- Radial Profiles ----

  inline double rC(double tau, double E) const {
      return rCHom(tau, E);
  }

  inline double rT(double tau, double E) const {
      return rTHom(tau, E);
  }

  inline double p(double tau, double E) const {
      return pHom(tau, E);
  }

  // Radial Profiles

  inline double rCHom(double tau, double E) const { return theRcfactor * (z1(E) + z2() * tau); }

  inline double rTHom(double tau, double E) const {
    return theRtfactor * k1() * (std::exp(k3() * (tau - k2())) + std::exp(k4(E) * (tau - k2())));
  }

  inline double pHom(double tau, double E) const {
    double arg = (p2() - tau) / p3(E);
    return p1() * std::exp(arg - std::exp(arg));
  }

  // ---- Fluctuations of the radial profiles ----

  inline double nSpots(double E) const {
      return nSpotsHom(E);
  }

  inline double meanTSpot(double T) const {
      return meanTSpotHom(T);
  }

  inline double meanAlphaSpot(double alpha) const {
      return meanAlphaSpotHom(alpha);
  }

  // Fluctuations of the radial profiles

  inline double nSpotsHom(double E) const { return 93. * std::log(theECAL->theZeff()) * std::pow(E, 0.876); }  /* / theECAL->theZeff()theECAL->theZeff() */

  inline double meanTSpotHom(double T) const { return T * (0.698 + 0.00212 * theECAL->theZeff()); }
     /* / theECAL->theZeff()theECAL->theZeff() */

  inline double meanAlphaSpotHom(double alpha) const { return alpha * (0.639 + 0.00334 * theECAL->theZeff()); }
/* / theECAL->theZeff()theECAL->theZeff() */
    
  inline const ECALProperties* ecalProperties() const { return theECAL; }

  inline const std::vector<double>& getCoreIntervals() const { return theCore; }
  inline const std::vector<double>& getTailIntervals() const { return theTail; }

private:
  const ECALProperties* theECAL;

  const std::vector<double> theCore;
  const std::vector<double> theTail;

  double theRcfactor;
  double theRtfactor;

  double p1() const { return 2.632 - 0.00094 * theECAL->theZeff(); }
    /* / theECAL->theZeff()theECAL->theZeff() */
  double p2() const { return 0.401 + 0.00187 * theECAL->theZeff(); }
    /* / theECAL->theZeff()theECAL->theZeff() */
  double p3(double E) const { return 1.313 - 0.0686 * std::log(E); }

  double z1(double E) const { return 0.0251 + 0.00319 * std::log(E); }
  double z2() const { return 0.1162 - 0.000381 * theECAL->theZeff() ; }
     /* / theECAL->theZeff()theECAL->theZeff() */

  double k1() const { return 0.6590 - 0.00309 * theECAL->theZeff() ; }
     /* / theECAL->theZeff()theECAL->theZeff() */
  double k2() const { return 0.6450; }
  double k3() const { return -2.59; }
  double k4(double E) const { return 0.3585 + 0.0421 * std::log(E); }
};

#endif