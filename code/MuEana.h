#ifndef MuEana_H
#define MuEana_H

///////////////////////////////////////////////
// Classes defining MuE analysis variables
//
// G.Abbiendi  4/Dec/2018
///////////////////////////////////////////////

#include "Math/Vector4D.h"

namespace MuE
{

typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;

// Analysis variables in the output tree

class MuonOut
{
public:
  Double_t Emu; // muon energy
  Double_t thmu; // muon theta (in mrad)
  Double_t phmu; // muon phi (from -pi to +pi)
  Double_t t13; // Mandelstam t (muon leg)
  Double_t x13; // Feynman x (muon leg)
  Double_t cooxmu;
  Double_t cooymu;
  Double_t tar;
  Double_t n_phot;
  Double_t def_angle_mu;
  Double_t def_angle_mu_smear;

  MuonOut() :
    Emu(0), thmu(0), phmu(0), t13(0), x13(0), cooxmu(0), cooymu(0), tar(0),
    def_angle_mu(0), def_angle_mu_smear(0) {};

  virtual ~MuonOut() {};
  ClassDef(MuonOut, 1)
};

class Electron
{
public:
  Double_t Ee; // electron energy
  Double_t the; // electron theta (in mrad)
  Double_t phe; // electron phi (from -pi to +pi)
  Double_t t24; // Mandelstam t (electron leg)
  Double_t x24; // Feynman x (electron leg)
  Double_t tt_e; // t computed from electron angle with LO formulas
  Double_t xt_e; // x computed from electron angle with LO formulas
  Double_t cooxe;
  Double_t cooye;
  Double_t def_angle_e;
  Double_t def_angle_e_smear;

  Double_t Ecell1;
  Double_t Ecell2;
  Double_t Ecell3;
  Double_t Ecell4;
  Double_t Ecell5;
  Double_t Ecell6;
  Double_t Ecell7;
  Double_t Ecell8;
  Double_t Ecell9;
  Double_t Ecell10;
  Double_t Ecell11;
  Double_t Ecell12;
  Double_t Ecell13;
  Double_t Ecell14;
  Double_t Ecell15;
  Double_t Ecell16;
  Double_t Ecell17;
  Double_t Ecell18;
  Double_t Ecell19;
  Double_t Ecell20;
  Double_t Ecell21;
  Double_t Ecell22;
  Double_t Ecell23;
  Double_t Ecell24;
  Double_t Ecell25;
  Double_t n_cell_e;

  Electron() :
    Ee(0), the(0), phe(0), t24(0), x24(0), tt_e(0), xt_e(0), cooxe(0),
    cooye(0), def_angle_e(0), def_angle_e_smear(0), Ecell1(0), Ecell2(0),
    Ecell3(0), Ecell4(0), Ecell5(0),  Ecell6(0), Ecell7(0), Ecell8(0),
    Ecell9(0), Ecell10(0), Ecell11(0), Ecell12(0), Ecell13(0), Ecell14(0),
    Ecell15(0), Ecell16(0), Ecell17(0), Ecell18(0), Ecell19(0), Ecell20(0),
    Ecell21(0), Ecell22(0), Ecell23(0), Ecell24(0), Ecell25(0), n_cell_e(0) {};

  virtual ~Electron() {};
  ClassDef(Electron, 1)
};

class MuEpair
{
public:
  Double_t deltaPhi; // acoplanarity (deltaPhi)
  Double_t openingAngle; // opening angle mu-e out in the Lab
  Double_t tripleProduct; // triple product btw normalized vectors i . mu x e
  MuEpair() : deltaPhi(0), openingAngle(0), tripleProduct(0) {};
  virtual ~MuEpair() {};
  ClassDef(MuEpair, 1)
};

class KineVars
{
public:
  Double_t t13; // Mandelstam t (muon leg)
  Double_t t24; // Mandelstam t (electron leg)
  Double_t x13; // Feynman x (muon leg)
  Double_t x24; // Feynman x (electron leg)
  Double_t tt_e; // t computed from electron angle with LO formulas
  Double_t xt_e; // x computed from electron angle with LO formulas

  Double_t Ee; // electron energy
  Double_t Emu; // muon energy
  Double_t the; // electron theta (in mrad)
  Double_t thmu; // muon theta (in mrad)
  Double_t phe; // electron phi (from -pi to +pi)
  Double_t phmu; // muon phi (from -pi to +pi)

  Double_t cooxmu;
  Double_t cooymu;
  Double_t cooxe;
  Double_t cooye;

  Double_t tar;
  Double_t n_phot;

  Double_t def_angle_mu;
  Double_t def_angle_e;

  Double_t def_angle_mu_smear;
  Double_t def_angle_e_smear;

  Double_t Ecell1;
  Double_t Ecell2;
  Double_t Ecell3;
  Double_t Ecell4;
  Double_t Ecell5;
  Double_t Ecell6;
  Double_t Ecell7;
  Double_t Ecell8;
  Double_t Ecell9;
  Double_t Ecell10;
  Double_t Ecell11;
  Double_t Ecell12;
  Double_t Ecell13;
  Double_t Ecell14;
  Double_t Ecell15;
  Double_t Ecell16;
  Double_t Ecell17;
  Double_t Ecell18;
  Double_t Ecell19;
  Double_t Ecell20;
  Double_t Ecell21;
  Double_t Ecell22;
  Double_t Ecell23;
  Double_t Ecell24;
  Double_t Ecell25;
  Double_t n_cell_e;

  Double_t deltaPhi; // acoplanarity (deltaPhi)
  Double_t openingAngle; // opening angle mu-e out in the Lab
  Double_t tripleProduct; // triple product btw normalized vectors i . mu x e

  void SetMuonOut(const MuonOut & m)
  {
    Emu = m.Emu;
    thmu = m.thmu;
    phmu = m.phmu;
    t13 = m.t13;
    x13 = m.x13;
    cooxmu = m.cooxmu;
    cooymu = m.cooymu;
    tar = m.tar;
    n_phot = m.n_phot;
    def_angle_mu = m.def_angle_mu;
    def_angle_mu_smear = m.def_angle_mu_smear;
  };

  void SetElectron(const Electron & e)
  {
    Ee = e.Ee;
    the = e.the;
    phe = e.phe;
    t24 = e.t24;
    x24 = e.x24;
    tt_e = e.tt_e;
    xt_e = e.xt_e;
    cooxe = e.cooxe;
    cooye = e.cooye;
    def_angle_e = e.def_angle_e;
    def_angle_e_smear = e.def_angle_e_smear;
    Ecell1 = e.Ecell1;
    Ecell2 = e.Ecell2;
    Ecell3 = e.Ecell3;
    Ecell4 = e.Ecell4;
    Ecell5 = e.Ecell5;
    Ecell6 = e.Ecell6;
    Ecell7 = e.Ecell7;
    Ecell8 = e.Ecell8;
    Ecell9 = e.Ecell9;
    Ecell10 = e.Ecell10;
    Ecell11 = e.Ecell11;
    Ecell12 = e.Ecell12;
    Ecell13 = e.Ecell13;
    Ecell14 = e.Ecell14;
    Ecell15 = e.Ecell15;
    Ecell16 = e.Ecell16;
    Ecell17 = e.Ecell17;
    Ecell18 = e.Ecell18;
    Ecell19 = e.Ecell19;
    Ecell20 = e.Ecell20;
    Ecell21 = e.Ecell21;
    Ecell22 = e.Ecell22;
    Ecell23 = e.Ecell23;
    Ecell24 = e.Ecell24;
    Ecell25 = e.Ecell25;
    n_cell_e = e.n_cell_e;
  };

  void SetMuEpair(const MuEpair & d)
  {
    deltaPhi = d.deltaPhi;
    openingAngle = d.openingAngle;
    tripleProduct = d.tripleProduct;
  };

  KineVars() :
    t13(0), t24(0), x13(0), x24(0), tt_e(0), xt_e(0), Ee(0), Emu(0), the(0),
    thmu(0), phe(0), phmu(0), tar(0), deltaPhi(0), openingAngle(0),
    tripleProduct(0) {};

  KineVars(const MuonOut & m, const Electron & e, const MuEpair & d) :
    t13(m.t13), t24(e.t24), x13(m.x13), x24(e.x24), tt_e(e.tt_e),
    xt_e(e.xt_e), Ee(e.Ee), Emu(m.Emu), the(e.the), thmu(m.thmu),
    phe(e.phe), phmu(m.phmu), cooxmu(m.cooxmu), cooymu(m.cooymu),
    cooxe(e.cooxe), cooye(e.cooye), tar(m.tar), n_phot(m.n_phot),
    def_angle_mu(m.def_angle_mu), def_angle_e(e.def_angle_e),
    def_angle_mu_smear(m.def_angle_mu_smear),
    def_angle_e_smear(e.def_angle_e_smear), Ecell1(e.Ecell1),
    Ecell2(e.Ecell2), Ecell3(e.Ecell3), Ecell4(e.Ecell4), Ecell5(e.Ecell5),
    Ecell6(e.Ecell6), Ecell7(e.Ecell7), Ecell8(e.Ecell8), Ecell9(e.Ecell9),
    Ecell10(e.Ecell10), Ecell11(e.Ecell11), Ecell12(e.Ecell1),
    Ecell13(e.Ecell13), Ecell14(e.Ecell14), Ecell15(e.Ecell15),
    Ecell16(e.Ecell16), Ecell17(e.Ecell17), Ecell18(e.Ecell18),
    Ecell19(e.Ecell19), Ecell20(e.Ecell20), Ecell21(e.Ecell21),
    Ecell22(e.Ecell22), Ecell23(e.Ecell23), Ecell24(e.Ecell24),
    Ecell25(e.Ecell25), n_cell_e(e.n_cell_e), deltaPhi(d.deltaPhi),
    openingAngle(d.openingAngle), tripleProduct(d.tripleProduct) {};

  virtual ~KineVars() {};
  ClassDef(KineVars, 2)
};

class Photon
{
public:
  Double_t energy;    // photon energy in the Lab frame
  Double_t theta;     //   "    theta in the Lab frame (in mrad)
  Double_t phi;       //   "    phi in the Lab frame (in rad)
  Double_t energyCoM; // photon energy in the Centre-of-Mass frame

  Double_t def_angle_ph; // deflection angle from the incoming muon direction

  Double_t cooxph;    // x coordinate photon
  Double_t cooyph;    // y coordinate photon
  Double_t n_cell_ph;

  Photon() :
    energy(-1), theta(-1), phi(0), energyCoM(-1), def_angle_ph(-1),
    cooxph(-1), cooyph(-1), n_cell_ph(0) {};

  virtual ~Photon() {};
  ClassDef(Photon, 1)
};

class MuEana
{

public:
  Int_t RunNr;
  Long_t EventNr;
  Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO, wgt_NLO;   // event weights
  PxPyPzEVector p_mu_in;
  PxPyPzEVector p_e_in;
  PxPyPzEVector p_mu_out;
  PxPyPzEVector p_e_out;
  std::vector<PxPyPzEVector> V_photons;
  KineVars genKin;    // kinematic variables at Generator-level for e and mu tracks
  KineVars detKin;    // kinematic variables at Detector-level for e and mu tracks
  //KineVars detKinFin; // kinematic variables at Detector-level after MCS, divergence and rotation for e and mu tracks
  Photon photon1;     // photons kinematic variables at Det level
  Photon photon2;     // initialized variables for no photon event
  //std::vector<Photon> load_photons;    // vector containing the two photons

  MuEana() :
    RunNr(0), EventNr(0), wgt_full(0), wgt_norun(0), wgt_lep(0), wgt_LO(0),
    wgt_NLO(0) {};

  virtual ~MuEana() {};
  ClassDef(MuEana, 3)
};
}

#endif
