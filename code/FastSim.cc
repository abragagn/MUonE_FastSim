#include <cmath>
#include <iostream>
#include "TRandom.h"
#include "TMath.h"
#include "ElasticState.h"
#include "ResolutionModels.h"
#include "FastSim.h"
#include "TMatrixF.h"
#include "TF2.h"
#include "TMatrixD.h"
#include "TMatrixFBase.h"
#include <TMatrixFSym.h>
#include "TString.h"
#include "TVector3.h"
#include "EMShower.h"



#include <TApplication.h>

using namespace MuE;
using namespace std;
using namespace ROOT::Math;

// from PDG Book 2018
const Double_t FastSim::mm_PDG = 105.6583745 * 0.001;
const Double_t FastSim::me_PDG = 0.5109989461 * 0.001;

// detector geometry
const Double_t FastSim::station_Length = 1.; // 1m
const Double_t FastSim::detector_Size = 0.1; // 10cm

FastSim::FastSim(const MuE::MCpara & pargen, const MuE::FS_Input & parsim,  bool _debug_) :
  mm(pargen.mass_mu), me(pargen.mass_e), Qbeam(pargen.charge_mu),
  Ebeam(pargen.Ebeam), EbeamRMS(pargen.EbeamRMS), model(parsim.model),
  MSopt(parsim.MSopt), twosteps(parsim.twosteps), thickness(parsim.thickness),
  intrinsic_resolution(parsim.resolution), zBias(parsim.zBias),
  zSigma(parsim.zSigma), debug(_debug_), Minv(0)
{
  if (std::abs(mm - mm_PDG) / mm_PDG > 1e-7)
  {
    cout << "\n" << "***WARNING: muon mass = " << mm
         << " is different from the PDG mass: " << mm_PDG << endl;
  }
  if (std::abs(me - me_PDG) / me_PDG > 1e-7)
  {
    cout << "\n" << "***WARNING: electron mass = " << me
         << " is different from the PDG mass: " << me_PDG << endl;
  }

  if (model == 0)
  {
    cout << "\n Simple detector model, DetModel = " << model
         << ", option = " << MSopt << ", twosteps = " << twosteps << endl;
    cout << " material thickness = " << thickness << " X0,"
         << " intrinsic angular resolution = " << intrinsic_resolution << " mrad " << endl;
  }
  else if (model == 1)
  {
    cout << "\n Antonios detector model, 3 parameters" << endl;
    twosteps = false;
  }
  else
  {
    cout << "\n" << "*** ERROR : FastSim, undefined detector model with DetModel = " << model << endl;
    std::exit(999);
  }

  thetaMax_ideal = detector_Size / 2 / station_Length;
  thetaMax_cor = thetaMax_ideal * (1. - zBias * 1e-6 / station_Length);

  if (zSigma > 1e-7)
  {
    zSigma_switch = true;
    cout << " ===> setting the zSigma_switch to true" << endl;
  }
  else
  {
    zSigma_switch = false;
    cout << " ===> setting the zSigma_switch to false" << endl;
  }

  cout << "\n Station Z length : bias = " << zBias << " (um),"
       << " sigma = " << zSigma << " (um) " << endl;
  cout << " thetaMax_ideal = " << thetaMax_ideal * 1e3 << " mrad,"
       << " thetaMax_cor = " << thetaMax_cor * 1e3 << " mrad" << endl;
}


// process an event
//
void FastSim::Process(const MuE::Event & event, GammaFunctionGenerator* & gamma,
                      EMECALShowerParametrization* const & myParam, ECAL* const & myGrid)
{

  if (debug) cout << "\n Process:  Run = " << event.RunNr << " , Event = " << event.EventNr << endl;

  //initial particle 4-momenta, only the muon is moving because the other don't exist yet
  p_mu_in.SetPxPyPzE(event.muin.px, event.muin.py, event.muin.pz, event.muin.E());
  p_e_in.SetPxPyPzE(0, 0, 0, me);

  p_system = p_mu_in + p_e_in;
  Double_t s = mm * mm + me * me + 2 * me * p_mu_in.E();
  Minv = sqrt(s);
  int i = event.EventNr;

  //  Double_t pcm = P_2bodies_CoM(Minv, mm, me);
  if (debug) cout << "\n" << "Incoming muon energy = " << std::fixed << setprecision(3) << p_mu_in.E() << " GeV" << endl;
  //      << " GeV, s = "<<s<<" GeV^2, sqrt(s) = "<<Minv<<" GeV, pcm = "<<pcm<<" GeV"<<endl;

  bool found_mu = false;
  bool found_e = false;
  UInt_t nphot = 0;
  photons.clear();

  // matrix of particle coordinates il the laboratory frame
  TMatrixD coordinates(4, 2);
  coordinates[0][0] = 0; // muon x coordinate
  coordinates[0][1] = 0; // muon y coordinate
  coordinates[1][0] = 0; // e x coordinate
  coordinates[1][1] = 0; // e y coordinate
  coordinates[2][0] = 0; // ph1 x coordinate
  coordinates[2][1] = 0; // ph1 y coordinate
  coordinates[3][0] = 0; // ph2 x coordinate
  coordinates[3][1] = 0; // ph2 y coordinate

  TMatrixD def_angle(2, 1);
  def_angle[0][0] = 0;
  def_angle[1][0] = 0;

  TMatrixD def_angle_smear(2, 1);
  def_angle_smear[0][0] = 0;
  def_angle_smear[1][0] = 0;


  myGrid->CreateGrid(5, -7.125, 7.125, 5, -7.125, 7.125);

  // particle 4-momentum after the scattering, at generation level
  for (const auto &p : event.fspart)
  {
    if (!found_mu && ((p.pdgId) == -Qbeam * 13))
    {
      found_mu = true;
      p_mu_out.SetPxPyPzE(p.px, p.py, p.pz, p.E());
    }
    else if (!found_e && (p.pdgId == 11))
    {
      found_e = true;
      p_e_out.SetPxPyPzE(p.px, p.py, p.pz, p.E());
    }
    else if (nphot < 2 && p.pdgId == 22)
    {
      nphot++;
      photons.push_back(PxPyPzEVector(p.px, p.py, p.pz, p.E()));
    }
    else
    {
      cout << "\n*** ERROR: unexpected particle in Event " << event.EventNr << endl;
      cout << "\t incoming muon: [id:" << event.muin.pdgId << "] (" << event.muin.px << ", " << event.muin.py << ", " << event.muin.pz << ")" << endl;
      cout << "\t final-state particles: " << endl;
      for (const auto &fsp : event.fspart)
        cout << "\t [id:" << fsp.pdgId << "] (" << fsp.px << ", " << fsp.py << ", " << fsp.pz << ")" << endl;
      std::exit(222);
    }
  }

  tar = gRandom->Integer(2);
  vertex = gRandom->Uniform(); // where inc. muon interacts in the Beryllium tar.

  //make a copy of the initial and after the scattering 4-momenta in order to modify them later without losing the generator level ones
  PxPyPzEVector p_mu_in_mod = p_mu_in;
  PxPyPzEVector p_mu_out_mod = p_mu_out;
  PxPyPzEVector p_e_in_mod = p_e_in;
  PxPyPzEVector p_e_out_mod = p_e_out;
  std::vector<PxPyPzEVector> photons_mod = photons;
  n_phot = photons.size();

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //TEST PER IL CALORIMETRO CON FASCIO DI ELETTRONI DA 500MEV/1GEV ECC//////////
  //p_e_out_mod.SetE(120);
  //FINE TEST RISOLUZIONE CALORIMETRO///////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  // Generator-level kinematics
  genKin = LoadKineVars(p_mu_in, p_e_in, p_mu_out, p_e_out, coordinates, def_angle, def_angle_smear, tar, n_phot, myGrid);


  // detector level cinematic with divergence, MCS and propagation
  Propagate(p_mu_in_mod, p_mu_out_mod, p_e_out_mod, photons_mod, coordinates, def_angle, def_angle_smear);
  detKin = LoadKineVars(p_mu_in_mod, p_e_in_mod, p_mu_out_mod, p_e_out_mod, coordinates, def_angle, def_angle_smear, tar, nphot, myGrid);


  LoadPhotons(photons_mod, coordinates, nphot, p_mu_in_mod, myGrid); // only rotation and propagation


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //TEST PER IL CALORIMETRO SCELGP POSIZIONE DELL'ELETTRONE
  //NELL'ANGOLO O AL CENTRO/////////////////////////////////////////////////////
  //coordinates[1][0]=0;
  //coordinates[1][1]=0;
  //FINE TEST RISOLUZIONE CALORIMETRO///////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////


  //creation of the shower and of the ECAL
  Double_t r_mu = sqrt((coordinates[0][0] * coordinates[0][0])
                        + (coordinates[0][1] * coordinates[0][1])) * 100; //rmu in cm

  PrepareEMShower(detKin, p_e_out_mod, photons_mod, coordinates, i, r_mu, gamma,
                  myParam, myGrid);

  photons_mod.clear();

}





KineVars FastSim::LoadKineVars(const PxPyPzEVector & p_mu_in,
                               const PxPyPzEVector & p_e_in,
                               const PxPyPzEVector & p_mu_out,
                               const PxPyPzEVector & p_e_out,
                               TMatrixD & coo,
                               TMatrixD & def_angle,
                               TMatrixD & def_angle_smear,
                               const Double_t & tar,
                               const Double_t & n_phot,
                               ECAL* const & myGrid)
{

  MuonOut m = p_mu_out.E() > 0 ?
    DetMuonOut(p_mu_in, p_mu_out, coo, def_angle, def_angle_smear, tar, n_phot) : MuonOut();

  Electron e = p_e_out.E() > 0 ?
    DetElectron(p_e_in, p_e_out, coo, def_angle, def_angle_smear, myGrid) : Electron();

  MuEpair d = (p_mu_out.E() > 0 && p_e_out.E() > 0) ?
    DetMuEpair(p_mu_in, p_e_in, p_mu_out, p_e_out) : MuEpair();

  return KineVars(m, e, d);
}


MuonOut FastSim::DetMuonOut(const PxPyPzEVector & p_mu_in,
                            const PxPyPzEVector & p_mu_out,
                            TMatrixD & coo,
                            TMatrixD & def_angle,
                            TMatrixD & def_angle_smear,
                            const Double_t & tar,
                            const Double_t & n_phot) const
{
  MuonOut m;
  m.Emu = p_mu_out.E();
  m.thmu = 1e3 * p_mu_out.Theta();
  m.phmu = p_mu_out.Phi();

  PxPyPzEVector q13 = p_mu_in - p_mu_out;
  m.t13 = q13.M2();
  ElasticState emu(Ebeam, mm, me);
  m.x13 = emu.X(m.t13);
  m.cooxmu = coo[0][0];
  m.cooymu = coo[0][1];
  m.def_angle_mu = def_angle[0][0];
  m.def_angle_mu_smear = def_angle_smear[0][0];
  m.tar = tar;
  m.n_phot = n_phot;


  return m;
}


Electron FastSim::DetElectron(const PxPyPzEVector & p_e_in,
                              const PxPyPzEVector & p_e_out,
                              TMatrixD & coo,
                              TMatrixD & def_angle,
                              TMatrixD & def_angle_smear,
                              ECAL* const & myGrid) const
{
  Electron e;
  e.Ee = p_e_out.E();
  e.the = 1e3 * p_e_out.Theta();
  //e.the = 1e3* def_angle[1][0];
  e.phe = p_e_out.Phi();

  // Note: here Ebeam is the average beam energy, so tt_e and xt_e are defined under this assumption
  MuE::ElasticState emu(Ebeam, mm, me, e.the);
  e.tt_e = emu.GetT();
  e.xt_e = emu.GetX();

  PxPyPzEVector q24 = p_e_in - p_e_out;
  // instead these are exact
  e.t24 = q24.M2();
  e.x24 = emu.X(e.t24);
  e.cooxe = coo[1][0];
  e.cooye = coo[1][1];
  e.def_angle_e = def_angle[1][0];
  e.def_angle_e_smear = def_angle_smear[1][0];
  e.n_cell_e = myGrid->GiveCentralCell(coo[1][0] * 100, coo[1][1] * 100); //in cm


  /*TH2F* Ecal=myGrid->GiveEcalGrid();
  if (Ecal->GetEntries()!=0)
  {
    double *Ecell=myGrid->EnergyContent();

    e.Ecell1=Ecell[0];
    e.Ecell2=Ecell[1];
    e.Ecell3=Ecell[2];
    e.Ecell4=Ecell[3];
    e.Ecell5=Ecell[4];
    e.Ecell6=Ecell[5];
    e.Ecell7=Ecell[6];
    e.Ecell8=Ecell[7];
    e.Ecell9=Ecell[8];
    e.Ecell10=Ecell[9];
    e.Ecell11=Ecell[10];
    e.Ecell12=Ecell[11];
    e.Ecell13 = Ecell[12];
    e.Ecell14 = Ecell[13];
    e.Ecell15 = Ecell[14];
    e.Ecell16 = Ecell[15];
    e.Ecell17 = Ecell[16];
    e.Ecell18 = Ecell[17];
    e.Ecell19 = Ecell[18];
    e.Ecell20 = Ecell[19];
    e.Ecell21 = Ecell[20];
    e.Ecell22 = Ecell[21];
    e.Ecell23 = Ecell[22];
    e.Ecell24 = Ecell[23];
    e.Ecell25 = Ecell[24];

    delete Ecal;
  }
  else
  {
    e.Ecell1 = 0;
    e.Ecell2 = 0;
    e.Ecell3 = 0;
    e.Ecell4 = 0;
    e.Ecell5 = 0;
    e.Ecell6 = 0;
    e.Ecell7 = 0;
    e.Ecell8 = 0;
    e.Ecell9 = 0;
    e.Ecell10 = 0;
    e.Ecell11 = 0;
    e.Ecell12 = 0;
    e.Ecell13 = 0;
    e.Ecell14 = 0;
    e.Ecell15 = 0;
    e.Ecell16 = 0;
    e.Ecell17 = 0;
    e.Ecell18 = 0;
    e.Ecell19 = 0;
    e.Ecell20 = 0;
    e.Ecell21 = 0;
    e.Ecell22 = 0;
    e.Ecell23 = 0;
    e.Ecell24 = 0;
    e.Ecell25 = 0;

    delete Ecal;
  }*/
  return e;
}


MuEpair FastSim::DetMuEpair(const PxPyPzEVector & p_mu_in,
                            const PxPyPzEVector & p_e_in,
                            const PxPyPzEVector & p_mu_out,
                            const PxPyPzEVector & p_e_out) const
{
  MuEpair d;

  // acoplanarity = deltaPhi as long as the incoming muon is collinear with z-axis
  Double_t deltaPhi = p_e_out.Phi() - p_mu_out.Phi();
  if (deltaPhi < 0.) deltaPhi = deltaPhi + 2.*TMath::Pi();
  deltaPhi = deltaPhi - TMath::Pi();
  d.deltaPhi = deltaPhi;

  XYZVector nvec_mu_in = p_mu_in.Vect().Unit();
  XYZVector nvec_mu_out = p_mu_out.Vect().Unit();
  XYZVector nvec_e_out = p_e_out.Vect().Unit();

  Double_t dotProduct = nvec_mu_out.Dot(nvec_e_out);
  d.openingAngle = std::abs(dotProduct) < 1. ? 1000.*acos(dotProduct) : 0.;

  XYZVector crossProduct = nvec_mu_out.Cross(nvec_e_out);
  d.tripleProduct = nvec_mu_in.Dot(crossProduct);

  return d;
}

void FastSim::LoadPhotons(const std::vector<PxPyPzEVector> & photons,
                          TMatrixD & coo,
                          Double_t n_photons1,
                          const PxPyPzEVector & p_mu_in,
                          ECAL* const & myGrid)
{ //, std::vector<Photon> load_photons) {
  //SHOULD WORK FOR TWO NOW

  auto n_photons = photons.size();

  //MuE::Photon photon_1;
  //MuE::Photon photon_2;
  //MuE::Photon no_photon;

  //photon1=no_photon;
  //photon2=no_photon;

  // we fill them as if there were no ohotons, so with the default values
  // load_photons.at(0)=no_photon;
  // load_photons.at(1)=no_photon;

  // if there are photons we fill it with the right number of photons
  if (n_photons > 0)
  {
    for (unsigned int i = 0; i < n_photons; i++)
    {
      if (i == 0)
      {
        PxPyPzEVector p_gamma_Lab = photons[0];
        PxPyPzEVector p_gamma_CoM = Lorentz_ToCoM(p_gamma_Lab);

        photon1.def_angle_ph = Def_angle_ph(p_mu_in, p_gamma_Lab);

        photon1.energy = p_gamma_Lab.E();
        photon1.theta = p_gamma_Lab.Theta() * 1e3;
        photon1.phi = p_gamma_Lab.Phi();
        photon1.energyCoM = p_gamma_CoM.E();
        photon1.cooxph = coo[2][0];
        photon1.cooyph = coo[2][1];
        photon1.n_cell_ph = myGrid->GiveCentralCell(coo[3][0] * 100, coo[3][1] * 100); //in cm
      }
      else if (i == 1)
      {
        PxPyPzEVector p_gamma_Lab = photons[1];
        PxPyPzEVector p_gamma_CoM = Lorentz_ToCoM(p_gamma_Lab);

        photon2.def_angle_ph = Def_angle_ph(p_mu_in, p_gamma_Lab);

        photon2.energy = p_gamma_Lab.E();
        photon2.theta = p_gamma_Lab.Theta() * 1e3;
        photon2.phi = p_gamma_Lab.Phi();
        photon2.energyCoM = p_gamma_CoM.E();
        photon2.cooxph = coo[3][0];
        photon2.cooyph = coo[3][1];
        photon2.n_cell_ph = myGrid->GiveCentralCell(coo[3][0] * 100, coo[3][1] * 100); //in cm
      }
    }
  }
}

// propagation of the particles inside the tracker, in necessary here it is possible to add other stations
//
void FastSim::Propagate(PxPyPzEVector & mu_in_mod, 
                        PxPyPzEVector & mu_out_mod,
                        PxPyPzEVector & e_out_mod,
                        std::vector<PxPyPzEVector> & phot_mod,
                        TMatrixD & c,
                        TMatrixD & def_angle,
                        TMatrixD & def_angle_smear)
{
  Double_t const before_scat = 1;
  Double_t const after_scat = 0;

  Double_t const first_half = 1;
  Double_t const second_half = 2;
  Double_t const through = 2;

  Double_t const Si1 = 1.;
  Double_t const Si6 = 6.;

  Double_t const nphotons = phot_mod.size();

  Double_t const mu = 0;
  Double_t const e = 1;

  //PxPyPzEVector mu_out_mod2 = mu_out_mod;
  //PxPyPzEVector e_out_mod2 = e_out_mod;

  Double_t const d_ini = 1 - 6 * 0.0064; //m lenght of one station without the 6 Si
  Double_t const dSS = 0.01; // m distance between the Si in the couple
  Double_t const d = 0.25 - 0.005; // distance between the couples of SI
  Double_t const dCAL = 0.10; // m distance between the last Si and the Ecal

  Double_t dph0 = 2.025; //m distance for photon interacting in the first target
  Double_t dph1 = 1.025; // distance for photon interacting in the second target


  // the first part is the same for both tar 0 and tar 1

  InitialCond(mu_in_mod , c); // we apply the initial divergence and beam size
  SmearSi(mu_in_mod, c, Si6, before_scat , mu); //we consider the three couples of SI all togethere
  BlankSpace(mu_in_mod, c, d_ini, before_scat , mu); // propagate the muon until the first Be target


  //TAR 0
  // here we distinguish between particles interacting in the first or second target
  if (tar == 0) //interact in the first target
  {

    //FIRST BERILLIUM TARGET - INTERACTION TARGET
    // the muon is attraversing the Be target untill the point of interaction (vertex)
    SmearBe(mu_in_mod, c, first_half, before_scat , mu);

    //rotation of the outgoing particles momenta
    mu_out_mod = RotDiv(mu_in_mod, mu_out_mod);
    e_out_mod = RotDiv(mu_in_mod, e_out_mod);
    def_angle = Def_angle(mu_in_mod, mu_out_mod, e_out_mod);

    // now both electron and muon are going out from the remaining part of the Be target
    SmearBe(mu_out_mod, c, second_half, after_scat , mu);
    SmearBe(e_out_mod, c, second_half, after_scat , e);

    // in one station we have tre couples of Si and a final distance before the next station, we can inplement it
    for (int i = 0; i < 3; i++) //entering the three couples of Be
    {
      //
      BlankSpace(mu_out_mod, c, d, after_scat , mu);
      BlankSpace(e_out_mod, c, d, after_scat , e);

      SmearSi(mu_out_mod, c, Si1, after_scat , mu);
      SmearSi(e_out_mod, c, Si1, after_scat , e);

      BlankSpace(mu_out_mod, c, dSS, after_scat , mu);
      BlankSpace(e_out_mod, c, dSS, after_scat , e);

      SmearSi(mu_out_mod, c, Si1, after_scat , mu);
      SmearSi(e_out_mod, c, Si1, after_scat , e);
    }

    def_angle_smear = Def_angle(mu_in_mod, mu_out_mod, e_out_mod); //smearing dopo le tre coppie di silici

    BlankSpace(mu_out_mod, c, d, after_scat , mu);
    BlankSpace(e_out_mod, c, d, after_scat , e);

    //SECOND BERILLIUM TARGET
    SmearBe(mu_out_mod, c, through, after_scat , mu);
    SmearBe(e_out_mod, c, through, after_scat , e);

    // in one station we have tre couples of Si and a final distance before the next station, we can inplement it
    for (int i = 0; i < 3; i++) //entering the three couples of Si
    {
      //
      BlankSpace(mu_out_mod, c, d, after_scat , mu);
      BlankSpace(e_out_mod, c, d, after_scat , e);

      SmearSi(mu_out_mod, c, Si1, after_scat , mu);
      SmearSi(e_out_mod, c, Si1, after_scat , e);

      BlankSpace(mu_out_mod, c, dSS, after_scat , mu);
      BlankSpace(e_out_mod, c, dSS, after_scat , e);

      SmearSi(mu_out_mod, c, Si1, after_scat , mu);
      SmearSi(e_out_mod, c, Si1, after_scat , e);
    }

    //last distance before the Ecal
    BlankSpace(mu_out_mod, c, dCAL, after_scat , mu);
    BlankSpace(e_out_mod, c, dCAL, after_scat , e);

    //PHOTONS ONLY ROTATION AND STRAIGHT PROPAGATION
    for ( int i = 0; i < nphotons; i++)
    {
      Double_t const ph_id = i + 2;
      phot_mod.at(i) = RotDiv(mu_in_mod, phot_mod.at(i));
      BlankSpace(phot_mod.at(i), c, dph0, after_scat , ph_id);
    }

  }

  //TAR 1
  if (tar == 1) //interact in the second target
  {

    //FIRST BERILLIUM TARGET, the muon passes through
    SmearBe(mu_in_mod, c, through, before_scat , mu);

    // in one station we have tre couples of Si and a final distance before the next station, we can inplement it
    for (int i = 0; i < 3; i++) //entering the three couples of Si only muon
    {
      //
      BlankSpace(mu_in_mod, c, d, before_scat , mu);

      SmearSi(mu_in_mod, c, Si1, before_scat , mu);

      BlankSpace(mu_in_mod, c, dSS, before_scat , mu);

      SmearSi(mu_in_mod, c, Si1, before_scat , mu);
    }

    BlankSpace(mu_in_mod, c, d, after_scat , mu);

    //SECOND BERILLIUM TARGET - INTERACTION TARGET
    // the muon is attraversing the Be target untill the point of interaction (vertex)
    SmearBe(mu_in_mod, c, through, before_scat , mu);

    //rotation of the outgoing particles momenta
    mu_out_mod = RotDiv(mu_in_mod, mu_out_mod);
    e_out_mod = RotDiv(mu_in_mod, e_out_mod);
    def_angle = Def_angle(mu_in_mod, mu_out_mod, e_out_mod);

    // now both electron and muon are going out from the remaining part of the Be target
    SmearBe(mu_out_mod, c, second_half, after_scat , mu);
    SmearBe(e_out_mod, c, second_half, after_scat , e);

    for (int i = 0; i < 3; i++) //entering the three couples of Si all particles
    {
      //
      BlankSpace(mu_out_mod, c, d, after_scat , mu);
      BlankSpace(e_out_mod, c, d, after_scat , e);

      SmearSi(mu_out_mod, c, Si1, after_scat , mu);
      SmearSi(e_out_mod, c, Si1, after_scat , e);

      BlankSpace(mu_out_mod, c, dSS, after_scat , mu);
      BlankSpace(e_out_mod, c, dSS, after_scat , e);

      SmearSi(mu_out_mod, c, Si1, after_scat , mu);
      SmearSi(e_out_mod, c, Si1, after_scat , e);
    }

    def_angle_smear = Def_angle(mu_in_mod, mu_out_mod, e_out_mod); //smearing dopo le tre coppie di silici

    //last distance before the Ecal
    BlankSpace(mu_out_mod, c, dCAL, after_scat , mu);
    BlankSpace(e_out_mod, c, dCAL, after_scat , e);

    //PHOTONS ONLY ROTATION AND STRAIGHT PROPAGATION
    for ( int i = 0; i < nphotons; i++)
    {
      Double_t const ph_id = i + 2;
      phot_mod.at(i) = RotDiv(mu_in_mod, phot_mod.at(i));
      BlankSpace(phot_mod.at(i), c, dph1, after_scat , ph_id);
    }

  }

  // THE PARTICLES HAVE ARRIVED AT THE CALORIMETER



}

// apply the divergence and the beam size of the incoming muon beam
//
void FastSim::InitialCond(PxPyPzEVector & k, TMatrixD & c)
{

  Double_t divthx = gRandom->Gaus(0., 0.00027); // sigmax'=0.00027 rad
  Double_t divthy = gRandom->Gaus(0., 0.00020); // digmay'=0.00020 rad


  Double_t dxdz = tan(divthx); // could approx tan ~ angle
  Double_t dydz = tan(divthy);

  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz * dxdz + dydz * dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;

  k.SetPxPyPzE(skx, sky, skz, k.E());

  Double_t xR = gRandom->Gaus(0, 0.026); // applying beam size in x
  Double_t yR = gRandom->Gaus(0, 0.027); // aplyinng beam size in y

  c[0][0] = xR; // x coordinate for incoming muon
  c[0][1] = yR; // y coordinate for incoming muon
  c[1][0] = xR; // x coordinate for electron, not existing yet
  c[1][1] = yR; // y coordinate for electron, not existing yet
  c[2][0] = xR; // x coordinate for photon 1, not existing yet
  c[2][1] = yR; // y coordinate for photon 1, not existing yet
  c[3][0] = xR; // x coordinate for photon 2, not existing yet
  c[3][1] = yR; // y coordinate for photon 2, not existing yet

}

// apply resolution smearing to particle momentum
//
void FastSim::SmearSi(PxPyPzEVector & k,
                      TMatrixD & c,
                      const Double_t & nSi,
                      const Double_t & before,
                      const Double_t & id )
{

  Double_t const sS  = nSi * 0.00064; //m Si thickness
  Double_t const x0S = 0.094; // m Si interaction lenght
  Double_t sigSI = (13.6 / (k.E() * 1000)) * sqrt(sS / x0S)
                    * (1 + 0.038 * log(sS / x0S)); //rad
  //Double_t thrms = istep==1 ? ThetaRMS(k) : intrinsic_resolution;

  Double_t smearx = gRandom->Gaus(0., sigSI);
  Double_t smeary = gRandom->Gaus(0., sigSI);

  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??

  // apply smearing
  anglex += smearx;
  angley += smeary;

  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);

  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz * dxdz + dydz * dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;

  k.SetPxPyPzE(skx, sky, skz, k.E());

  // how the coordinates change
  Double_t xid = gRandom->Gaus(c[id][0], (1 / sqrt(3)) * sS * sigSI);
  Double_t yid = gRandom->Gaus(c[id][1], (1 / sqrt(3)) * sS * sigSI);

  c[id][0] = xid;
  c[id][1] = yid;

  if (before == 1) // if we are before the scattering, all the particles have the same cooordinates
  {
    c[1][0] = xid;
    c[1][1] = yid;
    c[2][0] = xid;
    c[2][1] = yid;
    c[3][0] = xid;
    c[3][1] = yid;
  }

}


//
void FastSim::SmearBe(PxPyPzEVector & k,
                      TMatrixD & c,
                      const Double_t & f,
                      const Double_t & before,
                      const Double_t & id )
{

  // different possibilities of interaction with the berillium targetù
  //
  Double_t const sBe  = 0.015; //m Be thickness

  Double_t const x0B = 0.353; // m Be interaction lenght
  Double_t sigBe = (13.6 / (k.E() * 1000)) * sqrt(sBe / x0B)
                    * (1 + 0.038 * log(sBe / x0B)); //rad


  // smearing of the angles, gaussian ditribution
  Double_t smearx = gRandom->Gaus(0., sigBe);
  Double_t smeary = gRandom->Gaus(0., sigBe);

  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??

  // apply smearing
  anglex += smearx;
  angley += smeary;

  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);

  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz * dxdz + dydz * dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;

  k.SetPxPyPzE(skx, sky, skz, k.E());

  // how the coordinates change
  Double_t xid = 0;
  Double_t yid = 0;

  if (f == 0) //passing thorugh all the target
  {
    xid = gRandom->Gaus(c[id][0], (1 / sqrt(3)) * sBe * sigBe);
    yid = gRandom->Gaus(c[id][1], (1 / sqrt(3)) * sBe * sigBe);

  }
  else if (f == 1) // scattering target, before the vertex (only muon exist)
  {
    xid = gRandom->Gaus(c[id][0], (1 / sqrt(3)) * vertex * sBe * sigBe);
    yid = gRandom->Gaus(c[id][1], (1 / sqrt(3)) * vertex * sBe * sigBe);

  }
  else if (f == 2) // scattering target, the remaining part after the scattering
  {
    xid = gRandom->Gaus(c[id][0], (1 / sqrt(3)) * (1 - vertex) * sBe * sigBe);
    yid = gRandom->Gaus(c[id][1], (1 / sqrt(3)) * (1 - vertex) * sBe * sigBe);

  }


  c[id][0] = xid;
  c[id][1] = yid;

  if (before == 1) // if we are before the scattering, all the particles have the same cooordinates
  {
    c[1][0] = xid;
    c[1][1] = yid;
    c[2][0] = xid;
    c[2][1] = yid;
    c[3][0] = xid;
    c[3][1] = yid;
  }


}


//
void FastSim::BlankSpace(PxPyPzEVector & k, TMatrixD & c, const Double_t & d,
                         const Double_t & before, const Double_t & id )
{
  //first we have to calculate the angle from the momentum and then propagare the trajectory
  //in the space od distance d
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz());

  Double_t xid = c[id][0];
  Double_t yid = c[id][1];

  c[id][0] = xid + d * tan(anglex);
  c[id][1] = yid + d * tan(angley);

  if (before == 1) // if we are before the scattering, all the particles have the same cooordinates
  {
    c[1][0] = xid;
    c[1][1] = yid;
    c[2][0] = xid;
    c[2][1] = yid;
    c[3][0] = xid;
    c[3][1] = yid;
  }
}


PxPyPzEVector FastSim::RotDiv( PxPyPzEVector & k, PxPyPzEVector & out)
{
  Double_t pz = k.Pz();
  Double_t py = k.Py();
  Double_t px = k.Px();

  Double_t ptz = sqrt(px * px + pz * pz);


  // costruzione della matrice di rotazione
  Double_t psi = atan2(px, pz);
  Double_t phi = atan2(py, ptz);

  TMatrixD R(3, 3);
  R[0][0] = cos(psi);
  R[0][1] = 0;
  R[0][2] = sin(psi);
  R[1][0] = -sin(phi) * sin(psi);
  R[1][1] = cos(phi);
  R[1][2] = sin(phi) * cos(psi);
  R[2][0] = -cos(phi) * sin(psi);
  R[2][1] = -sin(phi);
  R[2][2] = cos(phi) * cos(psi);


  TMatrixD pO(3, 1);
  pO[0][0] = out.Px();
  pO[1][0] = out.Py();
  pO[2][0] = out.Pz();


  TMatrixD pN(R, TMatrixD::kMult, pO);

  PxPyPzEVector pnewdiv(pN[0][0], pN[1][0], pN[2][0], out.E());
  return pnewdiv;
}



TMatrixD FastSim::Def_angle(PxPyPzEVector & p_mu_in,
                            PxPyPzEVector & p_mu_out,
                            PxPyPzEVector & p_e_out)
{

  XYZVector p_mu_in_div3 = p_mu_in.Vect().Unit();
  XYZVector p_e_out_div3 = p_e_out.Vect().Unit();
  XYZVector p_mu_out_div3 = p_mu_out.Vect().Unit();

  Double_t DIR_mu = p_mu_in_div3.Dot(p_mu_out_div3);
  Double_t A_DIR_mu = abs(DIR_mu) < 1. ? 1000.*acos(DIR_mu) : 100.;

  Double_t DIR_e = p_mu_in_div3.Dot(p_e_out_div3);
  Double_t A_DIR_e = abs(DIR_e) < 1. ? 1000.*acos(DIR_e) : 100.;


  TMatrixD def_angle(2, 1);
  def_angle[0][0] = A_DIR_mu;
  def_angle[1][0] = A_DIR_e;

  return def_angle;
}

Double_t FastSim::Def_angle_ph(const PxPyPzEVector & p_mu_in_div,
                               const PxPyPzEVector & p_gamma_lab_div) const
{
  XYZVector p_mu_in_div3 = p_mu_in_div.Vect().Unit();
  XYZVector p_gamma_lab_div3 = p_gamma_lab_div.Vect().Unit();

  Double_t DIR_ph = p_mu_in_div3.Dot(p_gamma_lab_div3);
  Double_t A_DIR_ph = abs(DIR_ph) < 1. ? 1000.*acos(DIR_ph) : 100.;

  return A_DIR_ph;
}

// RMS Theta smearing due to Multiple scattering distribution and intrinsic resolution
//
Double_t FastSim::ThetaRMS(const PxPyPzEVector & k) const
{
  Double_t pmom = k.P();

  // resolution: gaussian sigma
  Double_t thrms(0);

  if (model == 0)
  {
    Double_t msc = thickness > 0 ? MuE::Target_thrms(pmom, thickness) : 0; // in rad
    // N.B. intrinsic_resolution is in mrad
    thrms = twosteps ? msc : sqrt(msc * msc + 1e-6 * intrinsic_resolution * intrinsic_resolution);
  }
  else if (model == 1)
  {
    thrms = MuE::Antonio_thrms(pmom); // in mrad
  }
  else
  {
    cout << "\n" << "***ERROR: Undefined smearing model = " << model << endl;
    exit(999);
  }

  return thrms;
}

// apply resolution smearing to particle momentum
//
PxPyPzEVector FastSim::Smear(const PxPyPzEVector & k, int istep) const
{

  Double_t thrms = istep == 1 ? ThetaRMS(k) : intrinsic_resolution;

  Double_t smearx = gRandom->Gaus(0., thrms);
  Double_t smeary = gRandom->Gaus(0., thrms);

  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??

  // apply smearing
  anglex += smearx;
  angley += smeary;

  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);

  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz * dxdz + dydz * dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;

  PxPyPzEVector psmeared(skx, sky, skz, k.E());

  return psmeared;
}


PxPyPzEVector FastSim::DetEffect(const PxPyPzEVector & p_sim) const
{
  PxPyPzEVector p = p_sim;

  // apply z-scale and uncertainty of station length
  if (zBias != 0 || zSigma != 0) p = z_StationLength(p_sim);

  // check geometric acceptance
  if (p.Theta() > thetaMax_cor)
  {
    //cout<<"theta = "<<p.Theta() <<" > "<<thetaMax_cor<<":  twosteps = "<< twosteps <<endl;
    if (twosteps)
    {
      //      cout<<"\t calling two randoms "<<endl;
      gRandom->Gaus(0., 1.); // dummy calls to preserve the random chain synchronization
      gRandom->Gaus(0., 1.);
    }
    return PxPyPzEVector(0, 0, 0, 0);
  }

  // apply intrinsic detector resolution for two-steps model 0
  if (twosteps)
  {
    p = Smear(p, 2);
    //   cout <<"intrinsic smearing " << endl;
  }

  return p;
}


PxPyPzEVector FastSim::z_StationLength(const PxPyPzEVector & k) const
{
  Polar3DVector p(k.P(), k.Theta(), k.Phi());
  double the = p.Theta();

  double relbias = (zBias * 1e-6) / station_Length; // zBias is in um
  double thebias = -the * relbias;
  double thesm = the + thebias;

  double thesigma = the * (zSigma * 1e-6) / station_Length;

  if (zSigma_switch)
  {
    thesm = thesm +  gRandom->Gaus(0., thesigma);
  }

  double phism = p.Phi();

  // note that theta is defined positive (0-pi)
  // going negative means changing the azimuth too and phi is defined in (-pi,+pi)
  if (thesm < 0)
  {
    thesm = -thesm;
    phism = phism > 0 ? phism - TMath::Pi() : phism + TMath::Pi();
    cout << "*** z_StationLength: thesm < 0  = " << thesm * 1e3 << " mrad" << endl;
  }

  p.SetTheta(thesm);
  p.SetPhi(phism);

  return PxPyPzEVector(p.X(), p.Y(), p.Z(), k.E());
}


// momentum for 2-body kinematics in the centre-of-mass system
//
Double_t FastSim::P_2bodies_CoM (Double_t M, Double_t mm, Double_t me) const
{
  Double_t msum = std::abs(mm + me);
  Double_t mdif = std::abs(mm - me);
  return 0.5 / M * sqrt((M + mdif) * (M - mdif) * (M + msum) * (M - msum));
}

// Lorentz transformation to the centre-of-mass system
//
PxPyPzEVector FastSim::Lorentz_ToCoM(const PxPyPzEVector & pLab) const
{
  if (p_system.E() == Minv) return pLab;
  else
  {
    Double_t ecm = (pLab.E() * p_system.E()
                    - pLab.Px() * p_system.Px() - pLab.Py() * p_system.Py() - pLab.Pz() * p_system.Pz()) / Minv;

    Double_t fn = (ecm + pLab.E()) / (p_system.E() + Minv);

    return PxPyPzEVector(pLab.Px() - fn * p_system.Px(),
                         pLab.Py() - fn * p_system.Py(),
                         pLab.Pz() - fn * p_system.Pz(),
                         ecm);
  }
}

// Lorentz transformation to the Laboratory system
//
PxPyPzEVector FastSim::Lorentz_ToLab(const PxPyPzEVector & pCoM) const
{
  if (p_system.E() == Minv) return pCoM;
  else
  {
    Double_t elab = (pCoM.Px() * p_system.Px() + pCoM.Py() * p_system.Py()
                     + pCoM.Pz() * p_system.Pz() + pCoM.E() * p_system.E()) / Minv;

    Double_t fn   = (elab + pCoM.E()) / (p_system.E() + Minv);

    return PxPyPzEVector(pCoM.Px() + fn * p_system.Px(),
                         pCoM.Py() + fn * p_system.Py(),
                         pCoM.Pz() + fn * p_system.Pz(),
                         elab);
  }
}


// apply resolution smearing to particle momentum (SMEAR only on the XZ plane)
//
PxPyPzEVector FastSim::SmearX(const PxPyPzEVector & k) const
{
  Double_t thrms = ThetaRMS(k);

  Double_t smearx = gRandom->Gaus(0., thrms);
  gRandom->Gaus(0., 1.); // dummy call to preserve the random chain synchronization
  Double_t smeary = 0;

  // angles in the xz and yz planes // defined in -pi, +pi, although should be small angles around zero
  Double_t anglex = atan2(k.Px(), k.Pz());
  Double_t angley = atan2(k.Py(), k.Pz()); // small-angle approx ??

  // apply smearing
  anglex += smearx;
  angley += smeary;

  Double_t dxdz = tan(anglex); // could approx tan ~ angle
  Double_t dydz = tan(angley);

  // assuming z-motion is always forward
  Double_t skz = sqrt(k.P2() / (dxdz * dxdz + dydz * dydz + 1));
  Double_t skx = skz * dxdz;
  Double_t sky = skz * dydz;

  PxPyPzEVector psmeared(skx, sky, skz, k.E());

  return psmeared;
}

// apply resolution smearing to particle momentum
// (theta smearing of rms*sqrt(2) in the plane defined by the vector and the z-axis)
//
PxPyPzEVector FastSim::SmearPolar(const PxPyPzEVector & k) const
{
  Polar3DVector p(k.P(), k.Theta(), k.Phi());

  // assumed smearing is sqrt(2) * smearing on a plane
  Double_t thrms = ThetaRMS(k);
  Double_t smearth = sqrt(2) * gRandom->Gaus(0., thrms);
  gRandom->Gaus(0., 1.); // dummy call to preserve the random chain synchronization
  Double_t thetasm = p.Theta() + smearth;
  Double_t phism = p.Phi();
  // note that theta is defined positive (0-pi)
  // going negative means changing the azimuth too and phi is defined in (-pi,+pi)
  if (thetasm < 0)
  {
    thetasm = -thetasm;
    phism = phism > 0 ? phism - TMath::Pi() : phism + TMath::Pi();
  }

  p.SetTheta(thetasm);
  p.SetPhi(phism);

  return PxPyPzEVector(p.X(), p.Y(), p.Z(), k.E());
}


// synchronize the random number chain to account for events with negligible weight
//  (skipped in the main event loop)
//
void FastSim::RandomNrSync()
{
  gRandom->Gaus(0., 1.); // fake smearing as for the outgoing muon theta X
  gRandom->Gaus(0., 1.); // fake smearing as for the outgoing muon theta Y
  gRandom->Gaus(0., 1.); // fake smearing as for the outgoing electron theta X
  gRandom->Gaus(0., 1.); // fake smearing as for the outgoing electron theta Y
}




void FastSim::PrepareEMShower(KineVars & kv, PxPyPzEVector & e_out_mod,
                              std::vector<PxPyPzEVector> & photons,
                              TMatrixD & c, const Double_t & i,
                              const Double_t & r_mu,
                              GammaFunctionGenerator* & gamma,
                              EMECALShowerParametrization* const & myParam,
                              ECAL* const & myGrid)
{
  bool bFixedLength = true;
  int nPart;
  double X0depth = 0.;
  std::vector<double> energy_in_el;
  std::vector<double> energy_in_ph;

  std::vector<double> coo_el;
  std::vector<double> coo_ph;
  std::vector<double> en_ph_sm;
  //std::vector<EMShower> TheShowerPh;


  auto n_photons = photons.size();

  // Energy and cell of Ecal
  double energy_sm_el = e_out_mod.E();
  double cellEL = myGrid->GiveCentralCell(c[1][0] * 100, c[1][1] * 100);
  double ECAL_E = energy_sm_el;
  vector <double> cellPH;
  //cout<< "pre ciclo"<<endl;
  if (n_photons > 0)
  {
    for (unsigned int j = 0; j < n_photons; j++)
    {
      cellPH.push_back(myGrid->GiveCentralCell(c[j + 2][0] * 100, c[j + 2][1] * 100));
      //cout<< "pre ciclo"<<endl;
      PxPyPzEVector phot = photons.at(j);
      en_ph_sm.push_back(phot.E());
      //cout<< "pre ciclo"<<endl;
      ECAL_E += en_ph_sm.at(j);
    }
  }
  //cout<<"quii3"<<endl;
  if (r_mu < 5 && ECAL_E > 0.2)
  {
    // for electrons
    if (energy_sm_el > 0.2 && cellEL != 0)
    {
      nPart = 1;
      X0depth = 0.;
      coo_el.push_back(c[1][0] * 100); // cm
      coo_el.push_back(c[1][1] * 100); // cm
      energy_in_el.push_back(energy_sm_el);
      myGrid->SetEnergy(energy_sm_el);
      EMShower TheShowerEl(gamma, myParam, myGrid, bFixedLength,
                           nPart, X0depth, energy_in_el, coo_el);
      TheShowerEl.compute();
      coo_el.clear();
      energy_in_el.clear();
      //cout<<"quii3"<<endl;

      if (n_photons > 20)//HO MODIFICATO PERCHE' NON VOGLIO VEDERE I FOTONI NEL CALO
      {
        for (unsigned int j = 0; j < n_photons; j++)
        {
          if ( en_ph_sm.at(j) > 0.2 && cellPH.at(j) != 0)
          {
            nPart = 2;
            X0depth = -log(gRandom->Uniform()) * (9. / 7.);
            //cout<< "pre ciclo"<<endl;
            coo_ph.push_back(c[j + 2][0] * 100); // cm
            coo_ph.push_back(c[j + 2][1] * 100); // cm
            //cout<< "pre ciclo"<<endl;
            energy_in_ph.push_back(en_ph_sm.at(j) / 2);
            energy_in_ph.push_back(en_ph_sm.at(j) / 2);
            //cout<< "pre ciclo"<<endl;
            myGrid->SetEnergy(en_ph_sm.at(j));
            EMShower TheShowerPh(gamma, myParam, myGrid, bFixedLength, nPart,
                                 X0depth, energy_in_ph, coo_ph);
            TheShowerPh.compute();
            coo_ph.clear();
            energy_in_ph.clear();
          }
        }
      }
    }
    LoadECAL(detKin, myGrid, i);
  }
  else
  {
    LoadECAL(detKin, myGrid, i);
  }
  cellPH.clear();
  en_ph_sm.clear();
}

void FastSim::LoadECAL(KineVars & kv, ECAL* const & myGrid, int j)
{
  TH2F* Ecal = myGrid->GiveEcalGrid();
  if (Ecal->GetEntries() != 0)
  {
    double *Ecell = myGrid->EnergyContent();
    /*double *E_clus = myGrid->Draw_ECAL(j);
    kv.n_max_Cell=E_clus[0];
    kv.E_clus3x3=E_clus[1];
    kv.E_1=E_clus[2];*/


    kv.Ecell1 = Ecell[0];
    kv.Ecell2 = Ecell[1];
    kv.Ecell3 = Ecell[2];
    kv.Ecell4 = Ecell[3];
    kv.Ecell5 = Ecell[4];
    kv.Ecell6 = Ecell[5];
    kv.Ecell7 = Ecell[6];
    kv.Ecell8 = Ecell[7];
    kv.Ecell9 = Ecell[8];
    kv.Ecell10 = Ecell[9];
    kv.Ecell11 = Ecell[10];
    kv.Ecell12 = Ecell[11];
    kv.Ecell13 = Ecell[12];
    kv.Ecell14 = Ecell[13];
    kv.Ecell15 = Ecell[14];
    kv.Ecell16 = Ecell[15];
    kv.Ecell17 = Ecell[16];
    kv.Ecell18 = Ecell[17];
    kv.Ecell19 = Ecell[18];
    kv.Ecell20 = Ecell[19];
    kv.Ecell21 = Ecell[20];
    kv.Ecell22 = Ecell[21];
    kv.Ecell23 = Ecell[22];
    kv.Ecell24 = Ecell[23];
    kv.Ecell25 = Ecell[24];

    delete Ecal;
  }
  else
  {
    /*kv.n_max_Cell=0;
    kv.E_clus3x3=0;
    kv.E_1=0.;*/

    kv.Ecell1 = 0;
    kv.Ecell2 = 0;
    kv.Ecell3 = 0;
    kv.Ecell4 = 0;
    kv.Ecell5 = 0;
    kv.Ecell6 = 0;
    kv.Ecell7 = 0;
    kv.Ecell8 = 0;
    kv.Ecell9 = 0;
    kv.Ecell10 = 0;
    kv.Ecell11 = 0;
    kv.Ecell12 = 0;
    kv.Ecell13 = 0;
    kv.Ecell14 = 0;
    kv.Ecell15 = 0;
    kv.Ecell16 = 0;
    kv.Ecell17 = 0;
    kv.Ecell18 = 0;
    kv.Ecell19 = 0;
    kv.Ecell20 = 0;
    kv.Ecell21 = 0;
    kv.Ecell22 = 0;
    kv.Ecell23 = 0;
    kv.Ecell24 = 0;
    kv.Ecell25 = 0;

    delete Ecal;
  }
}