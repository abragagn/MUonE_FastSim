/////////////////////////////////////////////////////////////////////////////
//
// Generate MESMER MuE MC events and process them through  
//  Fast Detector Simulation and Analysis in one step
//  saving relevant histograms, trees and plots
//
//   Input parameters: filenames for
//   1: input MESMER card file 
//   2: input MuE cfg file
//
// G.Abbiendi  24/Feb/2020 first version: read MESMER events from ascii input
//     updated 17/Apr/2020 
// G.A.+C.Carloni Calame 17/Mar/2022: embedded MESMER generation 
/////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TChain.h"
#include "TRint.h"
#include "TSystem.h"
#include "TTree.h"
#include "FastSim.h"
#include "Inputs.h"
#include "Analysis.h"
#include "TRandom.h"
#include "Mesmer.h"
#include "ECALProprieties.h"
#include "GammaFunctionGenerator.h"
#include "ECAL.h"

using namespace std;
using namespace MuE;

int main(int argc, char* argv[])
{

  if (argc != 3)
  {
    cerr << "Usage : "<< argv[0] << " input_mesmer_card_file  input_MuE_cfg_file \n";
    exit(100);
  }

  for(int i = 0; i < argc; i++)
    printf("%s \n", argv[i]);
    
  // set it to true to debug
  bool _debug_ = false;
  bool _debug_inp_ = true;

  //  read MESMER generator parameters
  MuE::MCpara pargen;
  pargen.InitandSetRunParams_mesmer(argv[1]);

  // read MuE configuration parameters
  MuE::One_Input cfg;
  cfg.configure(argv[2], _debug_inp_);
  //  cfg.configure("input.cfi", _debug_inp_);
  MuE::FS_Input fsi;
  fsi.configure(cfg.fastSim_ifname, _debug_inp_);
  MuE::AN_Input an_input;
  an_input.configure(cfg.analysis_ifname, _debug_inp_);

  TString workdir(cfg.output_dir);
  gSystem->mkdir(workdir);
  //  gSystem->cd(workdir);
  cerr << "Writing analysis results to output directory : " << workdir << endl;

  // Set it to kTRUE if you do not run interactively
  gROOT->SetBatch(kTRUE);
  // gROOT->SetBatch(kFALSE);

  // Initialize Root application
  auto app = new TRint("Root Application", &argc, argv);

  //................................................................
  // Initialise the muon and electron masses used in class Particle
  // consistently with those used by the generator
  Particle::set_mass_e(pargen.mass_e);
  Particle::set_mass_mu(pargen.mass_mu);
  //................................................................

  // Fast Simulation //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MuE::FastSim fs(pargen, fsi, _debug_);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ANALYSIS //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MuE::Analysis analyzer(cfg, pargen, fsi, an_input);
  analyzer.BeginJob();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //////////////////////////////////////////////////////////////////////////////
  ////////////////////////////   EVENT  LOOP   /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  MuE::Event event;

  // number of requested events
  const Long_t HowManyEvts = pargen.Nevreq;

  // Event generation loop
  cout << "Start generating MESMER events and processing them through MuE FastSim and Analysis" << endl;
  Long_t nev = 0;

  // number of events with negligible weight (skipped)
  int zero_wgt_events = 0;
  auto gamma = new GammaFunctionGenerator;
  auto ecalprop = new ECALProperties();
  auto myparam = new EMECALShowerParametrization(ecalprop,
                                                 {100.0, 0.1},
                                                 {1.0, 0.1, 100.0, 1.0},
                                                 1, 1);
  auto TheEcal = new ECAL(5, -7.125, 7.125, 5, -7.125, 7.125);

  while (nev < HowManyEvts)
  {

    double pmu[4];
    IncomingMuonMomentum_mesmer(pmu);
    // gets the initial state muon 4-momentum (E,px,py,pz) components
    // can be replaced by any 4-momentum "provider" function

    int ierr = event.GenerateEvent_mesmer(pmu);

    if (event.EventNr % 1000000 == 0 )
      cerr << "\n processing event : " << event.EventNr << "\r" << flush;

    //cout<< "creato ecalecc"<< endl;

    if (ierr == 0)
    {
      nev++;
      ////////////////////////////////////////
      fs.Process(event, gamma, myparam, TheEcal);
      //cout<< "process"<< endl;
      ////////////////////////////////////////
      analyzer.Analyze(event, fs);
      //cout<< "analize"<< endl;
      ////////////////////////////////////////
    }
    else
    {
      // synchronize the random number chain to the event loop to be reproducible
      zero_wgt_events++;
      fs.RandomNrSync();
      // event.Print();
    }
    //cout<< "predebu"<< endl;
    if (_debug_) event.Print();
    //cout<< "debug"<< endl;

  }


  cout << "End generating MESMER events." << endl;
  cout << "\n Processed " << nev << " events." << endl;
  cout << " number of zero-weight events = " << zero_wgt_events << endl;

  cout << endl << " last Random seed = " << gRandom->GetSeed() << endl;
  cout << "last call to Rndm() = " << setw(20) << setprecision(17) << gRandom->Rndm() << endl;

  // Fill summary MC statistics
  MCstat sums;
  sums.SetEndofRun_mesmer();

  // printout the summary MC statistics
  sums.Print();

  // fill the parameter tree (1 entry per mc run) and save it to output
  gFile->cd();
  cout << "\n" << "writing MC parameter tree to file: ";
  gDirectory->pwd();
  cout << endl;
  Setup setup(pargen, sums);
  Int_t splitBranches = 2;
  auto partree = new TTree("MuEsetup", "MuE MC parameters");
  partree->Branch("MuEparams", &setup, 64000, splitBranches);

  partree->Fill();
  if (_debug_) partree->Print();
  analyzer.EndJob(sums, fs);

  //delete event;//?

  cout << "Finished." << endl;

  if (!gROOT->IsBatch()) app->Run();

  return 0;

}
