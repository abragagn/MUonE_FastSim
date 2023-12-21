#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"

void ECAL_analysis(TString fname = "outtree.root")
{
  auto f = new TFile(fname);
  auto t = (TTree*)f->Get("atree");

  int ne = t->GetEntries();

  auto h = new TH1D("h", "Resolution;Resolution;[%]", 100, 0., 50.);

  double det_Ee, gen_Ee; // Dovrebber essere le energie dell'elettrone, non so bene come siano calcolate
  double Ecell[25]; // Crystals array

  for (int i = 0; i < 25; ++i) Ecell[i] = 0;

  t->SetBranchAddress("event.detKin.Ee", &det_Ee);
  t->SetBranchAddress("event.genKin.Ee", &gen_Ee);
  t->SetBranchAddress("event.detKin.Ecell1", &Ecell[0]);
  t->SetBranchAddress("event.detKin.Ecell2", &Ecell[1]);
  t->SetBranchAddress("event.detKin.Ecell3", &Ecell[2]);
  t->SetBranchAddress("event.detKin.Ecell4", &Ecell[3]);
  t->SetBranchAddress("event.detKin.Ecell5", &Ecell[4]);
  t->SetBranchAddress("event.detKin.Ecell6", &Ecell[5]);
  t->SetBranchAddress("event.detKin.Ecell7", &Ecell[6]);
  t->SetBranchAddress("event.detKin.Ecell8", &Ecell[7]);
  t->SetBranchAddress("event.detKin.Ecell9", &Ecell[8]);
  t->SetBranchAddress("event.detKin.Ecell10", &Ecell[9]);
  t->SetBranchAddress("event.detKin.Ecell11", &Ecell[10]);
  t->SetBranchAddress("event.detKin.Ecell12", &Ecell[11]);
  t->SetBranchAddress("event.detKin.Ecell13", &Ecell[12]);
  t->SetBranchAddress("event.detKin.Ecell14", &Ecell[13]);
  t->SetBranchAddress("event.detKin.Ecell15", &Ecell[14]);
  t->SetBranchAddress("event.detKin.Ecell16", &Ecell[15]);
  t->SetBranchAddress("event.detKin.Ecell17", &Ecell[16]);
  t->SetBranchAddress("event.detKin.Ecell18", &Ecell[17]);
  t->SetBranchAddress("event.detKin.Ecell19", &Ecell[18]);
  t->SetBranchAddress("event.detKin.Ecell20", &Ecell[19]);
  t->SetBranchAddress("event.detKin.Ecell21", &Ecell[20]);
  t->SetBranchAddress("event.detKin.Ecell22", &Ecell[21]);
  t->SetBranchAddress("event.detKin.Ecell23", &Ecell[22]);
  t->SetBranchAddress("event.detKin.Ecell24", &Ecell[23]);
  t->SetBranchAddress("event.detKin.Ecell25", &Ecell[24]);

  cout << "Analyzing " << ne << " events" << endl;

  for (int i = 0; i < ne; ++i)
  {
    t->GetEntry(i);

    // Create ecal matrix
    double ECALmatrix[5][5];
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 5; ++k)
        ECALmatrix[j][k] = Ecell[j * 5 + k];

    // Total energy deposited
    double etot = 0;
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 5; ++k)
        etot += ECALmatrix[j][k];

    // Evaluate energy resolution (?)
    if (etot != 0)
    {
      double reso = abs(etot-gen_Ee)/gen_Ee;
      h->Fill(100.*reso);
    } else
    {
      continue; // Ele not in calo
    }

    // Estimate point of impact
    // Search cell with higher energy
    double E_max = 0.;
    int x_max, y_max;

    for (int j = 0; j < 5; ++j)
    {
      for (int k = 0; k < 5; ++k)
      {
        if (ECALmatrix[j][k] > E_max)
        {
          E_max = ECALmatrix[j][k];
          y_max = j;
          x_max = k;
        }
      }
    }

    // Evaluate point of impact

  }

  h->Draw();
}