#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"

const int showerSize_ = 1;        // Shower size around the crystal with max energy deposited
                                  // for point of impact reconstruction (i.e. 1 -> 3x3 matrix)
const double w0_ = 4.;            // Parameter for point of impact reconstruction, 
                                  // It determine the minimum Ecell / Etot to be used
const double crystalSize = 2.85;  // Crystal size 


void ECAL_analysis(TString fname = "outtree.root")
{
  auto f = new TFile(fname);
  auto t = (TTree*)f->Get("atree");

  int ne = t->GetEntries();

  auto hReso = new TH1D("hReso", "Resolution;Resolution [%];", 100, 0., 50.);
  auto hR = new TH1D("hR", "r;r [cm];", 100, 0., 5 * crystalSize);

  double det_Ee, gen_Ee; // Dovrebbero essere le energie dell'elettrone generate 
                         // e da "detector", non so bene come siano calcolate
  double Ecell[25]; // Crystals array
  double n_cell_e;  // Cell hit by the electron

  for (int i = 0; i < 25; ++i) Ecell[i] = 0;

  t->SetBranchAddress("event.detKin.Ee", &det_Ee);
  t->SetBranchAddress("event.genKin.Ee", &gen_Ee);
  t->SetBranchAddress("event.detKin.n_cell_e", &n_cell_e);
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
    for (int y = 0; y < 5; ++y)
      for (int x = 0; x < 5; ++x)
        ECALmatrix[x][y] = Ecell[y * 5 + x];

    // Total energy deposited
    double etot = 0;
    for (int y = 0; y < 5; ++y)
      for (int x = 0; x < 5; ++x)
        etot += ECALmatrix[x][y];

    if (etot == 0) continue;

    // Evaluate energy resolution (?)
    double reso = abs(etot - gen_Ee) / gen_Ee;
    hReso->Fill(100.*reso);

    // Estimate point of impact (poi)
    double x_poi, y_poi;

    // Search cell with higher energy
    double E_max = 0.;
    int x_max, y_max;

    for (int y = 0; y < 5; ++y)
    {
      for (int x = 0; x < 5; ++x)
      {
        if (ECALmatrix[x][y] > E_max)
        {
          E_max = ECALmatrix[x][y];
          x_max = x;
          y_max = y;
        }
      }
    }

    // cout << "Emax cell: [" << x_max << ", " << y_max <<  "]" << endl;

    // Reconstruct PoI
    // Loops constrained by ECAL boundaries

    int x_low  = max(0, x_max - showerSize_);
    int x_high = min(x_max + showerSize_, 4);
    int y_low  = max(0, y_max - showerSize_);
    int y_high = min(y_max + showerSize_, 4);

    // cout << "Shower x size: " << x_low << " -> " << x_high << endl;
    // cout << "Shower y size: " << y_low << " -> " << y_high << endl;

    // First calculate total energy shower
    double showerEn = 0.;

    for (int y = y_low; y < y_high; ++y)
      for (int x = x_low; x < x_high; ++x)
        showerEn += ECALmatrix[x][y];

    // cout << "Shower energy = " << showerEn << endl;

    // Then estimate shower center
    double x_shower = 0.;
    double y_shower = 0.;
    double w_norm = 0.; // Sum of weights for normalization
    for (int y = y_low; y <= y_high; ++y)
    {
      for (int x = x_low; x <= x_high; ++x)
      {
        double en = ECALmatrix[x][y]; // energy deposited in the [j][k] cell
        double w = TMath::Max(0., w0_ + log(en / showerEn)); // weight of the cell
        w_norm += w;

        double cellX = (-2 + x) * crystalSize;
        double cellY = ( 2 - y) * crystalSize;

        x_shower += w * cellX; 
        y_shower += w * cellY;

        // cout << "Cell [" << x << ", " << y << "] Pos ("
        //      << cellX << ", " << cellY << ")"
        //      << " with En " << en << " and weight " << w << endl;
        // cout << "    Updated xPos = " << x_shower
        //      << ", yPos = " << y_shower << endl <<endl;
      }
    }

    x_shower /= w_norm; // Normalize to sum of weights
    y_shower /= w_norm;

    // Measure theta (I do not know what to use for the distance,
    // i.e. the z coordinate of the scattering)
    double r = sqrt(pow(x_shower, 2) + pow(y_shower, 2));
    // double dist = ...
    // double theta = TMath::ATan(r / dist);

    hR->Fill(r);

    // cout << " Est. PoI : (" << x_shower << ", " << y_shower << ") in "
    //      << "[" << x_max << ", " << y_max << "] corresponding to cell "
    //      << n_cell_e << endl;

  }

  hReso->Draw();
  hR->Draw();
}