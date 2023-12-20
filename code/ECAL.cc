#include <TH2F.h>
#include <map>
#include <cmath>

#include <iostream>
#include <sstream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>



#include "ECAL.h"


using namespace std;


ECAL::ECAL(double nbinsx, double xlow, double xup,
           double nbinsy, double ylow, double yup) :
            nbinX(nbinsx), nbinY(nbinsy), Xlow(xlow),
            Xup(xup), Ylow(ylow), Yup(yup) 
{
  //Queste mappe servono a mappare il numero di bin nel numero vero della cella e viceversa
  //perchè i numeri dei bin sono sballati a causa degli overflow e underflow bins
  number[36] = 1;
  number[37] = 2;
  number[38] = 3;
  number[39] = 4;
  number[40] = 5;
  number[29] = 6;
  number[30] = 7;
  number[31] = 8;
  number[32] = 9;
  number[33] = 10;
  number[22] = 11;
  number[23] = 12;
  number[24] = 13;
  number[25] = 14;
  number[26] = 15;
  number[15] = 16;
  number[16] = 17;
  number[17] = 18;
  number[18] = 19;
  number[19] = 20;
  number[8] = 21;
  number[9] = 22;
  number[10] = 23;
  number[11] = 24;
  number[12] = 25;

  Rev_number[1] = 36;
  Rev_number[2] = 37;
  Rev_number[3] = 38;
  Rev_number[4] = 39;
  Rev_number[5] = 40;
  Rev_number[6] = 29;
  Rev_number[7] = 30;
  Rev_number[8] = 31;
  Rev_number[9] = 32;
  Rev_number[10] = 33;
  Rev_number[11] = 22;
  Rev_number[12] = 23;
  Rev_number[13] = 24;
  Rev_number[14] = 25;
  Rev_number[15] = 26;
  Rev_number[16] = 15;
  Rev_number[17] = 16;
  Rev_number[18] = 17;
  Rev_number[19] = 18;
  Rev_number[20] = 19;
  Rev_number[21] = 8;
  Rev_number[22] = 9;
  Rev_number[23] = 10;
  Rev_number[24] = 11;
  Rev_number[25] = 12;
  

//  Energy_dist =new TH1F("Energy", "Energy",100,90,100);
// Energy_dist1 =new TH1F("Energy", "Energy 1 cell",200,0.20,1);
// Energy_dist3x3 =new TH1F("Energy", "Energy 3x3 cells",200,0.70,1);


  Array9 = 0;

}

// metodo che crea l'istogramma rappresentante il calorimetro

/*TH2F* ECAL::CreateGrid(double nbinsx,double xlow,double xup,double nbinsy,double ylow,double yup)
{
    TH2F* EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter with E in GeV",nbinsx,xlow,xup,nbinsy,ylow,yup);
    return EcalGrid;
};*/

void ECAL::CreateGrid(double nbinsx, double xlow, double xup,
                      double nbinsy, double ylow, double yup)
{
  EcalGrid = new TH2F("EcalGrid" , "EM Calorimeter with E in GeV",
                      nbinsx, xlow, xup, nbinsy, ylow, yup);
  // return EcalGrid;
};


TH2F* ECAL::GiveEcalGrid()
{
  return EcalGrid;
};

void ECAL::SetEnergy(double energy)
{
  energy_IN = energy;
}

// metodo che assegna il numero della cella che viene colpita dalla particella
double ECAL::GiveCentralCell(double coox, double cooy)
{
  int binx = EcalGrid->GetXaxis()->FindBin(coox);
  int biny = EcalGrid->GetYaxis()->FindBin(cooy);
  int nbin = EcalGrid->GetBin(binx, biny);

  //cout <<"Number of the cell:" << number[nbin] << endl;

  return number[nbin];
};

/*int* ECAL::GiveArray3x3(int n)
{
    if (n==1) {Array9= new int[9]{1,6,7,2,0,0,0,0,0};}
    if (n==2) {Array9= new int[9]{1,2,6,7,8,3,0,0,0};}
    if (n==3) {Array9= new int[9]{2,3,7,8,9,4,0,0,0};}
    if (n==4) {Array9= new int[9]{3,4,8,9,10,0,0,0};}
    if (n==5) {Array9= new int[9]{4,5,9,10,0,0,0};}
    if (n==6) {Array9= new int[9]{1,2,6,7,12,11,0,0,0};}
    if (n==7) {Array9= new int[9]{1,2,3,6,7,8,11,12,13};}
    if (n==8) {Array9= new int[9]{2,3,4,7,8,9,12,13,14};}
    if (n==9) {Array9= new int[9]{3,4,5,8,9,10,13,14,15};}
    if (n==10) {Array9= new int[9]{4,5,9,10,14,15,0,0,0};}
    if (n==11) {Array9= new int[9]{6,7,11,12,16,17,0,0,0};}
    if (n==12) {Array9= new int[9]{6,7,8,11,12,13,16,17,18};}
    if (n==13) {Array9= new int[9]{7,8,9,12,13,14,17,18,19};}
    if (n==14) {Array9= new int[9]{8,9,10,13,14,15,18,19,20};}
    if (n==15) {Array9= new int[9]{9,10,14,15,19,20,0,0,0};}
    if (n==16) {Array9= new int[9]{11,12,16,17,21,22,0,0,0};}
    if (n==17) {Array9= new int[9]{11,12,13,16,17,18,21,22,23};}
    if (n==18) {Array9= new int[9]{12,13,14,17,18,19,22,23,24};}
    if (n==19) {Array9= new int[9]{13,14,15,18,19,20,23,24,25};}
    if (n==20) {Array9= new int[9]{14,15,19,20,24,25,0,0,0};}
    if (n==21) {Array9= new int[9]{16,17,21,22,0,0,0,0,0};}
    if (n==22) {Array9= new int[9]{16,17,18,21,22,23,0,0,0};}
    if (n==23) {Array9= new int[9]{17,18,19,22,23,24,0,0,0};}
    if (n==24) {Array9= new int[9]{18,19,20,23,24,25,0,0,0};}
    if (n==25) {Array9= new int[9]{19,20,24,25,0,0,0,0,0};}
    
    return 0;
}*/

// metodo che aggiunge il punto di coo(x,y) all'istogramma, quindi al calorimetro e dà numero cella


double ECAL::AddHitCoo(double r, double phi, double xi, double yi, double w)
{
  r *= 2.19;
  double x = r * cos(phi) + xi; // coo x in cm
  double y = r * sin(phi) + yi; // coo y in cm
  EcalGrid->Fill(x, y, w);

  double number = ECAL::GiveCentralCell(x, y);
  return number;
};

void ECAL::AddHitCooDepth(double r, double phi, double xi, double yi,
                          double w, double depth, double X0depth)
{
  depth += X0depth;
  r *= 2.19;
  double x = r * cos(phi) + xi; // coo x in cm
  double y = r * sin(phi) + yi; // coo y in cm
  if (24.7 - X0depth > depth)
  {
    EcalGrid->Fill(x, y, w);
  }
};

// metodo che disegna l'evento nel calorimetro e le celle che vengono colpite


/*TCanvas * Ecal_= new TCanvas("Ecal_","Ecal_",1500,100,3500,2000);
Ecal_->Divide(2,1);
Ecal_->cd(1);
gStyle->SetPalette(kAquamarine);
//TColor::InvertPalette();
EcalGrid->SetXTitle("x (cm)");
EcalGrid->SetYTitle("y (cm)");
EcalGrid->Draw("COL");
EcalGrid->Draw("TEXT SAME");
Ecal_->cd(2);
EcalGrid->Draw("LEGO");
std::ostringstream name1;
name1 <<"/home/LHCB-T3/espedicato/tesi/ECALpng/Ecal"<< i << ".png";
TString name =name1.str();
Ecal_->SaveAs(name);*/

/*double* ECAL::Draw_ECAL(int i){
double* ECluster = new double[3];
// riempi celle    
int binMax=EcalGrid->GetMaximumBin();  
int CentralCell=number[binMax];
cout << "cella centrale rev " << Rev_number[CentralCell] <<" and vera " << CentralCell << endl;
//Energy_dist1->Fill(EcalGrid->GetBinContent(binMax)/energy_IN);
ECluster[0]=((double)CentralCell); 

double energy3x3=0.;    
ECAL::GiveArray3x3(CentralCell);
for (int i=0; i<9; ++i)
{
    if (Array9[i]>0 && Array9[i]<26) energy3x3+=EcalGrid->GetBinContent(Rev_number[Array9[i]]);
    cout << Rev_number[Array9[i]] << " and vera " << Array9[i]<< " c'è energia " << EcalGrid->GetBinContent(Rev_number[Array9[i]]) << endl;
}
//Energy_dist3x3->Fill(energy3x3/energy_IN);
ECluster[1]=(energy3x3); 
ECluster[2]=(EcalGrid->GetBinContent(binMax)); 

return ECluster;
};*/


double* ECAL::EnergyContent()
{
  double* E_cell = new double[25];
  for (int i = 1; i < 26 ; ++i)
  {
    E_cell[i - 1] = (EcalGrid->GetBinContent(Rev_number[i]));
  }
  return E_cell;
}

/*void ECAL::Print_()
{TCanvas * encell= new TCanvas("Energy cells","Energy cells",1000,100,2500,2000);
encell->Divide(1,2);
encell->cd(1);
Energy_dist1->SetLineWidth(3);
Energy_dist1->GetXaxis()->SetTitle("E_rec/E_in");
Energy_dist1->Draw();
encell->cd(2);
Energy_dist3x3->SetLineWidth(3);
Energy_dist3x3->GetXaxis()->SetTitle("E_rec/E_in");
Energy_dist3x3->Draw();
encell->SaveAs("/home/LHCB-T3/espedicato/tesi/EnCell.png");
}*/