#include "EMShower.h"

//#include "RandomEngineAndDistribution.h"
#include "GammaFunctionGenerator.h"
//#include <SpecFuncCephes.h>

#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TProfile.h>

//#include "FastSimulation/Utilities/interface/Histos.h"

using std::vector;

EMShower::EMShower(GammaFunctionGenerator *gamma,
                   EMECALShowerParametrization *const myParam,
                   ECAL *const myGrid,
                   bool bFixedLength,
                   int nPart,
                   double X0depth,
                   vector<double> energy_in,
                   vector<double> coo)

    : myGammaGenerator(gamma),
      theParam(myParam),
      // thePart(myPart),
      theGrid(myGrid),
      // random(engine),
      bFixedLength_(bFixedLength),
      nPart(nPart),
      X0depth(X0depth)
{

  stepsCalculated = false;

  theECAL = theParam->ecalProperties();

  double fotos = theECAL->photoStatistics() * theECAL->lightCollectionEfficiency();
  //--->double fotos = 50.E3 * 0.03;

  // nPart = thePart->size();
  totalEnergy = 0.;
  globalMaximum = 0.;
  double meanDepth = 0.;
  // Initialize the shower parameters for each particle
  Xi = coo[0];
  Yi = coo[1];

  for (int i = 0; i < nPart; ++i)
  {

    // The particle and the shower energy

    Etot.push_back(0.);
    E.push_back(energy_in[i]);
    totalEnergy += E[i];
    double lny = std::log(E[i] / theECAL->criticalEnergy());
    // double lny = std::log(E[i] / 8.74E-3);

    // Average and Sigma for T and alpha
    double theMeanT = myParam->meanT(lny);
    double theMeanAlpha = myParam->meanAlpha(lny);
    double theMeanLnT = myParam->meanLnT(lny);
    double theMeanLnAlpha = myParam->meanLnAlpha(lny);
    double theSigmaLnT = myParam->sigmaLnT(lny);
    double theSigmaLnAlpha = myParam->sigmaLnAlpha(lny);
    // The correlation matrix
    double theCorrelation = myParam->correlationAlphaT(lny);
    double rhop = std::sqrt((1. + theCorrelation) / 2.);
    double rhom = std::sqrt((1. - theCorrelation) / 2.);

    // The number of spots in ECAL / HCAL
    theNumberOfSpots.push_back(theParam->nSpots(E[i]));

    // Photo-statistics
    photos.push_back(E[i] * fotos);

    // The longitudinal shower development parameters
    // Fluctuations of alpha, T and beta
    double z1 = 0.;
    double z2 = 0.;
    double aa = 0.;

    // Protect against too large fluctuations (a < 1) for small energies
    while (aa <= 1.)
    {
      z1 = gRandom->Gaus(0., 1.);
      z2 = gRandom->Gaus(0., 1.);
      aa = std::exp(theMeanLnAlpha + theSigmaLnAlpha * (z1 * rhop - z2 * rhom));
    }

    a.push_back(aa);
    T.push_back(std::exp(theMeanLnT + theSigmaLnT * (z1 * rhop + z2 * rhom)));
    b.push_back((a[i] - 1.) / T[i]);
    maximumOfShower.push_back((a[i] - 1.) / b[i]);
    globalMaximum += maximumOfShower[i] * E[i];
    meanDepth += a[i] / b[i] * E[i];

    Ti.push_back(a[i] / b[i] * (std::exp(theMeanLnAlpha) - 1.) / std::exp(theMeanLnAlpha));

    // The parameters for the number of energy spots
    TSpot.push_back(theParam->meanTSpot(theMeanT));
    aSpot.push_back(theParam->meanAlphaSpot(theMeanAlpha));
    bSpot.push_back((aSpot[i] - 1.) / TSpot[i]);
  }

  globalMaximum /= totalEnergy;
  meanDepth /= totalEnergy;

  // Se voglio studiare cluster separati:
  // theGrid->CreateGrid(5,-7.125,7.125,5,-7.125,7.125);
  // Se voglio studiare cluster evento (e+gamma) insieme:
  // theGrid->GiveEcalGrid();
}

void EMShower::prepareSteps()
{
  double dt;
  double radlen;
  int stps;
  int first_Ecal_step = 0;
  int last_Ecal_step = 0;

  // The maximum is in principe 8 (with 5X0 steps in the ECAL)
  steps.reserve(24);

  // radlen = -theGrid->x0DepthOffset();

  // ECAL
  radlen = 24.7; // 22cm/0.89 cm

  if (radlen > 0.)
  {
    if (!bFixedLength_)
    {
      stps = (int)((radlen + 2.5) / 5.);
      if (stps == 0)
        stps = 1;
      dt = radlen / (double)stps;
      Step step(0, dt);
      first_Ecal_step = steps.size();
      for (int ist = 0; ist < stps; ++ist)
        steps.push_back(step);
      last_Ecal_step = steps.size() - 1;
      radlen = 0.;
    }
    else
    {
      dt = 1.0;
      stps = static_cast<int>(radlen);
      if (stps == 0)
        stps = 1;
      Step step(0, dt);
      first_Ecal_step = steps.size();
      for (int ist = 0; ist < stps; ++ist)
        steps.push_back(step);
      dt = radlen - stps;
      if (dt > 0)
      {
        Step stepLast(2, dt);
        steps.push_back(stepLast);
      }
      last_Ecal_step = steps.size() - 1;

      radlen = 0.;
    }
  }

  nSteps = steps.size();
  if (nSteps == 0)
    return;
  double ESliceTot = 0.;
  double MeanDepth = 0.;
  depositedEnergy.resize(nSteps);
  meanDepth.resize(nSteps);
  double t = 0.;

  int offset = 0;
  for (unsigned iStep = 0; iStep < nSteps; ++iStep)
  {
    ESliceTot = 0.;
    MeanDepth = 0.;
    double realTotalEnergy = 0;
    dt = steps[iStep].second;
    t += dt;
    for (int i = 0; i < nPart; ++i)
    {
      depositedEnergy[iStep].push_back(deposit(t, a[i], b[i], dt));
      // cout << " % energia depositata allo step " << iStep << " è " << depositedEnergy[iStep][i] << endl;
      ESliceTot += depositedEnergy[iStep][i];
      MeanDepth += deposit(t, a[i] + 1., b[i], dt) / b[i] * a[i];

      realTotalEnergy += depositedEnergy[iStep][i] * E[i];
    }

    if (ESliceTot > 0.) // can happen for the shower tails; this depth will be skipped anyway
      MeanDepth /= ESliceTot;
    else
      MeanDepth = t - dt;

    meanDepth[iStep] = MeanDepth;

    if (realTotalEnergy < 0.001)
    {
      offset -= 1;
    }
    // cout << " energia depositata allo step " << iStep << " è " << realTotalEnergy << " GeV"<< endl;

    Etot_step.push_back(0.);
  }

  innerDepth = meanDepth[first_Ecal_step];
  if (last_Ecal_step + offset >= 0)
    outerDepth = meanDepth[last_Ecal_step + offset];
  else
    outerDepth = innerDepth;

  stepsCalculated = true;
}

void EMShower::compute()
{

  double t = 0.;
  double dt = 0.;
  if (!stepsCalculated)
    prepareSteps();
  // Prepare the grids in Ecal
  bool status = false;

  // cout << "Step preparati" << endl;

  // Loop over all segments for the longitudinal development
  double totECalc = 0;

  for (unsigned iStep = 0; iStep < nSteps; ++iStep)
  {
    // The length of the shower in this segment

    dt = steps[iStep].second;

    // cout << "Lunghezza dello step " << iStep << " è dt = " << dt <<endl;
    // The elapsed length
    t += dt;
    // Build the grid of crystals at this ECAL depth
    // Actually, it might be useful to check if this grid is empty or not.
    // If it is empty (because no crystal at this depth), it is of no use
    // (and time consuming) to generate the spots

    // middle of the step
    double tt = t - 0.5 * dt;
    double realTotalEnergy = 0.;
    for (int i = 0; i < nPart; ++i)
    {
      realTotalEnergy += depositedEnergy[iStep][i] * E[i];
    }
    // cout << "Allo step " << iStep << " in tt = " << tt << " (metà step) ho E = " << realTotalEnergy << endl;

    // If the amount of energy is greater than 1 MeV, make a new grid
    // otherwise put in the previous one.
    bool usePreviousGrid = (realTotalEnergy < 0.001);

    // If the amount of energy is greater than 1 MeV, make a new grid
    // otherwise put in the previous one.

    // If less than 1 kEV. Just skip

    if (iStep > 2 && realTotalEnergy < 0.000001)
      continue;

    if (!usePreviousGrid)
    {
      if (tt > 24.5)
        status = false;
      else
        status = true;
    }

    if (!status)
      continue;

    bool detailedShowerTail = false;
    // check if a detailed treatment of the rear leakage should be applied
    if (!usePreviousGrid)
    {
      // E' UNA PROVA!!!
      // detailedShowerTail = (t - dt > theGrid->getX0back());
      detailedShowerTail = (t - dt > 24.7 - X0depth);
    }

    // The particles of the shower are processed in parallel
    for (int i = 0; i < nPart; ++i)
    {

      // cout << "questa è l'energi della particella " << E[i] << " VS critical energy " << theECAL->criticalEnergy() << endl;

      //  integration of the shower profile between t-dt and t
      double dE = depositedEnergy[iStep][i];
      // cout << " % di enrgia depositata dalla particella è E%= " << dE << endl;
      // no need to do the full machinery if there is ~nothing to distribute)

      if (dE * E[i] < 0.000001)
        continue;

      totECalc += dE;

      // The number of energy spots (or mips)
      double nS = 0;

      // ECAL case : Account for photostatistics and long'al non-uniformity

      dE = gRandom->Poisson(dE * photos[i]) / photos[i];
      double z0 = gRandom->Gaus(0., 1.);
      dE *= 1. + z0 * theECAL->lightCollectionUniformity();

      // Expected spot number
      nS = (theNumberOfSpots[i] * gam(bSpot[i] * tt, aSpot[i]) * bSpot[i] * dt / tgamma(aSpot[i]));
      // +

      if (detailedShowerTail)
        myGammaGenerator->setParameters(floor(a[i] + 0.5), b[i], t - dt);

      // The lateral development parameters

      // Energy of the spots
      double eSpot = (nS > 0.) ? dE / nS : 0.;
      double SpotEnergy = eSpot * E[i];

      //   cout << "La cui energia SpotEnergy " << SpotEnergy << endl;

      int nSpot = (int)(nS + 0.5);

      // double taui = t/T[i];
      double taui = tt / Ti[i];
      double proba = theParam->p(taui, E[i]);
      double theRC = theParam->rC(taui, E[i]);
      double theRT = theParam->rT(taui, E[i]);

      double dSpotsCore = gRandom->Gaus(proba * nSpot, std::sqrt(proba * (1. - proba) * nSpot));

      if (dSpotsCore < 0)
        dSpotsCore = 0;

      unsigned nSpots_core = (unsigned)(dSpotsCore + 0.5);
      unsigned nSpots_tail = ((unsigned)nSpot > nSpots_core) ? nSpot - nSpots_core : 0;

      // cout << "il numero di spot nel core " << nSpots_core << endl;
      // cout << "il numero di spot nella tail " << nSpots_tail << endl;

      for (unsigned icomp = 0; icomp < 2; ++icomp)
      {
        double theR = (icomp == 0) ? theRC : theRT;
        unsigned ncompspots = (icomp == 0) ? nSpots_core : nSpots_tail;

        // cout << "ora siamo in icomp = " << icomp << " quindi il raggio è " << theR << endl;
        // cout << "e la spot Energy " << SpotEnergy << endl;

        RadialInterval radInterval(theR, ncompspots, SpotEnergy /*, random*/);
        if (icomp == 0)
        {
          setIntervals(icomp, radInterval);
        }
        else
        {
          setIntervals(icomp, radInterval);
        }

        radInterval.compute();
        // irad = 0 : central circle; irad=1 : outside

        unsigned nrad = radInterval.nIntervals();

        for (unsigned irad = 0; irad < nrad; ++irad)
        {
          double spote = radInterval.getSpotEnergy(irad);
          theGrid->setSpotEnergy(spote);

          // cout << "Siamo nell'intervallo interno " << irad << " in cui la spot energy è " << spote << endl;

          unsigned nradspots = radInterval.getNumberOfSpots(irad);

          // cout << "Dopo aver fatto passaggi vari (moltiplicato per percentual etc..) trovo che il vero N spot = " << nradspots << endl;

          double umin = radInterval.getUmin(irad);
          double umax = radInterval.getUmax(irad);
          // Go for the lateral development
          for (unsigned ispot = 0; ispot < nradspots; ++ispot)
          {
            double z3 = gRandom->Uniform(umin, umax); //!!!!!!!!!
            double ri = theR * std::sqrt(z3 / (1. - z3));

            // Generate phi
            double phi = 2. * M_PI * gRandom->Uniform(); //!!!!!!!!!
                                                         // cout << "lo spot " << ispot << " si trova in (r,phi) = (" << ri << ", " << phi << ")" << endl;

            // Now the *moliereRadius is done in EcalHitMaker

            if (detailedShowerTail)
            {
              //			   std::// cout << "About to call addHitDepth " << std::endl;
              double depth;
              do
              {
                depth = myGammaGenerator->shoot();
              } while (depth > t);

              theGrid->AddHitCooDepth(ri, phi, Xi, Yi, spote, depth, X0depth);
              Etot[i] += spote;
              Etot_step[iStep] += spote;
              //			   std::// cout << " Done " << std::endl;
            }
            else

            {
              // This gives the number of the cell in which the particle impact
              // numberPart= theGrid->GiveCentralCell(Xi,Yi,EcalGrid);
              // This gives the number of the cell where the spot is set
              numberSpot = theGrid->AddHitCoo(ri, phi, Xi, Yi, spote);
              // if (numberSpot!=0)

              Etot[i] += spote;

              Etot_step[iStep] += spote;
            }
          }
        }
      }
      //  cout << " energia totale particella " << i << " è " << Etot[i]  << endl;
    }

    // cout << "-------> fine step numero " <<  iStep << " con Etotal_step = " << Etot_step[iStep] << " e con con Etot_step = " << Etot_step[iStep] << " con Etot 1 = " << Etot[0] << " e con Etot 2 = " << Etot[1] << endl;
  }
  double Etotal = 0.;

  for (int i = 0; i < nPart; ++i)
  {
    //      myHistos->fill("h10",Etot[i]);
    Etotal += Etot[i];
  }

  // theGrid->Draw_ECAL(EcalGrid);
}

double EMShower::gam(double x, double a) const
{

  // A stupid gamma function
  return std::pow(x, a - 1.) * std::exp(-x);
}

double EMShower::deposit(double t, double a, double b, double dt)
{

  double b1 = b * (t - dt);
  double b2 = b * t;
  double result = 0.;
  double rb1 = (b1 != 0.) ? myIncompleteGamma.MyGamma(a, b1) : 0.; // g->Eval(b1): 0.;
  double rb2 = (b2 != 0.) ? myIncompleteGamma.MyGamma(a, b2) : 0.; // g->Eval(b2): 0.;
  result = (rb2 - rb1);
  return result;
}

void EMShower::setIntervals(unsigned icomp, RadialInterval &rad)
{
  const std::vector<double> &myValues((icomp) ? theParam->getTailIntervals() : theParam->getCoreIntervals());
  unsigned nvals = myValues.size() / 2;
  for (unsigned iv = 0; iv < nvals; ++iv)
  {
    rad.addInterval(myValues[2 * iv], myValues[2 * iv + 1]);
  }
}

double EMShower::deposit(double a, double b, double t)
{
  double b2 = b * t;
  double result = 0.;
  if (fabs(b2) < 1.e-9)
    b2 = 1.e-9;
  result = myIncompleteGamma.MyGamma(a, b2);
  //  std::// cout << " deposit t = " << t  << " "  << result <<std::endl;
  return result;
}