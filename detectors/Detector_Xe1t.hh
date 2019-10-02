// Detector specifications for Xe1t (SR1)
//   created by JL Newstead
//   using info from 1906.04717 and 1902.11297

#ifndef Detector_Xe1tSH_hh
#define Detector_Xe1tSH_hh 1

#include "VDetector.hh"

using namespace std;

class Detector_Xe1t : public VDetector {
 public:
  Detector_Xe1t() {
    cout << "You are currently using the Xe1t detector."
         << endl;

    // Call the initialisation of all the parameters
    Initialization();
  };
  virtual ~Detector_Xe1t(){};

  // Do here the initialization of all the parameters that are not varying as a
  // function of time
  virtual void Initialization() {
    // Primary Scintillation (S1) parameters
    g1 = 0.1165;  // phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect (0.142/1.219) or .11813 for 21.9%
    sPEres = 0.28;   // single phe resolution (Gaussian assumed) (from 1906.00819)
    //sPEthr = 0.35;   // POD threshold in phe, usually used IN PLACE of sPEeff
    sPEeff = 0.93;   // actual efficiency, can be used in lieu of POD threshold
    noise[0] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    noise[1] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    P_dphe = 0.219;  // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT (from 1902.11297)

    coinWind = 100;  // S1 coincidence window in ns
    coinLevel = 3;   // how many PMTs have to fire for an S1 to count
    numPMTs = 248-36;    // For coincidence calculation 1906.04717

    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.0823;  // phd per S2 photon in gas, used to get SE size (set to match SE size of 28.2 PE/SE = 23.13 phd/SE from 1906.04717)
    s2Fano = 3.65;  // Fano-like fudge factor for SE width (set to match SE width of 7.2 PE = 5.9 phd from 1906.04717)
    s2_thr = 90.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
    E_gas = 11.045;    // field in kV/cm between liquid/gas border and anode (set to match ext_eff)
    eLife_us = 641.;  // the drift electron mean lifetime in micro-seconds 1906.04717

    // Thermodynamic Properties
    inGas = false;
    T_Kelvin = 273.-96.;  // for liquid drift speed calculation 1906.04717
    p_bar = 1.938;     // gas pressure in units of bars, it controls S2 size 1906.04717
    // if you are getting warnings about being in gas, lower T and/or raise p

    // Data Analysis Parameters and Geometry
    dtCntr = 238.;  // center of detector for S1 corrections, in usec.
    dt_min = 22.;  // minimum. Top of detector fiducial volume
    dt_max = 474.;  // maximum. Bottom of detector fiducial volume

    radius = 420.;  // millimeters (fiducial rad)
    radmax = 480.;  // actual physical geo. limit

    TopDrift = 970.;  // mm not cm or us (but, this *is* where dt=0)
    // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
    // PMTs)
    // In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not
    // anode!
    anode = 972.5;  // the level of the anode grid-wire plane in mm -Set to match g2 value
    // In a gas TPC, this is not TopDrift (top of drift region), but a few mm
    // above it
    gate = 967.5;  // mm. This is where the E-field changes (higher) -LUX value
    // in gas detectors, the gate is still the gate, but it's where S2 starts
    cathode = 5;  // mm. Defines point below which events are gamma-X

    // 2-D (X & Y) Position Reconstruction
    PosResExp = 0.015;     // exp increase in pos recon res at hi r, 1/mm
    PosResBase = 30.;  // baseline unc in mm, see NEST.cpp for usage
  }

  // S1 PDE custom fit for function of z
  // s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  virtual double FitS1(double xPos_mm, double yPos_mm, double zPos_mm) {
    return 1.;  // unitless, 1.000 at detector center
  }

  // Drift electric field as function of Z in mm
  // For example, use a high-order poly spline
  virtual double FitEF(double xPos_mm, double yPos_mm,
                       double zPos_mm) {  // in V/cm
    return 81.;
  }

  // S2 PDE custom fit for function of r
  // s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  virtual double FitS2(double xPos_mm, double yPos_mm) {
    return 1.;  // unitless, 1.000 at detector center
  }

  virtual vector<double> FitTBA(double xPos_mm, double yPos_mm,
                                double zPos_mm) {
    vector<double> BotTotRat(2);

    BotTotRat[0] = 0.6;  // S1 bottom-to-total ratio
    BotTotRat[1] = 0.4;  // S2 bottom-to-total ratio, typically only used for
                         // position recon (1-this)

    return BotTotRat;
  }

  virtual double OptTrans(double xPos_mm, double yPos_mm, double zPos_mm) {
    double phoTravT, approxCenter = (TopDrift + cathode) / 2.,
                     relativeZ = zPos_mm - approxCenter;

    double A = 0.048467 - 7.6386e-6 * relativeZ +
               1.2016e-6 * pow(relativeZ, 2.) - 6.0833e-9 * pow(relativeZ, 3.);
    if (A < 0.) A = 0.;  // cannot have negative probability
    double B_a = 0.99373 + 0.0010309 * relativeZ -
                 2.5788e-6 * pow(relativeZ, 2.) -
                 1.2000e-8 * pow(relativeZ, 3.);
    double B_b = 1. - B_a;
    double tau_a = 11.15;  // all times in nanoseconds
    double tau_b = 4.5093 + 0.03437 * relativeZ -
                   0.00018406 * pow(relativeZ, 2.) -
                   1.6383e-6 * pow(relativeZ, 3.);
    if (tau_b < 0.) tau_b = 0.;  // cannot have negative time

    // A = 0.0574; B_a = 1.062; tau_a = 11.1; tau_b = 2.70; B_b = 1.0 - B_a;
    // //LUX D-D conditions

    if (RandomGen::rndm()->rand_uniform() < A)
      phoTravT = 0.;  // direct travel time to PMTs (low)
    else {            // using P0(t) =
            // A*delta(t)+(1-A)*[(B_a/tau_a)e^(-t/tau_a)+(B_b/tau_b)e^(-t/tau_b)]
            // LUX PSD paper, but should apply to all detectors w/ diff #'s
      if (RandomGen::rndm()->rand_uniform() < B_a)
        phoTravT = -tau_a * log(RandomGen::rndm()->rand_uniform());
      else
        phoTravT = -tau_b * log(RandomGen::rndm()->rand_uniform());
    }

    double sig = RandomGen::rndm()->rand_gauss(
        3.84, .09);  // includes stat unc but not syst
    phoTravT += RandomGen::rndm()->rand_gauss(
        0.00, sig);  // the overall width added to photon time spectra by the
                     // effects in the electronics and the data reduction
                     // pipeline

    if (phoTravT > DBL_MAX) phoTravT = tau_a;
    if (phoTravT < -DBL_MAX) phoTravT = 0.000;

    return phoTravT;  // this function follows LUX (arXiv:1802.06162) not Xe10
                      // technically but tried to make general
  }

  virtual vector<double> SinglePEWaveForm(double area, double t0) {
    vector<double> PEperBin;

    double threshold = PULSEHEIGHT;  // photo-electrons
    double sigma = PULSE_WIDTH;      // ns
    area *= 10. * (1. + threshold);
    double amplitude = area / (sigma * sqrt(2. * M_PI)),
           signal;  // assumes perfect Gaussian

    double tStep1 = SAMPLE_SIZE / 1e2;  // ns, make sure much smaller than
                                        // sample size; used to generate MC-true
                                        // pulses essentially
    double tStep2 =
        SAMPLE_SIZE;  // ns; 1 over digitization rate, 100 MHz assumed here

    double time = -5. * sigma;
    bool digitizeMe = false;
    while (true) {
      signal = amplitude * exp(-pow(time, 2.) / (2. * sigma * sigma));
      if (signal < threshold) {
        if (digitizeMe)
          break;
        else
          ;  // do nothing - goes down to advancing time block
      } else {
        if (digitizeMe)
          PEperBin.push_back(signal);
        else {
          if (RandomGen::rndm()->rand_uniform() < 2. * (tStep1 / tStep2)) {
            PEperBin.push_back(time + t0);
            PEperBin.push_back(signal);
            digitizeMe = true;
          } else {
          }
        }
      }
      if (digitizeMe)
        time += tStep2;
      else
        time += tStep1;
      if (time > 5. * sigma) break;
    }

    return PEperBin;
  }
  // Vary VDetector parameters through custom functions
  virtual void ExampleFunction() { set_g1(0.0760); }
};

#endif
