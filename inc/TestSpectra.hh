/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   TestSpectra.hh
 * Author: brodsky3
 *
 * Created on December 11, 2017, 10:27 AM
 */

#ifndef TESTSPECTRA_HH
#define TESTSPECTRA_HH

#include <assert.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include "RandomGen.hh"

#define NEST_AVO \
  6.0221409e+23       // good to keep in sync w/ NEST.hh, can't define twice
#define ATOM_NUM 54.  // ibid.

#define RHO_NAUGHT 0.3  // local DM halo density in [GeV/cm^3]
#define V_EARTH \
  245.  // for LUX Run03; if you want Run04 use 230 km/s (arXiv:1705.03380)
#define V_WIMP 220.
#define V_ESCAPE 544.

#define NUMBINS_MAX 200

class TestSpectra {
 public:
  TestSpectra(){};  // private so that it cannot be manually called

  struct WIMP_spectrum_prep {
    double base[100] = {1.};
    double exponent[100] = {0.};
    double integral = 0.;
    double xMax = 0.;
    double divisor = 1.;
  };
  WIMP_spectrum_prep wimp_spectrum_prep;
  
  struct SPLINE_spectrum_prep {
    gsl_spline *spectrumSpline;
    gsl_interp_accel *accelSS;
    gsl_spline *specSplineICDF;
    gsl_interp_accel *accelICDF;
    string type;
    string subType;
    string filename;
    double xMax = 0;
    double xMin = 0;
    double yMax = 0;
    double totRate = 0;
    bool doMig = 0;
    double monoE = 0;
    double MFP = -1;
    double beamWidth = -1;
  };
  SPLINE_spectrum_prep spline_spectrum_prep;
  
  struct SPLINE_migdal_prep {
    gsl_spline *migdalSpline;
    gsl_interp_accel *accelMS;
    double xMax = 0;
    double xMin = 0;
    double yMax = 0;
    double totRate = 0;
    double ER=0;
    double monoE;
    int N=-1;
    SPLINE_spectrum_prep *NR_spec;
  };
  SPLINE_migdal_prep migdal_splines[6];
  bool doMigdal;
  bool isLshell=0;

  double CH3T_spectrum(double emin, double emax);
  double C14_spectrum(double emin, double emax);
  double B8_spectrum(double emin, double emax);
  double file_spectrum(SPLINE_spectrum_prep sSpec);
  double file_spectrum_invCDF(SPLINE_spectrum_prep sSpec);
  SPLINE_spectrum_prep SPLINE_read_spectrum_file(const string filename, SPLINE_spectrum_prep spectrum);
  double AmBe_spectrum(double emin, double emax);
  double Cf_spectrum(double emin, double emax);
  double DD_spectrum(double emin, double emax);
  double WIMP_dRate(double ER, double mWimp);
  WIMP_spectrum_prep WIMP_prep_spectrum(double mass, double eStep);
  double WIMP_spectrum(WIMP_spectrum_prep wprep, double mass);
  double ZeplinBackground();  // an example of how to do a better (non-flat) ER
                              // BG spectrum for a WS, from Henrique Araujo
  double radNC_spectrum();
};

#endif /* TESTSPECTRA_HH */
