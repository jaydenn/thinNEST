
#ifndef TARGET_HH
#define TARGET_HH

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <iostream>
#include <sys/time.h>
#include <assert.h>
#include <string>
 
//#include "Spectra.hh"

class Target {
    public:
        
        int migType;
        int migOptimize;
        double maxE,wimpMass;

        //target properties
        double Mt;
        double Leff;

        int elementZ;
        std::string target; 
        int isoN;
        int Norbits;
        std::vector<std::string> orbits;
        std::vector<double> isoTable;
        std::vector<double> isoFracTable;
        std::vector<double> orbits_energy;
        
        //interpolation objects
        gsl_spline *Znl_spline[17];
        gsl_interp_accel *Znl_accel[17];
        gsl_spline *Znl_int_spline[17];
        gsl_interp_accel *Znl_int_accel[17];
        gsl_spline *invZnl_int_spline[17];
        gsl_interp_accel *invZnl_int_accel[17];

        double Znl_y_max[17];
        double Znl_x_max[17];
        double Znl_x_min[17];
        double maxCumulativeProb;
        
        //random number generator
        const gsl_rng_type * T;
        gsl_rng * r;
        struct timeval tv;
        
        //integration workspace
        gsl_integration_workspace * W;
        
        //methods
        double Z_nl(int orbit_i, double v, double Ee);
        double Z_nl_max(int orbit_i, double v);
        double invZ_nl_integrated(int orbit_i, double v, double Znl);
        double Z_nl_integrated(int orbit_i, double v, double max_Ee);
        double totalMigProb(double maxE, double ENR, double MT);
        std::vector<double> rand_migdalE(double ERnr);
        double rand_isotope(double ENR);
        double ERmax(double Ein, double EM);
        double ERmax(double Ein, double EM, double MT);
        double EMmax(double Ein);
        double EMmax(double Ein, double MT);
        double EMmaxER(double Ein, double ER);
        double EMmaxER(double Ein, double ER, double MT);
//        double migdalIntegrand(double ER, void *par);
//       double dRdEmigdalNeutron(double Edet, Spectra *NRspec);
//        void calcMigdalSpectrum(Spectra NRspec);
        
        int init();
        
        //constructor
        Target(int Z, int type, double E, int op){
            elementZ = Z;
            migType = type;
            maxE = E;
            migOptimize = op;
            maxCumulativeProb =1;
            
            //initialize random number generator
            gsl_rng_env_setup();
            gettimeofday(&tv,0);     //use time as seed
            unsigned long mySeed = tv.tv_sec + tv.tv_usec;
            T = gsl_rng_default; // Generator setup
            r = gsl_rng_alloc (T);
            gsl_rng_set(r, mySeed);
            
            //initialize migdal calc
            if(init()<0)
                assert(0);
        };
};

#endif /* TARGET_HH */

