#include <unistd.h>
#include <getopt.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits.h>

#include "NEST.hh"

#include "TestSpectra.hh"
#include "Target.hh"
#include "PhysicalConstants.hh"
#include "trigger.hh"

#include "analysis_default.hh"
#include "Detector_default.hh"

using namespace std;
using namespace NEST;
void setDetectorPars(string detName, Detector_def *det);
void setAnalysisPars(string analysisFilename);

int main(int argc, char** argv) 
{
    //create detector object
    Detector_def* detector = new Detector_def();
    
    //set some parameters that nest uses as input for it's yield calculations (these are default values from nest)
    vector<double> signalE, vTable, FreeParam,
    freeParamDef = {1,1,0.1,0.5,0.19,2.25},
    freeParamBeta = {0.5,0.5,1.1,-5.,1.01,0.95,1.4e-2,1.8e-2};
    
    // Fi, Fex, and 3 non-binomial recombination fluctuation parameters
    //NuisParam = {11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1., 1.};
    // alpha,beta,gamma,delta,epsilon,zeta,eta,theta,iota for NR model
    // last 3 are the secret extra parameters for additional flexibility - all parameter values current as of 03/2022

    //version 2.3.9 update
    vector<double> NRYieldsParam =  {11., 1.1, 0.0480, -0.0533, 12.6, 0.3, 2., 0.3, 2., 0.5, 1., 1.},
    NRERWidthsParam = {1.0,1.0, 0.1,0.5,0.19,2.25,0.0015,  0.0553,0.205, 0.45,-0.2},
    LZ_NRERWidthsParam =     {0.4,0.4,0.04,0.5,0.19,2.25,0.0015,0.046452,0.205, 0.45,-0.2},
    ERWeightParam;
    
    //define a bunch of variables
    string position, delimiter, token;
    size_t loc;
    double g2=0, atomNum = 0, massNum = 0;
    YieldResult yieldsMax;

    unsigned long int numEvents = 0;
    double exposure=-1;
    string type = "NR";
    string detectorName = "-1";
    string spectrumFilename = "-1";
    string analysisFilename = "-1";
    string outputFilename = "-1";
    double eMin = -1;
    double eMax = -1;   
    double inField = -1;
    position = "-1";
    double fPos = -1;
    int migdal=0; int mig_type = 0; int migdalOptimize = 0;
    int seed=-1; 
    int verbose=0;
    int progress=0;
    int outputQuanta=0;
    int outputLindhard=0;
    int doBinning=0;
    int usePosition=0;
    double wimpMass=1;
    double wimpCS=1e-35;
    int ratioEvents = 0;
    
    //set up the options for the code
    const struct option longopts[] =
    {
        {"help",         no_argument,        0, 'h'},
        {"migdal",       no_argument,        0, 'm'},
        {"migdal (opt.)",no_argument,        0, 'M'},
        {"binned",       no_argument,        0, 'b'},
        {"timing",       no_argument,        0, 't'},
        {"verbose",      no_argument,        0, 'v'},
        {"outputQuanta", no_argument,        0, 'q'},
        {"outLindhard",  no_argument,        0, 'l'},
        {"position",     no_argument,        0, 'p'},
        {"progress",     no_argument,        0, 'r'},
        {"numEvents",    required_argument,  0, 'n'},
        {"exposure",     required_argument,  0, 'N'},
        {"spectra",      required_argument,  0, 's'},
        {"detector",     required_argument,  0, 'd'},
        {"analysis",     required_argument,  0, 'a'},
        {"eventPosition",required_argument,  0, 'P'},
        {"eMin",         required_argument,  0, 'e'},
        {"eMax",         required_argument,  0, 'E'},
        {"field",        required_argument,  0, 'f'},
        {"seed",         required_argument,  0, 'S'},
        {"output",       required_argument,  0, 'o'},
        {"ratioEvents",  required_argument,  0, 'R'},
        
        {0,0,0,0},
    };
    int iarg=0;
    opterr=1;     //turn off getopt error message
    int longindex;
    while(iarg != -1)
    {
        iarg = getopt_long(argc, argv, "hmbtvpqlrROMn:N:s:f:d:P:e:E:f:S:a:o:", longopts, &longindex);

        switch (iarg)
        {
            case 'm':
                migdal = 1;
                break;
            case 'M':
            {
                migdal = 1;
                migdalOptimize=1;
                break;
            }
            case 't':
                useTiming = 2;
                break;
            case 'b':
                doBinning = 1;
                break;        
            case 'v':
                verbose = 1;
                break;
            case 'q':
                outputQuanta = 1;
                break;
            case 'l':
                outputLindhard = 1;
                break;
            case 'p':
                usePosition = 1;
                break;
            case 'n':
                numEvents = atoi(optarg);
                break;
            case 'N':
                exposure = atof(optarg);
                break;
            case 'r':
                progress = 1;
                break;
            case 'R':
                ratioEvents = 1;
                break;
            case 's':
            {
                type = optarg;
                if (type.substr(0,4) == "file")
                {
                    spectrumFilename = type.substr(5,string::npos);
                    type = "file";
                }
                else if (type.substr(0,4) == "WIMP")
                {
                    wimpMass = stof(type.substr(5,type.find(",")-5));
                    wimpCS = stod(type.substr(type.find(",")+1,10));
                    type = "WIMP";
                }
                break;
            }            
            case 'a':
                analysisFilename = optarg;
                break;
            case 'd':
                detectorName = optarg;
                break;
            case 'P':
                position = optarg;
                fPos=1;
                break;
            case 'e':
                eMin = atof(optarg);
                break;
            case 'E':
                eMax = atof(optarg);
                break;
            case 'f':
                inField = atof(optarg);
                break;  
            case 'S':
                seed = atoi(optarg);
                break; 
            case 'o':
                outputFilename = optarg;
                break;
            case '?':
            case 'h':
                cout << "Usage: thinNEST [options]\n"
                     << "A lightweight and flexible code for running NEST simulations\n\n"
                     << "Options:\n"
                     << "\t-v, --verbose\tverbose output\n"
                     << "\t-t, --timing\tinclude s2-width calculation\n"
                     << "\t-m, --migdal\tinclude the migdal effect for NR\n"
                     << "\t-M, --migdal\tinclude the migdal effect for NR (optimized, relative NR/Mig rate will be unphysical)\n"
                     << "\t-n, --numEvents N\tsimulate N events\n"
                     << "\t-b, --binned \toutputs events binned in S1 and logS2 as specified in analysis file\n"
                     << "\t-q, --outputQuanta\tincludes quanta for each event in output\n"
                     << "\t-P, --eventPosition\tincludes position of each event in output\n"
                     << "\t-a, --analysis analysis_file\n"
                     << "\t-d, --detector detector_file\n"
                     << "\t-s, --spectrum {NR,ER,file=spectrum_file}\n"
                     << "\t-e, --eMin X\tstarts the spectrum at X keV, defaults to -1 (min in spectrum file)\n"
                     << "\t-E, --eMax X\tends the spectrum at X keV, defaults to -1 (max in spectrum file)\n"
                     << "\t-o, --output output_file, defaults to stdout\n"
                     << "\t-f, --field V\tset the drift field to V volts/cm\n"
                     << "\t-S, --seed S\tuse the random seed S, defaults to sys clock\n"               
                     << "\t-r, --progress \tgives percent done (out of number of samples) to stdout\n";               
                return 0;
        }

    }
    
    //set random seed or take it from the clock (default)
    if (seed == -1)
        RandomGen::rndm()->SetSeed(time(NULL));
    else
        RandomGen::rndm()->SetSeed(seed);

    //open file for output or use cout (default)
    ofstream outputFile; 
    if(outputFilename!="-1")
        outputFile.open (outputFilename);
    ostream & outStream = (outputFilename!="-1" ? outputFile : cout);

    // detector parameter modifications
    setDetectorPars(detectorName, detector);
    outStream << "using the " << (detectorName!="-1" ? detectorName : "default" ) << " detector file\n";

    if(detectorName.find("LZ_SR1") != std::string::npos)
    {    
        NRERWidthsParam = LZ_NRERWidthsParam;
        cout << "USING LZ WIDTH PARAMETERS\n";
    }
    
    // analysis parameter modifications
    setAnalysisPars(analysisFilename);
    outStream << "using the " << (analysisFilename!="-1" ? analysisFilename : "default" ) << " analysis file" << std::endl;
    if (verbose==1) verbosity=true; //overwrite analysis with command line arg

    //set up array for binned data storage
    int**** s1s2RZbins;
    if(usePosition!=1 && doBinning==1)
        numBinsZ = numBinsR = 1;
    s1s2RZbins = new int***[numBinsS1];
    for(int i = 0; i < numBinsS1; ++i)
    {
        s1s2RZbins[i] = new int**[numBinsS2];
        for(int j = 0; j < numBinsS2; ++j)
        {
            s1s2RZbins[i][j] = new int*[numBinsR];
            for(int k = 0; k < numBinsR; ++k)
            {
                 s1s2RZbins[i][j][k] = new int[numBinsZ]();
            }
        }
    }
    
    //calculate bin widths
    double s1binWidth = (maxS1-minS1)/numBinsS1;
    if (doFiducialCut)
        maxR*=detector->get_radius(); 
    else
        maxR*=detector->get_radmax();
    minR*=detector->get_radius();
    double RbinWidth = (maxR-minR)/numBinsR;
    double s2binWidth;
    if(logS2==1)
        s2binWidth = (log10(maxS2)-log10(minS2))/numBinsS2;
    else
        s2binWidth = (maxS2-minS2)/numBinsS2;

    // Construct NEST class using detector object
    NESTcalc n(detector);
    
    INTERACTION_TYPE type_num;
    TestSpectra spec;

    //do some basic checks of detector geometry
    if (detector->get_TopDrift() <= 0. || detector->get_anode() <= 0. ||
        detector->get_gate() <= 0.) 
    {
        cerr << "ERROR, unphysical value(s) of position within the detector "
            "geometry.";  // negative or 0 for cathode position is OK (e.g., LZ)
        return 0;
    }
    
    //print out type of spectrum being simulated
    outStream << "spectrum type: " << (type=="file" ? "file ("+spectrumFilename+")" : type ) << endl;
    
    //set up the spectrum when using a precomputed recoil spectrum from a file
    if (type == "file")
    {
    
        spec.spline_spectrum_prep.xMin = eMin;
        spec.spline_spectrum_prep.xMax = eMax;
        spec.spline_spectrum_prep = spec.SPLINE_read_spectrum_file(spectrumFilename,spec.spline_spectrum_prep);
        //eMin and eMax may be adjusted to account for the endpoints in the file        
        eMin = spec.spline_spectrum_prep.xMin;
        eMax = spec.spline_spectrum_prep.xMax;

        //get the type of spectrum from first line
        if (spec.spline_spectrum_prep.type == "NR" || spec.spline_spectrum_prep.type == "neutronB")
            type_num = NR;
        else if (spec.spline_spectrum_prep.type == "ER")
            type_num = beta;
        else if (spec.spline_spectrum_prep.type == "ERl")
        {
            spec.isLshell=1;
            type_num = beta;
        }
        else
        {
            cerr << "type " << spec.spline_spectrum_prep.type << " not recognized\n";
            assert(0);
        }
        
        //extra options of the type of spectrum which can also affect calculation of the Migdal effect
        if (spec.spline_spectrum_prep.subType == "monoE")
            mig_type = NEUTRON;
        else if (spec.spline_spectrum_prep.subType == "neutrino")
            mig_type = NEUTRINO;
        else if (spec.spline_spectrum_prep.subType == "wimpNR")
            mig_type = WIMPmig;

        if ( numEvents == 0 )
        {
            double expectedN = spec.spline_spectrum_prep.totRate * exposure;
            if (expectedN > 2.e9) //I thought it was supposed to be 4.29e9, but was getting problems at lower numbers
            {
                cout << "exposure expects more events than can be stored in an unsigned long int, reduce exposure\n";
                return 0;
            }            
            numEvents = RandomGen::rndm()->poisson_draw(expectedN);
            cout << "simulating " << numEvents << " events (" << spec.spline_spectrum_prep.totRate * exposure << " expected from an exposure of " << exposure << ")\n";
        }
        else
            cout << "simulating " << numEvents << " events\n";

    }    
    else if (type == "NR" || type == "neutron" || type == "-1" || type == "neutronM") //-1: default particle type is also NR
       type_num = NR;  
    else if (type == "WIMP")
    {
        if (wimpMass < 0.44) 
        {
            cerr << "WIMP mass too low, you're crazy!" << endl;
            return 0;
        }
    
        type_num = WIMP;
        spec.wimp_spectrum_prep = spec.WIMP_prep_spectrum(wimpMass, E_step);
        numEvents =
            RandomGen::rndm()->poisson_draw(spec.wimp_spectrum_prep.integral * 1.0 *
                                            exposure * wimpCS / 1e-36);
        cout << "simulating " << numEvents << " events\n";
    }
    else if (type == "B8" || type == "Boron8" || type == "8Boron" ||
             type == "8B" || type == "Boron-8") 
    {
        type_num = B8;
        numEvents = RandomGen::rndm()->poisson_draw(0.0026 * stof(argv[1]));
    }
    else if (type == "DD" || type == "D-D")
        type_num = DD;
    else if (type == "AmBe")
        type_num = AmBe;
    else if (type == "Cf" || type == "Cf252" || type == "252Cf" ||
           type == "Cf-252")
        type_num = Cf;
    else if (type == "ion" || type == "nucleus" || type == "alpha") 
    {
        type_num = ion;
        if (type == "alpha") 
        {
            atomNum = 2;
            massNum = 4;
        } 
        else
        {
            cerr << "Atomic Number: ";
            cin >> atomNum;
            cerr << "Mass Number: ";
            cin >> massNum;
        }
        if (atomNum == ATOM_NUM) 
            type_num = NR;
    }
    else if (type == "gamma" || type == "gammaRay" || type == "x-ray" ||
             type == "xray" || type == "xRay" || type == "X-ray" ||
             type == "Xray" || type == "XRay")
        type_num = gammaRay;  // includes photo-absorption and electron capture
    else if (type == "Kr83m" || type == "83mKr" || type == "Kr83")
        type_num = Kr83m;
    else if (type == "CH3T" || type == "tritium")
        type_num = CH3T;
    else if (type == "C14" || type == "Carbon14" || type == "14C")
        type_num = C14;
    else if (type == "beta" || type == "ER" || type == "Compton" ||
           type == "compton" || type == "electron" || type == "e-" ||
           type == "muon" || type == "MIP" || type == "LIP" || type == "mu" ||
           type == "mu-")
        type_num = beta;  // default electron recoil model
    else if (type == "RNC")
    {
        type_num = gammaRay;
    }
    else 
    {
        cerr << "UNRECOGNIZED PARTICLE TYPE!! VALID OPTIONS ARE:" << endl;
        cerr << "NR or neutron," << endl;
        cerr << "WIMP," << endl;
        cerr << "B8 or Boron8 or 8Boron or 8B or Boron-8," << endl;
        cerr << "DD or D-D," << endl;
        cerr << "AmBe," << endl;
        cerr << "Cf or Cf252 or 252Cf or Cf-252," << endl;
        cerr << "ion or nucleus," << endl;
        cerr << "alpha," << endl;
        cerr << "gamma or gammaRay," << endl;
        cerr << "x-ray or xray or xRay or X-ray or Xray or XRay," << endl;
        cerr << "Kr83m or 83mKr or Kr83," << endl;
        cerr << "CH3T or tritium," << endl;
        cerr << "Carbon14 or 14C or C14," << endl;
        cerr << "beta or ER or Compton or compton or electron or e-, and" << endl;
        cerr << "muon or MIP or LIP or mu or mu-" << endl;
        return 0;
    }

    //set NEST parameters based on recoil type (ER or NR)
    if (type_num == beta)
        FreeParam = freeParamBeta;
    else
        FreeParam = freeParamDef;

    //checks on energy range of recoils
    if (eMin == -1.)
        eMin = 0.;
    if (eMax == -1. && eMin == 0.)
        eMax = 1e4;  // the default energy max is 10 MeV
    if (eMax == 0.) 
    {
        cerr << "ERROR: The maximum energy cannot be 0 keV!" << endl;
        return 0;
    }

    //if including the Migdal effect do some setup
    Target myTarget( 54, mig_type, spec.spline_spectrum_prep.monoE, migdalOptimize);
    if (migdal==1)
    {
        spec.doMigdal = 1;
        cout << "initializing migdal calc.." << endl;
        //init_Znl(spec.spline_spectrum_prep.monoE, mig_type, migdalOptimize);
        //optional output theoretical spectrum [currently turned off]
        //calcMigdalSpectrum(&(spec.spline_spectrum_prep));
    }
    else
        spec.doMigdal = 0;

    if (type_num == Kr83m)
    {
        if (eMin == 9.4 && eMax == 9.4) 
        {
        }
        else if (eMin == 32.1 && eMax == 32.1) 
        {
        }
        else 
        {
            cerr << "ERROR: For Kr83m, put both energies as 9.4 or both as 32.1 keV "
              "please."
            << endl;
            return 0;
        }
    }

    if ((eMin < 10. || eMax < 10.) && type_num == gammaRay) 
    {
        cerr << "WARNING: Typically beta model works better for ER BG at low "
            "energies as in a WS."
         << endl;
        cerr << "ER data is often best matched by a weighted average of the beta & "
            "gamma models."
         << endl;
    }

    //get density of the xenon based on temp and pressure
    double rho = n.SetDensity(detector->get_T_Kelvin(), detector->get_p_bar());
    if (rho <= 0. || detector->get_T_Kelvin() <= 0. ||
        detector->get_p_bar() <= 0.) 
    {
        cerr << "ERR: Unphysical thermodynamic property!";
        return 0;
    }
    if (rho < 1.75) 
        detector->set_inGas(true);
    
    //get work function based on density
    double Wq_eV = NESTcalc::WorkFunction(rho,detector->get_molarMass()).Wq_eV;  // out-of-sync danger: copied from NEST.cpp
    if(verbose==1)
	cout << "Work function: " << Wq_eV << " eV\n"; 
    // Calculate and print g1, g2 parameters (once per detector)
    vector<double> g2_params = n.CalculateG2(verbosity);
    g2 = fabs(g2_params[3]);
    double g1 = detector->get_g1();
    double centralZ =
      (detector->get_gate() - 130. + detector->get_cathode() + 20) /
      2.;  // fid vol def usually shave more off the top, because of gas
           // interactions (100->10cm)
    double centralField = detector->FitEF(0.0, 0.0, centralZ);
    massNum = detector->get_molarMass();
    
    if (type_num == WIMP)
    {
        yieldsMax = n.GetYields(NR, 25.0, rho, centralField, double(massNum),
                            double(atomNum), NRYieldsParam);
    }
    else if (type_num == B8)
    {
        yieldsMax = n.GetYields(NR, 4.00, rho, centralField, double(massNum),
                            double(atomNum), NRYieldsParam);
    }
    else
    {
        double energyMaximum;
        if (eMax < 0.)
            energyMaximum = 1. / fabs(eMax);
        else
            energyMaximum = eMax;
        if (type_num == Kr83m)
        {            
            yieldsMax = n.GetYields(beta, energyMaximum, rho, centralField,
                              double(massNum), double(atomNum),
                              NRYieldsParam);  // the reason for this: don't do the
                                           // special Kr stuff when just
                                           // checking max
        }
        else
        {
            yieldsMax = n.GetYields(type_num, energyMaximum, rho, centralField,
                              double(massNum), double(atomNum), NRYieldsParam);
        }
    }
    if ((g1 * yieldsMax.PhotonYield) > (2. * maxS1) && eMin != eMax)
    {
        cerr << "\nWARNING: Your energy maximum may be too high given your maxS1.\n";
    }
    
    double vD_middle;
    //for now dont support variable field
    if(inField == -1)
        inField = detector->FitEF(0,0,0);
    
    vD_middle = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, inField);
    //if (inField == -1.)
    //{
        // build a vD table for non-uniform field, but if field varies in XY not
        // just Z you need to do more coding
    //    vTable = n.SetDriftVelocity_NonUniform(rho, z_step, 0, 0);
    //    vD_middle = vTable[int(floor(centralZ / z_step + 0.5))];
        // for ( int jj = 0; jj < vTable.size(); jj++ ) //DEBUG
        // cerr << double(jj)*z_step << "\t" << vTable[jj] << endl;
    //}
    //else
    

    //setup fiducial volume depths now that we have velocity - ASSUMING CONSTANT FIELD
    if(doFiducialCut)
    {
        maxZ=detector->get_TopDrift()-detector->get_dt_min()*vD_middle;
        minZ=detector->get_TopDrift()-detector->get_dt_max()*vD_middle;
    }
    else
    {
        maxZ=detector->get_TopDrift();
        minZ=0;
    }
    //set bin width in z direction
    double ZbinWidth = (maxZ-minZ)/numBinsZ; 


    if(verbose==1)
        cout << "drift velocity in middle: " << vD_middle << " mm/us";
    cout << endl;
    double keV = -999.;
    vector<double> migdalE = {0,0,0,massNum};

    //this code generates header for output columns
    string corr="";
    string unit="[phd]";
    string header="E[keV]\t\t";
    int outputPars;
    if(useCorrected==1)
        corr="c";
    if(usePD==0)
        unit="[phe]";
    if(outputQuanta == true)
    {
        outputPars=8;
        header.append("Nph\tNe-\tNhits\tNpe\tNeExt\t"+corr+"S1"+unit+"\t"+corr+"S2"+unit+"\t");
    }
    else
    {
        outputPars=3;
        header.append(corr+"S1"+unit+"\t"+corr+"S2"+unit+"\t");
    }
    if(spec.doMigdal == 1 && verbosity == true)
    {
        outputPars+=2;
        header.append("migE[keV]\t\tmigA[keV]\t\tmigN\t\t");
    }
    if(useTiming==2)
    {
        outputPars+=1;
        header.append("t80[ns]\t\t");
    }
    if(usePosition==1)
    {
        outputPars+=2;
        header.append("x[mm]\ty[mm]\tz[mm]\t\t");
    }
    if(MCtruthE == false && verbosity == true)
    {
        outputPars+=1;
        header.append("Etrue[keV]\t");
    }
    if(outputLindhard == 1)
    {
        outputPars+=1;
        header.append("lindhard");
    }

    if(doBinning==1)
        outStream <<  "Binned events, with bin spec: " << minS1 << " < "+corr+"S1"+unit+" < " << maxS1 <<", "<< minS2 << (logS2==1 ? " < log":" < ")+corr+"S2"+unit+" < " << maxS2 << "\n";
    else 
        outStream << header << "\n";

    // *** Main code loop for MC begins here ***
    unsigned long int numTrials=0;
        for ( unsigned long int j = 0; j < numEvents; j++) 
        {
            double signal1=0, signal2=0, smearRad=0,pos_x=0, pos_y=0, pos_z=0, r=0, phi=0, driftTime=0, field=0, vD=0;
            int index=0,indexR=0,indexZ=0,indexS1=0,indexS2=0;
            numTrials++; //keep track for when calculating effective exposure based off fixed event number simulation
            if (eMin == eMax && eMin >= 0. && eMax > 0.) 
                keV = eMin;
            else 
            {
                //this returns a random recoil energy from the given spectra
                if ( type == "file" )
                    keV = spec.file_spectrum_invCDF(spec.spline_spectrum_prep); //keV = spec.file_spectrum(spec.spline_spectrum_prep);
                else
                {
                    switch (type_num)
                    {
                        case CH3T:
                            keV = spec.CH3T_spectrum(eMin, eMax);
                            break;
                        case C14:
                            keV = spec.C14_spectrum(eMin, eMax);
                            break;
                        case B8:  // normalize this to ~3500 / 10-ton / year, for E-threshold of
                                  // 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
                            keV = spec.B8_spectrum(eMin, eMax);
                            break;
                        case AmBe:  // for ZEPLIN-III FSR from HA (Pal '98)
                            keV = spec.AmBe_spectrum(eMin, eMax);
                            break;
                        case Cf:
                            keV = spec.Cf_spectrum(eMin, eMax);
                            break;
                        case DD:
                            keV = spec.DD_spectrum(eMin, eMax);
                            break;
                        case WIMP: 
                            keV = spec.WIMP_spectrum(spec.wimp_spectrum_prep, wimpMass);
                            break;
                        default:
                            if (eMin < 0.) 
                                assert(0);
                            if (eMax > 0.)
                                keV = eMin + (eMax - eMin) * RandomGen::rndm()->rand_uniform();
                            else
                            {
                                // negative eMax signals to NEST that you want to use an
                                // exponential energy spectrum profile
                                if (eMin == 0.) 
                                    assert(0);
                                keV = 1e100;  // eMin will be used in place of eMax as the maximum
                                              // energy in exponential scenario
                                while (keV > eMin)
                                    keV = eMax * log(RandomGen::rndm()->rand_uniform());
                            }
                            break;
                    }
                    if(type == "RNC")
                        keV = spec.radNC_spectrum();
                }
            }
            if (migdal == 1 && type_num == NR)
                migdalE = myTarget.rand_migdalE(keV);  // Returns tuple of [electron energy, binding energy of shell, shell index, mass number of recoiling xe atom]
        //EDITED for testing
        //if (migdalE[0]>0)
        //    outStream << keV << "  " << migdalE[0] << "  " << migdalE[1] << "  " << migdalE[2] << endl;
        //continue;
        //cout << "testing Migdal energy, do not use for real results\n";
        //migdalE[0]+=migdalE[1];
        //migdalE[1]=0;
        
            if (type_num != WIMP && type_num != B8 && eMax > 0.) 
            {
                if (keV > eMax) 
                    keV = eMax;
                if (keV < eMin) 
                    keV = eMin;
            }

            Z_NEW:
            if (fPos == -1.)  // -1 means default, random location mode
            {  
                pos_z = 0. +
                      (detector->get_TopDrift() - 0.) *
                          RandomGen::rndm()->rand_uniform();  // initial guess
                r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
                phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
                pos_x = r * cos(phi);
                pos_y = r * sin(phi);
            }
            else
            {
                delimiter = ",";
                loc = 0;
                int i = 0;
                while ((loc = position.find(delimiter)) != string::npos) 
                {
                    token = position.substr(0, loc);
                    if (i == 0)
                        pos_x = stof(token);
                    else
                        pos_y = stof(token);
                    position.erase(0, loc + delimiter.length());
                    i++;
                }
                pos_z = stof(position);
                if (stof(position) == -1.) //choose a random position
                    pos_z =
                        0. +
                        (detector->get_TopDrift() - 0.) * RandomGen::rndm()->rand_uniform();
                if (stof(token) == -999.) 
                {
                    r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
                    phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
                    pos_x = r * cos(phi);
                    pos_y = r * sin(phi);
                }
                // if ( j == 0 ) { origX = pos_x; origY = pos_y; }
            }
            if (spec.spline_spectrum_prep.MFP!=-1) //if this is a beam with small MFP then place collision along x axis, MFP is calculated at density of 3g/cm3 so we correct for that here
            {
                r=1.e10;
                while(r>detector->get_radmax() || pos_z > detector->get_TopDrift() || pos_z < 0)
                {
                    //gaussian beam profile:
                    pos_y = RandomGen::rndm()->rand_gauss(0,spec.spline_spectrum_prep.beamWidth);
                    pos_z = RandomGen::rndm()->rand_gauss(0.5*detector->get_TopDrift(),spec.spline_spectrum_prep.beamWidth);
                    pos_x = spec.spline_spectrum_prep.MFP*rho/3*log(1/(1-RandomGen::rndm()->rand_uniform()))-sqrt(pow(detector->get_radmax(),2)+pow(pos_y,2));
                    r = sqrt(pow(pos_x,2)+pow(pos_y,2));
                }
            }

            if (inField == -1.) 
            {  // -1 means use poly position dependence
                field = detector->FitEF(pos_x, pos_y, pos_z);
            }
            else
                field = inField;  // no fringing

            if (field < 0. || detector->get_E_gas() < 0.) 
            {
                cerr << "\nERROR: Neg field is not permitted. We don't simulate field "
                      "dir (yet). Put in magnitude.\n";
                assert(0);
            }
            if (field == 0. || std::isnan(field))
                cerr << "\nWARNING: A LITERAL ZERO (or undefined) FIELD MAY YIELD WEIRD "
                      "RESULTS. USE A SMALL VALUE INSTEAD.\n";
            if (field > 12e3 || detector->get_E_gas() > 17e3)
                cerr << "\nWARNING: Your field is >12,000 V/cm. No data out here. Are "
                      "you sure about this?\n";
            
            if (inField == -1.) 
            {
                index = int(floor(pos_z / z_step + 0.5));
                vD = vTable[index];
            }
            else
                 vD = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, field);
            driftTime =
                (detector->get_TopDrift() - pos_z) / vD;  // (mm - mm) / (mm / us) = us
            if (inField != -1. &&
                detector->get_dt_min() > (detector->get_TopDrift() - 0.) / vD &&
                field >= FIELD_MIN) 
            {
                cerr << "ERROR: dt_min is too restrictive (too large)" << endl;
                assert(0);
            }
            if ((driftTime > detector->get_dt_max() ||
                 driftTime < detector->get_dt_min()) &&
                (fPos == -1. || stof(position) == -1.) && field >= FIELD_MIN && doFiducialCut)
                goto Z_NEW;
            if (detector->get_dt_max() > (detector->get_TopDrift() - 0.) / vD && !j &&
                field >= FIELD_MIN) 
            {
                cerr << "WARNING: dt_max is greater than max possible" << endl;
            }
            // The following should never happen: this is simply a just-in-case
            // code-block dealing with user error
            if (pos_z <= 0.) 
            {
                cerr << "ERROR: unphysically low Z coordinate (vertical axis of "
                          "detector) of "
                       << pos_z << " mm" << endl;
                assert(0);
            }
            if ((pos_z > (detector->get_TopDrift() + z_step) || driftTime < 0.0) &&
                field >= FIELD_MIN) 
            {
                cerr << "ERROR: unphysically big Z coordinate (vertical axis of "
                      "detector) of "
                   << pos_z << " mm" << endl;
                assert(0);
            }

            YieldResult yields;
            QuantaResult quanta;

                if (keV > .001 * Wq_eV || migdalE[0] > 0) 
                {
                    yields = n.GetYields(type_num, keV, rho, field, double(migdalE[3]),
                                     double(atomNum), NRYieldsParam);

                    if (migdal == 1 && type_num == NR && migdalE[0] > 0)
                    { 
                        YieldResult yieldsMigE;   //migdal: ejected electron
                        YieldResult yieldsMigDex; //migdal: de-excite atom
                
                        yieldsMigE = n.GetYields(beta, migdalE[0], rho, field, double(migdalE[3]),
                                     double(atomNum), NRYieldsParam);
                        yieldsMigDex = n.GetYields(beta, migdalE[1], rho, field, double(migdalE[3]),
                                     double(atomNum), NRYieldsParam);
                        
                        //As of version 2.3 this doesn't appear to be necessary
                        //if (int(migdalE[2]) == 2) //include this line for altering quanta from L-shell as observed in arXiv:2109.11487 
                        //{
                        //     double q = 0.91;      //this parameter may be weakly field dependent, this is for V=258V/cm
                        //     yieldsMigDex.PhotonYield += (1-q)*yieldsMigDex.ElectronYield;            
                        //     yieldsMigDex.ElectronYield *= q;     
                        //}
                        
                        yields.PhotonYield += yieldsMigE.PhotonYield + yieldsMigDex.PhotonYield;
                        yields.ElectronYield += yieldsMigE.ElectronYield + yieldsMigDex.ElectronYield;
                    }
                    quanta = n.GetQuanta(yields, rho, NRERWidthsParam, false, -999.);
                    if (spec.isLshell==1)                 
                    {    
                        cout << "not yet implemented\n";
                        return 0;
                    }
                }
                else 
                {
                    yields.PhotonYield = 0.;
                    yields.ElectronYield = 0.;
                    yields.ExcitonRatio = 0.;
                    yields.Lindhard = 0.;
                    yields.ElectricField = 0.;
                    yields.DeltaT_Scint = 0.;
                    quanta.photons = 0;
                    quanta.electrons = 0;
                    quanta.ions = 0;
                    quanta.excitons = 0;
                }

            // If we want the smeared positions (non-MC truth), then implement
            // resolution function
            double truthPos[3] = {pos_x, pos_y, pos_z};
            double smearPos[3] = {pos_x, pos_y, pos_z};
        
            smearRad = sqrt(pow(smearPos[0],2) + pow(smearPos[1],2));
            double Nphd_S2 =
                g2 * quanta.electrons * exp(-driftTime / detector->get_eLife_us());
            if (!MCtruthPos && Nphd_S2 > PHE_MIN) 
            {
                vector<double> xySmeared(2);
                xySmeared   = n.xyResolution(pos_x, pos_y, Nphd_S2);
                smearPos[0] = xySmeared[0];
                smearPos[1] = xySmeared[1];
                smearRad = sqrt(pow(smearPos[0],2) + pow(smearPos[1],2));
            }

            vector<long int> wf_time;
            vector<double> wf_amp;
            vector<double> scint =
                n.GetS1(quanta, truthPos[0], truthPos[1], truthPos[2], smearPos[0],smearPos[1],smearPos[2], vD, vD_middle, type_num, j, field,
                        keV, NEST::S1CalculationMode::Full, verbosity, wf_time, wf_amp);
            if (truthPos[2] < detector->get_cathode()) 
                quanta.electrons = 0;
            vector<double> scint2 = 
                n.GetS2(quanta.electrons, truthPos[0], truthPos[1], truthPos[2], smearPos[0],smearPos[1],smearPos[2], driftTime, vD, j, field,
                        NEST::S2CalculationMode::Full, verbosity, wf_time, wf_amp, g2_params);

            if (usePD == 0 && fabs(scint[2+useCorrected]) > minS1 && scint[2+useCorrected] < maxS1)
                signal1=scint[2+useCorrected];
            else if (usePD == 1 && fabs(scint[4+useCorrected]) > minS1 && scint[4+useCorrected] < maxS1)
                signal1=scint[4+useCorrected];
            else if (usePD >= 2 && fabs(scint[6+useCorrected]) > minS1 && scint[6+useCorrected] < maxS1)
                signal1=scint[6+useCorrected];
            else
                signal1=-999.;
            
            if (usePD == 0 && fabs(scint2[4+useCorrected]) > minS2 && scint2[4+useCorrected] < maxS2)
                signal2=scint2[4+useCorrected];
            else if (usePD >= 1 && fabs(scint2[6+useCorrected]) > minS2 && scint2[6+useCorrected] < maxS2)
                signal2=scint2[6+useCorrected];  // no spike option for S2
            else
                signal2=-999.;

            if(signal1 > 0 && signal2 > 0 ) //&&  migdalE[0] > 0
            {
                //inside fiducial vol?
                if( smearRad<maxR && smearPos[2]<maxZ && smearPos[2] > minZ)
                {
                    double keVtrue = keV;
                    if (!MCtruthE)
                    {
                        double Nph, Ne;
                        if (usePD == 0)
                            Nph = fabs(scint[3]) / (g1 * (1. + detector->get_P_dphe()));
                        else if (usePD == 1)
                            Nph = fabs(scint[5]) / g1;
                        else
                            Nph = fabs(scint[7]) / g1;
                        if (usePD == 0)
                            Ne = fabs(scint2[5]) / (g2 * (1. + detector->get_P_dphe()));
                        else
                            Ne = fabs(scint2[7]) / g2;
                        if (signal1 <= 0.)
                            Nph = 0.;
                        if (signal2 <= 0.) 
                            Ne = 0.;
                        if (yields.Lindhard > DBL_MIN && Nph > 0. && Ne > 0.) 
                        {
                            if(eeEnergy == true)
                                keV = (Nph + Ne) * Wq_eV * 1e-3; 
                            else
                                keV = (Nph + Ne) * Wq_eV * 1e-3 / yields.Lindhard;
                            
                        }
                        else
                            keV = 0.;
                    }
                    if(doBinning == 1)
                    {
                        indexS1 = (int)floor((signal1-minS1)/s1binWidth);
                        if(logS2 == 1)
                            indexS2 = (int)floor((log10(signal2)-log10(minS2))/s2binWidth);
                        else
                            indexS2 = (int)floor((signal2-minS2)/s2binWidth);
                        if(usePosition==1)
                        {
                            indexR = (int)floor( (sqrt(pow(smearPos[0],2)+pow(smearPos[1],2))-minR)/RbinWidth);
                            indexZ = (int)floor((smearPos[2]-minZ)/ZbinWidth);
                        }
                        s1s2RZbins[indexS1][indexS2][indexR][indexZ]+=1;

                    }
                    else            
                    {
                        stringstream tempString;
                        tempString << keV << "\t\t";
                        if(outputQuanta == true)
                            tempString << quanta.photons << "\t" << quanta.electrons << "\t" << (int)scint[0] << "\t" << (int)scint[1] << "\t" << (int)scint2[0] << "\t" << signal1 << "\t\t" << signal2 << "\t\t";
                        else
                            tempString << signal1 << "\t\t" << signal2 << "\t\t";
                        if(spec.doMigdal == 1 && verbosity == true)
                            tempString << migdalE[0] << "\t\t\t" << migdalE[1] << "\t\t\t" << int(migdalE[2]) << "\t\t";
                        if(useTiming==2)
                            tempString << scint2[9] << "\t\t";
                        if(usePosition==1)
                            tempString << fixed << setprecision(1) << smearPos[0] << "\t" << smearPos[1] << "\t" << smearPos[2] << "\t\t";//tempString << sqrt(pow(smearPos[0],2)+pow(smearPos[1],2)) << "\t" << smearPos[2] << "\t\t";
                        if(MCtruthE == false && verbosity == true)
                            tempString << keVtrue << "\t\t";
                        if(outputLindhard == 1)
                            tempString << yields.Lindhard << "\n";                
                        else
                            tempString << "\n";
                        outStream << tempString.str();
                    }
                    if (progress == 1 && fmod(((double) j/numEvents * 100),10) ==0)
                        cout << (double) j/numEvents * 100 << "%\n";
                }
                else
                {
                    j--;
                    numTrials--;
                }
            }
            else
                j--;

        }

    //for spectra from a file we have the option of outputting exposure or effective exposure for given even sample
    if( type == "file" )
    {
        if( exposure == -1 && ratioEvents == 0)
            outStream << "effective exposure: " << numTrials/spec.spline_spectrum_prep.totRate << endl;
        else if( exposure == -1 && type == "file" && ratioEvents == 1)
            outStream << "event efficiency ratio: " << (double)numEvents/numTrials << endl;
        else if( exposure == -1 && type == "RNC" )
            outStream << "event efficiency ratio: " << (double)numEvents/numTrials << endl;
        else
            outStream << "exposure: " << exposure << endl;
    }

    //output events binned in r and Z
    if(doBinning==1)
    {
        for(int l=0;l<numBinsZ;l++)
        {
            for(int k=0;k<numBinsR;k++)
            {
                if(usePosition==1)
                    outStream << "(rBin, ZBin) = (" << k << ", " << l << ")\n";
                for(int i=0;i<numBinsS1;i++)
                {
                    for(int j=0;j<numBinsS2;j++)
                        outStream << s1s2RZbins[i][j][k][l] << "\t";
                    outStream << "\n";     
                }
            }
        }    
    }

    outputFile.close();
    return 1;
}


//function for setting the detector parameters from a detector file
void setDetectorPars(string detectorName, Detector_def *detector)
{
    if(detectorName=="-1")
        return;
    std::ifstream RFF;

    RFF.open(detectorName,ifstream::in);
    if(!RFF)
    {
        RFF.open("detectors/"+detectorName,ifstream::in);
        if(!RFF)
        {
            cerr << "ERROR: file not found for detector " << detectorName << std::endl;
            assert(0);
        }
    }
  
    string line;
    vector<string> lineElements;
    double noise[4] = {0,0,0,0};
    while (!(RFF.eof()))
    {
        getline(RFF, line);
        istringstream linestream(line);
        for(string s; linestream >> s; )
            lineElements.push_back(s);
        if(lineElements.size()>3)
        {   
            if(lineElements[1] == "=")
            {
                lineElements[2].pop_back();
                if(lineElements[0] == "g1")
                    detector->set_g1(stof(lineElements[2]));
                if(lineElements[0] == "sPEres")
                    detector->set_sPEres(stof(lineElements[2]));
                if(lineElements[0] == "sPEthr")
                    detector->set_sPEthr(stof(lineElements[2]));
                if(lineElements[0] == "sPEeff")
                    detector->set_sPEeff(stof(lineElements[2]));
                if(lineElements[0] == "P_dphe")
                    detector->set_P_dphe(stof(lineElements[2]));
                if(lineElements[0] == "noise[0]")
                    noise[0]=stof(lineElements[2]);
                if(lineElements[0] == "noise[1]")
                    noise[1]=stof(lineElements[2]);
                if(lineElements[0] == "noise[2]")
                    noise[2]=stof(lineElements[2]);
                if(lineElements[0] == "noise[3]")
                    noise[3]=stof(lineElements[2]);
                if(lineElements[0] == "coinWind")
                    detector->set_coinWind(stof(lineElements[2]));
                if(lineElements[0] == "coinLevel")
                    detector->set_coinLevel(stoi(lineElements[2]));
                if(lineElements[0] == "numPMTs")
                    detector->set_numPMTs(stoi(lineElements[2]));
                if(lineElements[0] == "g1_gas")
                    detector->set_g1_gas(stof(lineElements[2]));
                if(lineElements[0] == "s2Fano")
                    detector->set_s2Fano(stof(lineElements[2]));
                if(lineElements[0] == "s2_thr")
                    detector->set_s2_thr(stof(lineElements[2]));
                if(lineElements[0] == "E_gas")
                    detector->set_E_gas(stof(lineElements[2]));
                if(lineElements[0] == "eLife_us")
                    detector->set_eLife_us(stof(lineElements[2]));
                if(lineElements[0] == "T_Kelvin")
                    detector->set_T_Kelvin(stof(lineElements[2]));
                if(lineElements[0] == "p_bar")
                    detector->set_p_bar(stof(lineElements[2]));
                if(lineElements[0] == "dtCntr")
                    detector->set_dtCntr(stof(lineElements[2]));
                if(lineElements[0] == "dt_min")
                    detector->set_dt_min(stof(lineElements[2]));
                if(lineElements[0] == "dt_max")
                    detector->set_dt_max(stof(lineElements[2]));
                if(lineElements[0] == "radius")
                    detector->set_radius(stof(lineElements[2]));
                if(lineElements[0] == "radmax")
                    detector->set_radmax(stof(lineElements[2]));
                if(lineElements[0] == "TopDrift")
                    detector->set_TopDrift(stof(lineElements[2]));
                if(lineElements[0] == "anode")
                    detector->set_anode(stof(lineElements[2]));
                if(lineElements[0] == "gate")
                    detector->set_gate(stof(lineElements[2]));
                if(lineElements[0] == "cathode")
                    detector->set_cathode(stof(lineElements[2]));
                if(lineElements[0] == "PosResExp")
                    detector->set_PosResExp(stof(lineElements[2]));
                if(lineElements[0] == "PosResBase")
                    detector->set_PosResBase(stof(lineElements[2]));
                if(lineElements[0] == "driftField")
                    detector->set_driftField(stof(lineElements[2]));
            }        
        }
        detector->set_noiseBaseline(noise[0],noise[1],noise[2],noise[3]);
        lineElements.clear();
    }
    RFF.close();

}

//function for setting analysis parameters from an analysis file
void setAnalysisPars(string analysisFilename)
{
    if(analysisFilename=="-1")
        return;

    std::ifstream RFF;

    RFF.open(analysisFilename,ifstream::in);
 
    if(!RFF)
    {
        RFF.open("inc/"+analysisFilename,ifstream::in);
        if(!RFF)
        {
            cerr << "file not found for analysis " << analysisFilename << std::endl;
            assert(0);
        }
    }
    
    string line;
    vector<string> lineElements;
    while (!(RFF.eof()))
    {
        getline(RFF, line);
        istringstream linestream(line);
        for(string s; linestream >> s; )
            lineElements.push_back(s);
        if(lineElements.size()>3)
        {   
            
            if(lineElements[2] == "=")
            {
                std::istringstream is(lineElements[3]);
                bool b;
                is >> std::boolalpha >> b;
                if(lineElements[1] == "verbosity")
                    verbosity = b;
                if(lineElements[1] == "MCtruthE")
                    MCtruthE = b;
                if(lineElements[1] == "MCtruthPos")
                    MCtruthPos = b;
                if(lineElements[1] == "doFiducialCut")
                    doFiducialCut = b;
                if(lineElements[1] == "useTiming")
                    useTiming = stoi(lineElements[3]);
                if(lineElements[1] == "usePD")
                    usePD = stoi(lineElements[3]);
                if(lineElements[1] == "useS2")
                    useS2 = stoi(lineElements[3]);
                if(lineElements[1] == "minS1")
                    minS1 = stof(lineElements[3]);
                if(lineElements[1] == "maxS1")
                    maxS1 = stof(lineElements[3]);
                if(lineElements[1] == "numBinsS1")
                    numBinsS1 = stoi(lineElements[3]);
                if(lineElements[1] == "numBinsS2")
                    numBinsS2 = stoi(lineElements[3]);
                if(lineElements[1] == "numBinsZ")
                    numBinsZ = stoi(lineElements[3]);
                if(lineElements[1] == "numBinsR")
                    numBinsR = stoi(lineElements[3]);
                if(lineElements[1] == "minZ")
                    minZ = stoi(lineElements[3]);
                if(lineElements[1] == "maxZ")
                    maxZ = stoi(lineElements[3]);
                if(lineElements[1] == "minR")
                    minR = stoi(lineElements[3]);
                if(lineElements[1] == "maxR")
                    maxR = stoi(lineElements[3]);
                if(lineElements[1] == "minS2")
                    minS2 = stof(lineElements[3]);
                if(lineElements[1] == "maxS2")
                    maxS2 = stof(lineElements[3]);
                if(lineElements[1] == "logS2")
                    logS2 = stof(lineElements[3]);
                if(lineElements[1] == "z_step")
                    z_step = stof(lineElements[3]);
                if(lineElements[1] == "E_step")
                    E_step = stof(lineElements[3]);
                if(lineElements[1] == "useCorrected")
                    useCorrected = stoi(lineElements[3]);
            }
        }
        lineElements.clear();
    }
    RFF.close();
}
