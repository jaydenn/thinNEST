using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "TestSpectra.hh"

gsl_integration_workspace * W;
double Leff = 0.15;
double MtXe = 131.29*931493;
double Me = 510.999;
double nl_energy[6][5] = {{    -1,    -1,    -1,-1},
                          {    35,    -1,    -1,-1},
                          {   5.4,   4.9,    -1,-1},
                          {   1.1,  0.93,  0.66,-1},
                          {   0.2,  0.14,6.1e-2,-1},
                          {2.1e-2,9.8e-3,    -1,-1},
}; //energy levels for xenon
double Pnl[6][3] = {  {     0,     0,     0},
                      {4.6e-6,     0,     0},
                      {2.9e-5,1.3e-4,     0},
                      {8.7e-5,5.2e-4,3.5e-3},
                      {3.4e-4,1.4e-3,3.4e-2},
                      {4.1e-4,   0.1,     0},
}; //total ionization prob at qe=0.001*Me

int Nlmax[6] = {-1,1,2,3,3,2};

gsl_spline *Znl_spline[6][3];
gsl_interp_accel *Znl_accel[6][3];
double Znl_y_max[6][3];
double Znl_x_max[6][3];
double Znl_x_min[6][3];

void init_Znl()
{
    
    ifstream RFF("src/Xe.dat");
    
    double*** EMnl = new double**[6];
    double*** Znl = new double**[6];
    for (int ni=1; ni < 6; ni++)
    {
        EMnl[ni] = new double*[3]; 
        Znl[ni] = new double*[3];
    }

    int n,l;
    string line;
    while (getline(RFF,line))
    {
        stringstream ss(line);
        ss >> n >> l;
        EMnl[n][l] = new double[251];
        Znl[n][l] = new double[251];
        Znl_y_max[n][l] = 0;

        int i = 0;
        while (i<251)
        {
            getline(RFF,line);
            stringstream ss(line);
            ss >> EMnl[n][l][i] >> Znl[n][l][i];
            EMnl[n][l][i]/=1000;   //convert units
            Znl[n][l][i]*=1000/2/M_PI;
            if ( Znl[n][l][i] > Znl_y_max[n][l] )
                Znl_y_max[n][l] = Znl[n][l][i];
            i++;         
        }
        Znl_spline[n][l] = gsl_spline_alloc(gsl_interp_linear, 251);
        Znl_accel[n][l] = gsl_interp_accel_alloc();
        gsl_spline_init(Znl_spline[n][l], EMnl[n][l], Znl[n][l], 251);
        Znl_x_min[n][l] = EMnl[n][l][0];
        Znl_x_max[n][l] = EMnl[n][l][250];

    }
    if(RFF.bad())
        perror("error while reading file ");
    RFF.close();
}

double qe(double ERnr)
{
    return Me*sqrt(2 * ERnr / MtXe);
}

double Z_nl(int N, int l, double Qe, double Ee)
{
    if ( Ee > 70 )
    {
        return 0;
    }
    else if ( Ee < .001 )
        return pow((Qe / 0.001),2) * gsl_spline_eval(Znl_spline[N][l], 0.001, Znl_accel[N][l]);
    else
        return pow((Qe / 0.001),2) * gsl_spline_eval(Znl_spline[N][l], Ee, Znl_accel[N][l]);
}

vector<int> rand_nl(double ERnr)
{
    double r = RandomGen::rndm()->rand_uniform();
    double p = 0;
    vector <int> NL = {-1,-1};

    for(int N=1; N<6; N++)
    {
        for (int L=0; L<Nlmax[N]; L++)
        {
            p += pow(qe(ERnr)/(Me*0.001),2) * Pnl[N][L];
            if(r<p)
            {
                NL[0]=N;
                NL[1]=L;
                return NL;
            }
        }
    }
    return NL;
}

vector<double> rand_migdalE(double ERnr)
{
    
    vector<double> energies = {0,0};
    //get a random electron (proportional to prob)
    vector<int> NL = rand_nl(ERnr);

    if (NL[0]>0)
    {
        energies[1] = nl_energy[NL[0]][NL[1]];
        
        double yMax = pow((qe(ERnr) / 0.001),2)*Znl_y_max[NL[0]][NL[1]];
        double xMin = Znl_x_min[NL[0]][NL[1]];
        double xMax = Znl_x_max[NL[0]][NL[1]];
        double FuncValue;

        vector<double> xyTry = {
          xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
          yMax * RandomGen::rndm()->rand_uniform(), 1.};
        while (xyTry[2] > 0.)
        {
            FuncValue = Z_nl( NL[0], NL[1], qe(ERnr),  xyTry[0]);
            xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, 0., yMax, xyTry[0],
                                              xyTry[1], FuncValue);
        }
        energies[0] = xyTry[0];
    }

    return energies;

}

double migdalIntegrand(double ErNR, void *pars)
{
    TestSpectra::SPLINE_migdal_prep *migSpec = (TestSpectra::SPLINE_migdal_prep *) pars;
    double Edet = migSpec->ER;
    double integrand=0;
    double dRnr;
    int N = migSpec->N;
    
    if (ErNR > migSpec->NR_spec->xMax )
        return 0;
    else if(ErNR < migSpec->NR_spec->xMin)
        dRnr = gsl_spline_eval(migSpec->NR_spec->spectrumSpline, migSpec->NR_spec->xMin, migSpec->NR_spec->accelSS);
    else
    {
        dRnr = gsl_spline_eval(migSpec->NR_spec->spectrumSpline, ErNR, migSpec->NR_spec->accelSS);
    }
    
    for( int l=0; l<Nlmax[N]; l++ )
    {
        if ( Edet - ErNR*Leff - nl_energy[N][l] > 0 )
            integrand += Z_nl( N, l, qe(ErNR), Edet - ErNR*Leff - nl_energy[N][l]);
    }
    return integrand * dRnr;
}


double dRdEmigdal(double Edet, TestSpectra::SPLINE_migdal_prep *migSpec)
{
    
    gsl_function F;
    F.function = &migdalIntegrand;
    F.params = migSpec;
    migSpec->ER = Edet;
    int limit = 1000;
	double integral,absErr,tol;
    tol=1e-6;

    W = gsl_integration_workspace_alloc (1000);

    double ErMin = migSpec->NR_spec->xMin;
    double ErMax = migSpec->NR_spec->xMax;
    if(ErMax > Edet/Leff)
        ErMax = Edet/Leff;

    gsl_integration_qag(&F, ErMin, ErMax, tol, 1e-3, limit, 6, W, &integral, &absErr);
    
    if(absErr/integral > .01)
    {
        gsl_integration_qag(&F, ErMin, ErMax, integral/100, 1e-3, limit, 6, W, &integral, &absErr);
        if(absErr/integral > .01)
            cerr << "warning, migdal integral may not be accurate ~" << absErr/integral*100 << "%" << endl;
    }    
    return integral; 

}

void calcMigdalSpectrum(TestSpectra::SPLINE_spectrum_prep *NRspec)
{
    TestSpectra::SPLINE_migdal_prep migdalSpec;
    migdalSpec.NR_spec = NRspec;
    migdalSpec.migdalSpline = gsl_spline_alloc(gsl_interp_linear,2000);
    migdalSpec.accelMS = gsl_interp_accel_alloc();

    double Er[2000];
    double Rate[2000];
    
    std::fstream RFF;
    
    for (int N=2; N<6; N++)
    {
        migdalSpec.N = N;

        string filename = NRspec->filename;
        string ins("mig"+to_string(N));
        //filename.erase( filename.find("s"),8);
        filename.replace( filename.find("NR"),2,ins);
        RFF.open(filename,fstream::out);

        double e_min[6] = {-1,35,4.9,0.66,.061,.0098};
        double interval = pow(50/e_min[N],0.0005);
        Er[0]=e_min[N];
        for (int i=0; i<2000; i++)
        {
            Er[i] = interval*Er[(i-1>=0 ? i-1 : 0)];
            Rate[i] = dRdEmigdal( Er[i], &migdalSpec);
            RFF << Er[i] << " " << Rate[i] << endl;
        }
        RFF.close();

        migdalSpec.xMin = Er[0];
        migdalSpec.xMax = Er[1999];

        //gsl_spline_init( migdalSpec.migdalSpline, Er, Rate, 2000);
        //migdalSpec.totRate=gsl_spline_eval_integ(migdalSpec.migdalSpline, migdalSpec.xMin, migdalSpec.xMax, migdalSpec.accelMS);
    }

}
