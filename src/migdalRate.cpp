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
double Mn = 931493;
double MtXe = 131.29*931493;
double Me = 510.999;
double isoFracTable[7] = {0.0191421,0.283724,0.324514,0.536981,0.806574,0.911205,1};
double isoTable[7] = {128,129,130,131,132,134,136};
double nl_energy[6][5] = {{    -1,    -1,    -1,-1},
                          {    35,    -1,    -1,-1},    
                          {   5.4,   4.9,    -1,-1},
                          {   1.1,  0.93,  0.66,-1},
                          {   0.2,  0.14,6.1e-2,-1},
                          {2.1e-2,9.8e-3,    -1,-1},
}; //energy levels for xenon

int Nlmax[6] = {-1,1,2,3,3,2};

gsl_spline *Znl_spline[6][3];
gsl_interp_accel *Znl_accel[6][3];
gsl_spline *Znl_int_spline[6][3];
gsl_interp_accel *Znl_int_accel[6][3];
double Znl_y_max[6][3];
double Znl_x_max[6][3];
double Znl_x_min[6][3];

double mu(double x, double y)
{
    return x*y/(x+y);
}
   
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
    
    //create integrated probabilities to a maximum Ee
    double Znl_int[251];
    for(n=1; n<6; n++)
    {
        for(l=0; l<Nlmax[n]; l++)
        {
            for(int j=0; j<250; j++)
            {
                Znl_int[j] = gsl_spline_eval_integ(Znl_spline[n][l], Znl_x_min[n][l], EMnl[n][l][j], Znl_accel[n][l]);
            }
            Znl_int_spline[n][l] = gsl_spline_alloc(gsl_interp_linear, 251);
            Znl_int_accel[n][l] = gsl_interp_accel_alloc();
            gsl_spline_init(Znl_int_spline[n][l], EMnl[n][l], Znl_int, 251);
        }
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

//returns total prob below max_Ee
double Z_nl_integrated(int N, int l, double Qe, double max_Ee)
{ 
    if ( max_Ee > 69.9 )
    {
        return pow((Qe / 0.001),2) * gsl_spline_eval(Znl_int_spline[N][l], 69.9, Znl_int_accel[N][l]);
    }
    else if ( max_Ee < .001 )
        return 0;
    else
        return pow((Qe / 0.001),2) * gsl_spline_eval(Znl_int_spline[N][l], max_Ee, Znl_int_accel[N][l]);
}

//returns maximum diff probability of ionization for a given level
double Znl_max(int N, int l, double Qe)
{
    return pow((Qe / 0.001),2) * Znl_y_max[N][l];
}

//returns a randomly selected xenon isotope
double rand_xe_isotope() 
{
    double r = RandomGen::rndm()->rand_uniform();
    for (int i=0;i<8;i++)
    {   
        if (r < isoFracTable[i])
            return isoTable[i];
    }
    return 131.3; //in case something goes wrong
}

//eject a random electron? (proportional to prob)
vector<double> rand_migdalE(double ERnr, int mig_type, double monoE)
{
    
    MtXe = Mn*rand_xe_isotope();
    double maxEM=0;
    if (mig_type == 1)
        maxEM = 2*sqrt(monoE*ERnr*MtXe/Mn)-ERnr*(MtXe+Mn)/Mn;
    if (mig_type == 2)
        maxEM = 100;

    double Ee,EeMax;

    double t,prob=0;
    t = RandomGen::rndm()->rand_uniform();
    for(int N=1; N<6; N++)
    {
        for (int L=0; L<Nlmax[N]; L++)
        {
            EeMax = maxEM - nl_energy[N][L];
            if (EeMax < 0)
                continue;
            prob+=Z_nl_integrated(N,L,qe(ERnr),EeMax);
            if(t<prob)
            {
                Ee = EeMax*RandomGen::rndm()->rand_uniform();
                while(RandomGen::rndm()->rand_uniform()*Znl_max(N, L, qe(ERnr)) > Z_nl(N, L, qe(ERnr), Ee))
                { 
                    Ee = EeMax*RandomGen::rndm()->rand_uniform();
                }
                return {Ee,nl_energy[N][L],(double)N};
            }
        }
    }
    return {0,0,0};

}

double ERmaxNeutronEM(double En, double EM)
{
    return pow(mu(MtXe,Mn),2)/(MtXe*Mn)*En*(pow(1-sqrt(1-EM/(mu(MtXe,Mn)*En/Mn)),2)+4*sqrt(1-EM/(mu(MtXe,Mn)*En/Mn)));
}

double migdalIntegrand(double ErNR, void *pars)
{
    TestSpectra::SPLINE_migdal_prep *migSpec = (TestSpectra::SPLINE_migdal_prep *) pars;
    double Edet = migSpec->ER;
    double integrand=0;
    double dRnr;
    int N = migSpec->N;
    
    if( ErNR > ERmaxNeutronEM(migSpec->monoE, Edet - ErNR*Leff))
        return 0;
    
    dRnr = gsl_spline_eval(migSpec->NR_spec->spectrumSpline, ErNR, migSpec->NR_spec->accelSS);
    
    for( int l=0; l<Nlmax[N]; l++ )
    {
        if ( Edet - ErNR*Leff - nl_energy[N][l] > 0 )
            integrand += Z_nl( N, l, qe(ErNR), Edet - ErNR*Leff - nl_energy[N][l]);
    }
    return integrand * dRnr;
}

double ERmaxNeutron(double Edet, double En)
{
    return 4*En*Mn*MtXe/pow(Mn+MtXe,2);
}

double dRdEmigdalNeutron(double Edet, TestSpectra::SPLINE_migdal_prep *migSpec)
{
    
    gsl_function F;
    F.function = &migdalIntegrand;
    F.params = migSpec;
    migSpec->ER = Edet;
    int limit = 1000;
	double integral,absErr,tol;
    tol=1e-6;

    W = gsl_integration_workspace_alloc (1000);

    double ErMin = migSpec->NR_spec->xMin; //This is effectively zero for all regions we care about
    double ErMax = ERmaxNeutron( Edet, migSpec->monoE);//migSpec->NR_spec->xMax;
    if(ErMax > Edet/Leff)
        ErMax = Edet/Leff;
    if(ErMax > migSpec->NR_spec->xMax)
        ErMax = migSpec->NR_spec->xMax;
    
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
    migdalSpec.monoE = NRspec->monoE;
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
            Rate[i] = dRdEmigdalNeutron( Er[i], &migdalSpec);
            RFF << Er[i] << " " << Rate[i] << endl;
        }
        RFF.close();

        migdalSpec.xMin = Er[0];
        migdalSpec.xMax = Er[1999];

        //gsl_spline_init( migdalSpec.migdalSpline, Er, Rate, 2000);
        //migdalSpec.totRate=gsl_spline_eval_integ(migdalSpec.migdalSpline, migdalSpec.xMin, migdalSpec.xMax, migdalSpec.accelMS);
    }

}
