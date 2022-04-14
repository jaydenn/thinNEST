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

double Z_nl_integrated(int N, int l, double Qe, double max_Ee);

gsl_integration_workspace * W;

//set some physical constants
double Leff = 0.15;          //approximate quenching factor
double Mn = 939565.42;       //neutron mass in keV
double AMU = Mn/1.00727;     //1 atomic mass unit in keV
double MtXe = 131.29*931493; //RAM of xenon in keV
double Me = 510.999;         //electron mass
double isoFracTable[7] = {0.0191421,0.283724,0.324514,0.536981,0.806574,0.911205,1}; //naturally occuring isotopic abundances of xenon
double isoTable[7] = {128,129,130,131,132,134,136}; //atomic numbers of isotopes
double nl_energy[6][5] = {{    -1,    -1,    -1,-1}, //the atomic binding energies of shells indexed by [n][l]
                          {    35,    -1,    -1,-1},    
                          {   5.4,   4.9,    -1,-1},
                          {   1.1,  0.93,  0.66,-1},
                          {   0.2,  0.14,6.1e-2,-1},
                          {2.1e-2,9.8e-3,    -1,-1},
}; 
int Nlmax[6] = {-1,1,2,3,3,2}; 

//setup interpolation objects
gsl_spline *Znl_spline[6][3];
gsl_interp_accel *Znl_accel[6][3];
gsl_spline *Znl_int_spline[6][3];
gsl_interp_accel *Znl_int_accel[6][3];
gsl_spline *invZnl_int_spline[6][3];
gsl_interp_accel *invZnl_int_accel[6][3];

double Znl_y_max[6][3];
double Znl_x_max[6][3];
double Znl_x_min[6][3];
double maxCumulativeProb=1;

//reduced mass definition
double mu(double x, double y)
{
    return x*y/(x+y);
}

//electron recoil momentum
double qe(double ERnr, double Mt)
{
    return Me*sqrt(2 * ERnr / Mt);
}

//max recoil energy for a neutron of energy En
double ERmaxNeutron(double En, double Mt)
{
    return 4*En*Mn*Mt/pow(Mn+Mt,2);
}

//max electromagnetic energy for a neutron of energy En
double EMmaxNeutron(double En, double Mt)
{
    return mu(Mn,Mt)*En/Mn;
}

//max electromagnetic energy for a given nuclear recoil energy ENR and neutron energy En
double EMmaxNeutronER(double En, double ENR, double Mt)
{
    return 2.*sqrt(En*ENR*Mt/Mn) - ENR*(Mt+Mn)/Mn;
}

double EMmaxNeutrinoER(double Enu, double ENR, double Mt)
{
    return Enu-ENR;
}

double ERmaxNeutronEM(double En, double EM, double Mt)
{
    return pow(mu(Mt,Mn),2)/(Mt*Mn)*En*(pow(1-sqrt(1-EM/(mu(Mt,Mn)*En/Mn)),2)+4*sqrt(1-EM/(mu(Mt,Mn)*En/Mn)));
}

double ERmaxNeutrinoEM(double Enu, double EM, double Mt)
{
    return pow(2.*Enu - EM,2)/(2.*(2.*Enu+Mt));
}

//calculate the total migdal probability for a given recoil
double totalMigProb(double maxE, double ENR, double Mt, int mig_type)
{
    double maxEM=0, maxEe=0;
    if (mig_type == 1)
        maxEM  = EMmaxNeutronER(maxE, ENR, Mt);
    if (mig_type == 2)
        maxEM = EMmaxNeutrinoER(maxE, ENR, Mt);
    if (mig_type == 3)
    {
        maxEM = 0; //need to define
    }

    double prob=0;
    for(int N=1; N<6; N++)
    {
        for (int L=0; L<Nlmax[N]; L++)
        {
            maxEe = maxEM - nl_energy[N][L];
            if (maxEe < 0)
                continue;
            prob += Z_nl_integrated(N,L,qe(ENR,Mt),maxEe);
        }
    }
    return prob;
}

//initialize the Migdal ionization probabilities (read in from file obtained from Ibe et al.)
void init_Znl(double maxE, int mig_type, int migdalOptimize)
{
    
    ifstream RFF("src/Xe_new.dat");
    if(!RFF.is_open())
    {
        cout << "error opening migdal data file" << endl;
        return;
    }
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
            for(int j=0; j<251; j++)
            {
                Znl_int[j] = gsl_spline_eval_integ(Znl_spline[n][l], Znl_x_min[n][l], EMnl[n][l][j], Znl_accel[n][l]);
            }
            Znl_int_spline[n][l] = gsl_spline_alloc(gsl_interp_linear, 251);
            Znl_int_accel[n][l] = gsl_interp_accel_alloc();
            gsl_spline_init(Znl_int_spline[n][l], EMnl[n][l], Znl_int, 251);
            
            invZnl_int_spline[n][l] = gsl_spline_alloc(gsl_interp_linear, 251);
            invZnl_int_accel[n][l] = gsl_interp_accel_alloc();
            gsl_spline_init(invZnl_int_spline[n][l], Znl_int, EMnl[n][l], 251);
        }
    }
    
    if(RFF.bad())
        perror("error while reading file ");
    RFF.close();
    
    double maxENR = 0;
    if (mig_type == 1)
    {
        maxENR = ERmaxNeutron(maxE,128*AMU);
    }
    if (mig_type == 2)
    {
        maxENR =  ERmaxNeutrinoEM(maxE,0,128*AMU);
    }    
    if (mig_type == 3)
    {
        maxENR = 0; //need to define
    }

    //find maximum probability and normalize to it
    if(migdalOptimize==1)
    {
        cout << "optimizing migdal probabilities..";
        double maxProb = 0;
        double nr = .001;
        while (totalMigProb(maxE, nr, 128*AMU, mig_type) > maxProb && nr < maxENR)
        {
            maxProb = totalMigProb(maxE, nr, 128*AMU, mig_type);
            nr*=1.001;
        }
        maxCumulativeProb = maxProb;
        cout << " Migdal probabilities *~= " << 1/maxCumulativeProb << endl;
    }

}

//function for the ionization probabilities
double Z_nl(int N, int l, double Qe, double Ee)
{
    if ( Ee > 70 )
    {
        return 0;
    }
    else if ( Ee < .001 )
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_spline[N][l], 0.001, Znl_accel[N][l]);
    else
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_spline[N][l], Ee, Znl_accel[N][l]);
}

//returns total prob below max_Ee
double Z_nl_integrated(int N, int l, double Qe, double max_Ee)
{ 
    if ( max_Ee > 69.9 )
    {
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_int_spline[N][l], 69.9, Znl_int_accel[N][l]);
    }
    else if ( max_Ee < .001 )
        return 0;
    else
        return (pow((Qe / 0.001),2) ) * gsl_spline_eval(Znl_int_spline[N][l], max_Ee, Znl_int_accel[N][l]);
}

//returns total prob below max_Ee
double invZ_nl_integrated(int N, int l, double Qe, double Znl)
{ 
    if(Znl>Z_nl_integrated( N, l, Qe, 70))
    {    cout << "Prob. too high result may be skewed" << endl; return 70; }
    return gsl_spline_eval(invZnl_int_spline[N][l], Znl/pow((Qe / 0.001),2), invZnl_int_accel[N][l]);
}

//returns maximum diff probability of ionization for a given level
double Znl_max(int N, int l, double Qe)
{
    return pow((Qe / 0.001),2) * Znl_y_max[N][l];
}

//returns a randomly selected xenon isotope
double rand_xe_isotope(double ENR, double maxE, int mig_type) 
{
    isoLoop:
        double r = RandomGen::rndm()->rand_uniform();
        double ERmax=0;
        for (int i=0;i<8;i++)
        {   
            if (r < isoFracTable[i])
            {
                if(mig_type==1)
                    ERmax = ERmaxNeutron(maxE, Mn*isoTable[i]);
                if(mig_type==2)
                    ERmax = ERmaxNeutrinoEM(maxE, 0, Mn*isoTable[i]);
                if( ENR < ERmax ) //ensures kinematic consistency
                    return isoTable[i];
                if( ENR > ERmaxNeutron(maxE, Mn*129) && mig_type == 1) //speeds up case where only 128 is consistent
                    return 128;
                else
                    goto isoLoop;
            }
        }
        return 131.3; //in case something goes wrong
}

//eject a random electron? (proportional to prob) - main Migdal MC function
vector<double> rand_migdalE(double ERnr, int mig_type, double maxE)
{
    MtXe = AMU*131.3;
    double maxEM=0;
    if (mig_type == 1)
        maxEM = EMmaxNeutronER(maxE,ERnr,MtXe);
    if (mig_type == 2)
        maxEM = EMmaxNeutrinoER(maxE,ERnr,MtXe);
    if (mig_type == 3)
        maxEM = 100;

    double Ee,EeMax;

    double t,prob=0;
    double maxProb=totalMigProb(maxE, ERnr, MtXe, mig_type)/maxCumulativeProb;

    t = RandomGen::rndm()->rand_uniform();
    if ( t < maxProb )
    {
        for(int N=1; N<6; N++)
        {
            for (int L=0; L<Nlmax[N]; L++)
            {
                EeMax = maxEM - nl_energy[N][L];
                if (EeMax < 0)
                    continue;
                prob+=Z_nl_integrated(N,L,qe(ERnr,MtXe),EeMax)/maxCumulativeProb;
                if( t < prob )
                {
                    t = RandomGen::rndm()->rand_uniform()*Z_nl_integrated(N,L,qe(ERnr,MtXe),EeMax);
                    Ee = invZ_nl_integrated(N, L, qe(ERnr,MtXe), t);
                    return {Ee,nl_energy[N][L],(double)N,MtXe/AMU};
                }
            }
        }
    }
    return {0,0,0,0};

}

//the following 3 functions are for numerical calculation of the Migdal spectrum
double migdalIntegrand(double ErNR, void *pars)
{
    TestSpectra::SPLINE_migdal_prep *migSpec = (TestSpectra::SPLINE_migdal_prep *) pars;
    double Edet = migSpec->ER;
    double integrand=0;
    double dRnr;
    int N = migSpec->N;
    
    if( ErNR > ERmaxNeutronEM(migSpec->monoE, Edet - ErNR*Leff,MtXe))
        return 0;
    
    dRnr = gsl_spline_eval(migSpec->NR_spec->spectrumSpline, ErNR, migSpec->NR_spec->accelSS);
    
    for( int l=0; l<Nlmax[N]; l++ )
    {
        if ( Edet - ErNR*Leff - nl_energy[N][l] > 0 )
            integrand += Z_nl( N, l, qe(ErNR,MtXe), Edet - ErNR*Leff - nl_energy[N][l]);
    }
    return integrand * dRnr;
}

//differential Migdal spectrum
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
    double ErMax = ERmaxNeutron( migSpec->monoE, MtXe);//migSpec->NR_spec->xMax;
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

//calculate spectrum and write to file
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
