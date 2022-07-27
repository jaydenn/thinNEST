#include <string>

#include "Target.hh"
//#include "Spectra.hh"
#include "PhysicalConstants.hh"

using namespace std;
using std::vector;

//reduced mass definition
double mu(double x, double y)
{
    return x*y/(x+y);
}

//electron recoil momentum
double Ve(double ERnr, double MT)
{
    return sqrt(2 * ERnr / MT);
}

double Target::ERmax(double Ein, double EM, double MT)
{
    switch(migType)
    {
        case NEUTRON:
            return pow(mu(MT,Mn),2)/(MT*Mn)*Ein*(pow(1-sqrt(1-EM/(mu(MT,Mn)*Ein/Mn)),2)+4*sqrt(1-EM/(mu(MT,Mn)*Ein/Mn)));
        case NEUTRINO:
            return pow(2.*Ein - EM,2)/(2.*(2.*Ein+MT));
        case WIMPmig:
            return pow(mu(wimpMass,MT),2)/MT * pow(wimpVMAX,2) *(1 - EM/(mu(wimpMass,MT) * pow(wimpVMAX,2)) + sqrt(1 - 2*EM/(mu(wimpMass,MT) * pow(wimpVMAX,2))));
        default:
            return -1;
    }
    
}
double Target::ERmax(double Ein, double EM)
{
    return ERmax( Ein, EM, Mt);
}

double Target::EMmax(double Ein, double MT)
{
    switch(migType)
    {
        case NEUTRON:
            return mu(Mn,Mt)*Ein/Mn;
        case NEUTRINO:
            return Ein;
        case WIMPmig:
            return mu(wimpMass,MT)*pow(wimpVMAX,2)/2.0;
        default:
            return -1;
    }
    
}
double Target::EMmax(double Ein)
{
    return EMmax( Ein, Mt);
}
    
double Target::EMmaxER(double Ein, double ENR, double MT)
{
    switch(migType)
    {
        case NEUTRON:
            return 2.*sqrt(Ein*ENR*MT/Mn) - ENR*(MT+Mn)/Mn;
        case NEUTRINO:
            return Ein-ENR;
        case WIMPmig:
            return sqrt(2*MT*pow(wimpVMAX,2))-ENR*(MT+wimpMass)/wimpMass;
        default:
            return -1;
    }
}

double Target::EMmaxER(double Ein, double ENR)
{
    return EMmaxER( Ein, ENR, Mt);
}
    
//returns total prob below max_Ee
double Target::Z_nl_integrated(int orbit_i, double v, double max_Ee)
{ 
    if ( max_Ee > 19.99 )
    {
        return (pow((v / 0.0001),2) ) * gsl_spline_eval(Znl_int_spline[orbit_i], 19.99, Znl_int_accel[orbit_i]);
    }
    else if ( max_Ee < .0001 )
        return 0;
    else
        return (pow((v / 0.0001),2) ) * gsl_spline_eval(Znl_int_spline[orbit_i], max_Ee, Znl_int_accel[orbit_i]);
}

//returns total prob below max_Ee
double Target::invZ_nl_integrated(int orbit_i, double v, double Znl)
{ 
    if(Znl>Z_nl_integrated( orbit_i, v, 20))
    {    cout << "Prob. too high result may be skewed" << endl; return 20; }
    return gsl_spline_eval(invZnl_int_spline[orbit_i], Znl/pow((v / 0.0001),2), invZnl_int_accel[orbit_i]);
}

//returns maximum diff probability of ionization for a given level
double Target::Z_nl_max(int orbit_i, double v)
{
    return pow((v / 0.0001),2) * Znl_y_max[orbit_i];
}


int Target::init()
{

    //Set target element
    if (elementZ == XENON)
    {
        target = "Xe";
        Mt = MtXe;
        isoN = isoNXe;
        isoTable = isoTableXe;
        isoFracTable = isoFracTableXe;
        orbits = orbitsXe;
        orbits_energy = orbits_energyXe;
        Norbits = NorbitsXe;
        Leff = LeffXe;
    }
    else if (elementZ == ARGON)
    {
        target = "Ar";
        Mt = MtAr;
        isoN = isoNAr;
        isoTable = isoTableAr;
        isoFracTable = isoFracTableAr;
        orbits = orbitsAr;
        orbits_energy = orbits_energyAr;
        Norbits = NorbitsAr;
        Leff = LeffAr;
    }
    else
    {
        cout << "target not found\n";
        return -1;
    }
    
    switch(migType)
    {
        case WIMPmig:
        {
            wimpMass = maxE;
            maxE = 0.5*wimpMass*pow(wimpVMAX,2); //kinetic energy in keV
        }
    }
    
    double** EMnl = new double*[Norbits];
    double** Znl = new double*[Norbits];
    
    int orbit_i=0;
    string line, orbitFile;
    while (orbit_i < Norbits)
    {    
        orbitFile = "./dat/"+orbits[orbit_i]+"_"+target+".dat";
        //cout << "reading: " << orbitFile << endl;
        std::fstream RFF(orbitFile);

        if(!RFF.is_open())
        {
            cout << "error opening file " << orbitFile << endl;
            return -1;
        }
        
        EMnl[orbit_i] = new double[100];
        Znl[orbit_i] = new double[100];
        Znl_y_max[orbit_i] = 0;
        int i=0;
        while (i<100)
        {
            getline(RFF,line);
            stringstream ss(line);
            ss >> EMnl[orbit_i][i] >> Znl[orbit_i][i];
            //EMnl[orbit_i][i]/=1000;   //convert units
            //Znl[orbit_i][i]*=1000/2/M_PI;
            if ( Znl[orbit_i][i] > Znl_y_max[orbit_i] )
                Znl_y_max[orbit_i] = Znl[orbit_i][i];
           i++;        
        }
        Znl_spline[orbit_i] = gsl_spline_alloc(gsl_interp_linear, 100);
        Znl_accel[orbit_i] = gsl_interp_accel_alloc();
        gsl_spline_init(Znl_spline[orbit_i], EMnl[orbit_i], Znl[orbit_i], 100);
        Znl_x_min[orbit_i] = EMnl[orbit_i][0];
        Znl_x_max[orbit_i] = EMnl[orbit_i][99];
        orbit_i++;  
        
        if(RFF.bad())
            perror("error while reading file ");
        RFF.close();
    }
    
    //create integrated probabilities to a maximum Ee
    double Znl_int[100];
   
    for(int orbit_i=0; orbit_i<Norbits; orbit_i++)
    {
            for(int j=0; j<100; j++)
            {
                Znl_int[j] = gsl_spline_eval_integ(Znl_spline[orbit_i], Znl_x_min[orbit_i], EMnl[orbit_i][j], Znl_accel[orbit_i]);
            }
            Znl_int_spline[orbit_i] = gsl_spline_alloc(gsl_interp_linear, 100);
            Znl_int_accel[orbit_i] = gsl_interp_accel_alloc();
            gsl_spline_init(Znl_int_spline[orbit_i], EMnl[orbit_i], Znl_int, 100);

            invZnl_int_spline[orbit_i] = gsl_spline_alloc(gsl_interp_linear, 100);
            invZnl_int_accel[orbit_i] = gsl_interp_accel_alloc();
            gsl_spline_init(invZnl_int_spline[orbit_i], Znl_int, EMnl[orbit_i], 100);
    }
    
    double maxENR = 0;
    maxENR = ERmax(maxE,0,isoTable[0]*AMU);
   
    //find maximum probability and normalize to it
    if(migOptimize==1)
    {
        cout << "optimizing migdal probabilities..";
        double nr = .001;
        double maxProb = totalMigProb(maxE, nr, AMU*isoTable[0]);        
        
        while (totalMigProb(maxE, nr, AMU*isoTable[0]) > maxProb && nr < maxENR)
        {
            maxProb = totalMigProb(maxE, nr, AMU*isoTable[0]);
            nr*=1.001;
        }
        maxCumulativeProb = maxProb;
        cout << " Migdal probabilities *~= " << 1/maxCumulativeProb << endl;
    }
    
    return 0;
}


//calculate the total migdal probability for a given recoil
double Target::totalMigProb(double maxE, double ENR, double MT)
{
    double maxEM=0, maxEe=0;
    maxEM  = EMmaxER(maxE, ENR);

    double prob=0;
    for(int orbit_i=0; orbit_i<Norbits; orbit_i++)
    {
        maxEe = maxEM - orbits_energy[orbit_i];
        
        if (maxEe < 0)
            continue;
        prob += Z_nl_integrated(orbit_i,Ve(ENR,MT),maxEe);
    }

    return prob;
}

//function for the ionization probabilities
double Target::Z_nl(int orbit_i, double v, double Ee)
{
    if ( Ee > 19.99 || Ee < .0001001)
        return 0;
    else
        return (pow((v / 0.0001),2) ) * gsl_spline_eval(Znl_spline[orbit_i], Ee, Znl_accel[orbit_i]);
}


//eject a random electron? (proportional to prob) - main Migdal MC function
vector<double> Target::rand_migdalE(double ENR)
{
    //select random isotope
    double MT = -1;
    if (isoN>1)
        MT = AMU*rand_isotope(ENR);
    else
        MT = Mt;
    
    double maxEM=0;
    maxEM = EMmaxER( maxE, ENR, MT);

    double Ee,EeMax;

    double t,prob=0;
    double maxProb=totalMigProb(maxE, ENR, MT)/maxCumulativeProb;

    t = gsl_rng_uniform(r);
    if ( t < maxProb )
    {
        for(int orbit_i=0; orbit_i<Norbits; orbit_i++)
        {
            EeMax = maxEM - orbits_energy[orbit_i];
           
            if (EeMax < 0)
                continue;
            prob+=Z_nl_integrated(orbit_i,Ve(ENR,MT),EeMax)/maxCumulativeProb;
            if( t < prob )
            {
                t = gsl_rng_uniform(r)*Z_nl_integrated(orbit_i,Ve(ENR,MT),EeMax);
                Ee = invZ_nl_integrated(orbit_i, Ve(ENR,MT), t);
                return {Ee,orbits_energy[orbit_i],orbits[orbit_i][0]-'0',MT/AMU};
            }
        }
    }
    return {0,0,0,0};

}


//returns a randomly selected isotope consistent with ENR
double Target::rand_isotope(double ENR)
{
    double u = gsl_rng_uniform(r);
    double maxER=0;
    for (int i=0;i<isoN;i++)
    {   
        if (u < isoFracTable[i])
        {
            maxER = ERmax(maxE, 0, AMU*isoTable[i]);
            //cout << u << " " << maxER << " " << ENR << " " << AMU*isoTable[i] << endl;
            if( ENR < maxER ) //ensures kinematic consistency
                return isoTable[i];
            else
                continue;
        }
    }
    return Mt; //in case something goes wrong
}


/*
//the following 3 functions are for numerical calculation of the Migdal spectrum
double Target::migdalIntegrand(double ErNR, void *pars)
{
    Spectra *NRspec = (Spectra *) pars;
    double Edet = NRspec->Edet;
    double integrand=0;
    double dRnr;
    int N = NRspec->Natom;

    if( ErNR > ERmax(maxE, Edet - ErNR*Leff))
        return 0;
    
    dRnr = gsl_spline_eval(NRspec->spectrumSpline, ErNR, NRspec->accelSS);
    
    for( int orbit_i=0; orbit_i<Norbits; orbit_i++ )
    {
        if ( Edet - ErNR*Leff - nl_energy[orbit_i] > 0 )
            integrand += Z_nl( orbit_i, Ve(ErNR,Mt), Edet - ErNR*Leff - nl_energy[orbit_i]);
    }
    return integrand * dRnr;
}

typedef double func(double , void*);

//differential Migdal spectrum
double Target::dRdEmigdalNeutron(double Edet, Spectra *NRspec)
{
    cout << maxE;
    gsl_function F;
    F.function = (func*)&this->migdalIntegrand;
    F.params = NRspec;
    NRspec->Edet = Edet;
    int limit = 1000;
	double integral,absErr,tol;
    tol=1e-6;

    W = gsl_integration_workspace_alloc (1000);

    double ErMin = NRspec->xMin; //This is effectively zero for all regions we care about
    double ErMax = ERmax(maxE, Edet);
    if(ErMax > Edet/Leff)
        ErMax = Edet/Leff;
    if(ErMax > NRspec->xMax)
        ErMax = NRspec->xMax;
    
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
void Target::calcMigdalSpectrum(Spectra NRspec)
{
    //TestSpectra::SPLINE_migdal_prep migdalSpec;
    //migdalSpec.NR_spec = NRspec;
    //migdalSpec.migdalSpline = gsl_spline_alloc(gsl_interp_linear,2000);
    //migdalSpec.accelMS = gsl_interp_accel_alloc();
    //migdalSpec.monoE = NRspec->monoE;
    double Er[2000];
    double Rate[2000];
    
    std::fstream RFF;
    for (int N=3; N<6; N++)
    {
        NRspec.Natom = N;

        string filename = NRspec.filename;
        string ins("mig"+to_string(N));
        //filename.erase( filename.find("s"),8);
        filename.replace( filename.find("NR"),2,ins);
        RFF.open(filename,fstream::out);

        double e_min[6] = {-1,35,4.9,0.66,.061,.0098};
        double interval = pow(5/e_min[N],0.005);
        Er[0]=e_min[N];
        for (int i=0; i<200; i++)
        {
            Er[i] = interval*Er[(i-1>=0 ? i-1 : 0)];
            Rate[i] = dRdEmigdalNeutron( Er[i], &NRspec);
            RFF << Er[i] << " " << Rate[i] << endl;
        }
        RFF.close();

        //migdalSpec.xMin = Er[0];
        //migdalSpec.xMax = Er[1999];

        //gsl_spline_init( migdalSpec.migdalSpline, Er, Rate, 2000);
        //migdalSpec.totRate=gsl_spline_eval_integ(migdalSpec.migdalSpline, migdalSpec.xMin, migdalSpec.xMax, migdalSpec.accelMS);
    }

}
*/

