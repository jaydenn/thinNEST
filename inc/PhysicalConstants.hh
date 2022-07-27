#include <vector>
#include <string>
using std::vector;

#ifndef PHYSICALCONSTANTS_HH
#define PHYSICALCONSTANTS_HH

/////////////////////////////////////////////////
//     Define some physical constants          //
/////////////////////////////////////////////////
const double Mn = 939565.42;       //neutron mass in keV
const double AMU = Mn/1.00727;     //1 atomic mass unit in keV
const double Me = 510.999;         //electron mass

const double wimpVMAX = 0.00266; //vesc = 544km/s, v0 = 238km/s, vEavg = 16.7km/s

//Atomic numbers for targets
const int XENON = 54;
const int ARGON = 18;

//Migdal type
#define NEUTRON 1
#define NEUTRINO 2
#define WIMPmig 3

//mass numbers for targets
extern double MtXe;    //RAM of xenon in keV
extern double MtAr;    //RAM of argon in keV

//isotope data
const int isoNXe=7;
const vector<double> isoTableXe = {128,129,130,131,132,134,136}; //mass numbers of isotopes
const vector<double> isoFracTableXe = {0.0191421,0.283724,0.324514,0.536981,0.806574,0.911205,1}; //Cumulative abundance of naturally occuring isotopes of xenon

const int isoNAr=1;
const vector<double> isoTableAr = {40}; //atomic numbers of isotopes
const vector<double> isoFracTableAr = {1}; //naturally occuring isotopic abundances of argon

//atomic binding energies
const int NorbitsXe = 17;
const vector<std::string> orbitsXe = {"1s","2s","2p-","2p","3s","3p-","3p","3d-","3d","4s","4p-","4p","4d-","4d","5s","5p-","5p"};
const vector<double> orbits_energyXe = {34.7559, 5.50935, 5.16145, 4.83559, 1.17037, 1.02478, 0.961249, \
                                        0.708132, 0.6949, 0.22939, 0.175581, 0.1628, 0.0737791, 0.0716683, \
                                        0.0274873, 0.0134036, 0.0119677};

const int NorbitsAr = 2; 
const vector<std::string> orbitsAr = {"1s","2s"};
const vector<double> orbits_energyAr = {1,1}; //the atomic binding energies of shells 


extern double LeffXe;          //approximate quenching factors
extern double LeffAr;


/////////////////////////////////////////////////

#endif /* PHYSICALCONSTANTS_HH */

