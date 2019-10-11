
// Verbosity flag (for limiting output to S1/S2; nothing else
bool verbosity = false;

// controls whether quanta is output or just S1/S2
bool outputQuanta = false;

// General parameters of importance changing global behavior
bool MCtruthE = false;    // false means reconstructed energy
bool MCtruthPos = true;  // false means reconstructed position

int useTiming = 0;  // calculates S2 width, should rename this

int usePD = 1; // 0 means PE, 1 means phd (PE/~1.2), 2 means spike count
int useS2 = 1;  // xtra feature: 2 means S2 x-axis energy scale
int useCorrected = 1;  //print corrected S1/S2 values

double minS1 = 1.;  // units are controlled by the usePD flag
// this is separate from S1 thresholds controlled by detector
double maxS1 = 1e3;
int numBinsS1 = 100;

// minS2 need not match S2 threshold in detector.hh
// you can treat as trigger vs. analysis thresholds
double minS2 = 100.0;
double maxS2 = 1e9;
int numBinsS2 = 100;
int logS2 = 0;

//spatial bins as a fraction of size defined in detector file
double minR=0;
double maxR=1;
int numBinsR=3;
double minZ=0;
double maxZ=1;
int numBinsZ=5;

// some numbers for fine-tuning the speed vs. the accuracy
double z_step = 100.0;  // mm, for integrating non-uniform field
double E_step = 2.0;  // keV, for integrating WIMP spectrum
