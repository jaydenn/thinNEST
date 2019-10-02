#include "NEST.hh"

int trigger(double S2, int trigerIndex)
{
    switch(triggerIndex)
    {
        case 1; //XENON1t
        {
            if(S2<90)
                return 0; 
            if(S2<=142 && RandomGen::rndm()->rand_uniform()<0.5*erfc(0.036855*(142-S2)) )
                return 1;
            else if(RandomGen::rndm()->rand_uniform()<0.5*erfc(0.030066*(142-S2)) )
                return 1;
            else
                return 0;
            break;
        }
    }
    return 1;
}

