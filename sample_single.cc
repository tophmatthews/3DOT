#include <math.h>

#include "sample_single.h"

//create object, inheret from sampleBase, xyz = w[3] = size of simulation
sampleSingle::sampleSingle( double x, double y, double z , boundaryCondition b)  : sampleBase( x, y, z, b)
{
}

// look if we are within dr of a cluster
// dr == 0.0 means looking if we are inside the cluster
materialBase* sampleSingle::lookupMaterial( double* pos ) 
{
  double R; // distance from center
  
  for (int i=0; i<3; ++i)
    R += sqr( ( w[i]/2 ) - pos[i] );
  R = sqrt(R);
  
  if ( R < simconf->bub_rad )
    return material[1];
  else
    return material[0];
}
