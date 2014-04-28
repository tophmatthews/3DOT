#include "sample.h"
#include <iostream>
sampleBase::sampleBase( double x, double y, double z, sampleBoundary b )
{
  w[0] = x;
  w[1] = y;
  w[2] = z;
  
  for( int i = 0; i < 3; i++ )
    bc[i] = b;
}

void sampleBase::averages( const ionBase *pka )
{
  for( int i = 0; i < material.size(); i++ )
    material[i]->average( pka );
}
