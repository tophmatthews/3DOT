#include "sample.h"
#include "material.h"
#include "element.h"

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

void sampleBase::make_fuel( std::string fueltype, sampleBase *sample, double smear_den )
{
  materialBase *material;
  elementBase *element;

  if( fueltype == "uc" ) // uranium carbide
  {
    material = new materialBase( 13.63 * smear_den ); // rho
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
    element = new elementBase;
    element->z = 6;
    element->m = 12.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else if( fueltype == "un" ) // uranium nitride
  {
    material = new materialBase( 14.33 * smear_den ); // rho
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
    element = new elementBase;
    element->z = 7;
    element->m = 14.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else if( fueltype == "um" ) // uranium metal
  {
    material = new materialBase( 19.0 * smear_den ); // rho
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else
  {
    fprintf( stderr, "Invalid fuel type specified");
  }
  
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample
}

void sampleBase::make_fg( std::string fueltype, sampleBase *sample, double bub_den, bool xe_only )
{
  materialBase *material;
  elementBase *element;
  
  if( xe_only )
  {
    material = new materialBase( bub_den );
    element = new elementBase;
    element->z = 54;
    element->m = 131.3;
    element->t = 1.0;
    material->element.push_back( element );
    material->prepare();
    sample->material.push_back( material );
  }
  else
  {
    // relative concentrations from Matzke sci. adv. lmfbr fuels, pg 433
    
    // Add Xe
    material = new materialBase( bub_den );
    element = new elementBase;
    element->z = 54;
    element->m = 131.3;
    element->t = 22.8;
    material->element.push_back( element );
    
    // Add Kr
    element = new elementBase;
    element->z = 36;
    element->m = 83.8;
    element->t = 2.5;
    material->element.push_back( element );
    
    // Add cesium
    element = new elementBase;
    element->z = 55;
    element->m = 132.9;
    element->t = 19.2;
    material->element.push_back( element );
    
    material->prepare();
    sample->material.push_back( material );
  }
}