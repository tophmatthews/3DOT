#include "sample.h"
#include "material.h"
#include "element.h"

#include <iostream>

sampleBase::sampleBase( double x, double y, double z, boundaryCondition b )
{
  w[0] = x;
  w[1] = y;
  w[2] = z;
  
  for( int i = 0; i < 3; i++ )
    bc[i] = simconf->bounds;
}

void sampleBase::averages( const ionBase *pka )
{
  for( int i = 0; i < material.size(); i++ )
    material[i]->average( pka );
}

void sampleBase::make_fuel( std::string fueltype, sampleBase *sample, double frac_den )
{
  materialBase *material;
  elementBase *element;

  if( fueltype == "uc" ) // uranium carbide
  {
    material = new materialBase( 13.63 * frac_den ); // density in [g/cc]
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
    material = new materialBase( 14.33 * frac_den ); // density in [g/cc]
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
    material = new materialBase( 19.0 * frac_den ); // density in [g/cc]
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else if( fueltype == "uo2" ) // uranium nitride
  {
    material = new materialBase( 10.97 * frac_den ); // density in [g/cc]
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
    element = new elementBase;
    element->z = 8;
    element->m = 16.0;
    element->t = 2.0;
    material->element.push_back( element );
  }
  else
  {
    fprintf( stderr, "Invalid fuel type specified");
  }
  
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample
}

void sampleBase::make_fg( sampleBase *sample, double bub_den, bool xe_only )
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
  }
  material->prepare();
  sample->material.push_back( material );
}

void sampleBase::make_FF( std::queue<ionBase*> &recoils, int fsn_num )
{
  massInverter *m = new massInverter;
  energyInverter *e = new energyInverter;
  ionBase *ff1, *ff2;
  
  double norm;
  double A1, A2, Etot, E1, E2; // inputs for ff creation
  int Z1, Z2;
  
  // generate fission fragment data
  A1 = m->x( dr250() ); // Randomize first mass from double hump probability
  A2 = 235.0 - A1;
  e->setMass(A1);
  Etot = e->x( dr250() ); // This E is in MeV
  E1 = Etot * A2 / ( A1 + A2 );
  E2 = Etot - E1;
  Z1 = round( ( A1 * 92.0 ) / 235.0 );
  Z2 = 92 - Z1;
  
  // -- Spawn 1st fission fragment -- //
  ff1 = new ionBase;
  ff1->prep_FF();
  ff1->z1 = Z1;
  ff1->m1 = A1;
  ff1->e  = E1 * 1.0e6; // Change energy units to eV
  ff1->ionId = simconf->ionId++;
  
  // Random direction
  do
  {
    for( int i = 0; i < 3; ++i )
      ff1->dir[i] = dr250() - 0.5;
    norm = v_dot( ff1->dir, ff1->dir );
  } while( norm <= 0.0001 );
  v_scale( ff1->dir, 1.0 / sqrtf( norm ) );
  
  // random origin
  double R;
  do {
    R = 0;
    for( int i = 0; i < 3; ++i )
    {
      ff1->pos[i] = dr250() * simconf->length;
      R += sqr( simconf->length/2.0 - ff1->pos[i] );
    }
    R = sqrt(R);
  } while (R < simconf->bub_rad);
  
  for( int i = 0; i < 3; ++i )
    ff1->pos0[i] = ff1->pos[i]; // set orinal position
  recoils.push( ff1 );
  
  // -- Spawn 2nd fission fragments -- //
  ff2 = new ionBase( *ff1 ); // copy constructor
  ff2->prep_FF();
  for( int i = 0; i < 3; ++i )
    ff2->dir[i] *= -1.0; // reverse direction
  ff2->z1 = Z2;
  ff2->m1 = A2;
  ff2->e  = E2 * 1.0e6;
  ff2->ionId = simconf->ionId++;
  recoils.push( ff2 );
  
  printf( "%s Fsn %i/%.0f Z1=%d (%.2f MeV)\t Z2=%d (%.2f MeV)\n", simconf->run_name.c_str(), fsn_num, simconf->fissions, Z1, E1, Z2, E2 );
}