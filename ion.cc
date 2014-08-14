#include <cmath>
#include <stdio.h>
#include <iostream>

#include "ion.h"
#include "invert.h"
#include "r250.h"
#include "functions.h"

ionBase::ionBase()
{
  t = 0.0; //clock
  reset();
}

ionBase::ionBase( ionBase* prototype ) // copy constructor
{ 
  ef = prototype->ef; // final energy
  t = prototype->t;   // clock
  z1 = prototype->z1; // atomic number
  m1 = prototype->m1; // mass
  e = prototype->e;   // energy
  reset();
}

ionBase::ionBase( int _z1, double _m1, double _e ) : z1(_z1), m1(_m1), e(_e)
{
  t = 0.0; //clock
  reset();
}

void ionBase::set_ef()
{
  ef = fmax( 5.0, 0.00001 * e ); // final energy -> either 5 or 0.00001 e
}

void ionBase::reset()
{
  pass = 0;
  punch = 0;
  escapee = false;
  ef = 5.0;
  Ehit = 0;
  Eout = 0;
  travel = 0;
}

void ionBase::parent( ionBase *parent )
{
  gen = parent->gen + 1;
  t = parent->t;
  
  famtree = parent->famtree; // inhereits family tree
  fam_parent = parent->z1;
  famtree.push_back( fam_parent ); // add parent to family tree

  fam_fuel = parent->fam_fuel;
  fam_fg = parent->fam_fg;
  
  if( fam_parent <= 90 && fam_parent >= 10 )
  {
    if( gen > 1 )
      ++fam_fg;
  }
  else
    ++fam_fuel;
  
  for( int i = 0; i < 3; i++ )
  {
    pos[i] = parent->pos[i];
    pos0[i] = pos[i];
  }
  
  if (z1 <=90 && z1 >= 10)
    type = LAT;
  else
    type = FG;
}

ionBase* ionBase::spawnRecoil()
{
  ionBase *recoil = new ionBase;
  recoil->parent(this);
  return recoil;
}

double ionBase::RangeInFuel( std::string fueltype )
{
  double c[5];
  if( fueltype == "uc" )
  {
    if( z1 == 54 ) // if xenon
    {
      c[0]=-0.00462; c[1]=0.08158; c[2]=-0.4515; c[3]=1.47142; c[4]=-0.96243;
    }
    else if( z1 == 55 ) // if cesium
    {
      c[0]=-0.00631; c[1]=0.10921; c[2]=-0.5889; c[3]=1.69355; c[4]=-0.98078;
    }
    else if( z1 == 36 ) // if krypton
      c[0]=-0.00567; c[1]=0.09466; c[2]=-0.48892; c[3]=1.48024; c[4]=-0.85407;
  }
  else if( fueltype == "um" )
  {
    if( z1 == 54 ) // if xenon
    {
      c[0]=-0.00412; c[1]=0.0705; c[2]=-0.36317; c[3]=1.18349; c[4]=-0.74233;
    }
    else if( z1 == 55 ) // if cesium
    {
      c[0]=-0.00587; c[1]=0.09994; c[2]=-0.51886; c[3]=1.4794; c[4]=-0.87568;
    }
    else if( z1 == 36 ) // if krypton
    {
      c[0]=-0.00513; c[1]=0.08295; c[2]=-0.39975; c[3]=1.20048; c[4]=-0.65063;
    }
  }
  else if( fueltype == "un" )
  {
    if( z1 == 54 ) // if xenon
    {
      c[0]=-0.00454; c[1]=0.07926; c[2]=-0.42869; c[3]=1.38063; c[4]=-0.84672;
    }
    else if( z1 == 55 ) // if cesium
    {
      c[0]=-0.00634; c[1]=0.10921; c[2]=-0.58467; c[3]=1.66738; c[4]=-0.95908;
    }
    else if( z1 == 36 ) // if krypton
    {
      c[0]=-0.00569; c[1]=0.09428; c[2]=-0.48186; c[3]=1.44725; c[4]=-0.8276;
    }
  }
  else
  {
    fprintf( stderr, "Invalid fuel type specified");
  }
  
  double leng = log10( e );
  double lrange = c[4];
  for( int i=0; i<4; ++i)
    lrange += c[i] * pow( leng, 4-i );
  return pow(10, lrange);
}

void ionBase::prep_FF(ionBase *ff)
{
  ff->gen = 0;
  ff->fam_fuel = 0;
  ff->fam_fg = 0;
  ff->tag = -1; // -1 = born in fuel
  ff->type = FF;
}

void ionBase::make_FF( std::queue<ionBase*> &recoils, int fsn_num )
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
  ff1->prep_FF(ff1);
  ff1->z1 = Z1;
  ff1->m1 = A1;
  ff1->e  = E1 * 1.0e6; // Change energy units to eV
  ff1->ionId = simconf->ionId++;
  
  // Random direction
  do
  {
    for( int i = 0; i < 3; ++i ) ff1->dir[i] = dr250() - 0.5;
    norm = v_dot( ff1->dir, ff1->dir );
  } while( norm <= 0.0001 );
  v_scale( ff1->dir, 1.0 / sqrtf( norm ) );
  
  // random origin
  for( int i = 0; i < 3; ++i )
  {
    ff1->pos[i] = dr250() * simconf->length;
    ff1->pos0[i] = ff1->pos[i]; // set orinal position
  }
  recoils.push( ff1 );
  
  // -- Spawn 2nd fission fragments -- //
  ff2 = new ionBase( *ff1 ); // copy constructor
  ff2->prep_FF(ff2);
  for( int i = 0; i < 3; ++i )
    ff2->dir[i] *= -1.0; // reverse direction
  ff2->z1 = Z2;
  ff2->m1 = A2;
  ff2->e  = E2 * 1.0e6;
  ff2->ionId = simconf->ionId++;
  recoils.push( ff2 );
  
  printf( "%s Fsn %i/%.0f Z1=%d (%.2f MeV)\t Z2=%d (%.2f MeV)\n", simconf->run_name.c_str(), fsn_num, simconf->fissions, Z1, E1, Z2, E2 );
}
