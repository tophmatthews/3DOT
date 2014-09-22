#include <cmath>
#include <stdio.h>
#include <iostream>

#include "ion.h"

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

void ionBase::reset()
{
  pass = 0.0;
  punch = 0.0;
  escapee = false;
  Ehit = 0.0;
  Eout = 0.0;
  travel = 0.0;
  
  if ( type == FG)
    ef = simconf->fg_min_e;
  else
    ef = simconf->ion_min_e;
  switch (type)
  {
    case FF:
      pot = simconf->pot_ff;
      break;
    case LAT:
      pot = simconf->pot_lat;
      break;
    case FG:
      pot = simconf->pot_fg;
      break;
  }
}

void ionBase::parent( ionBase *parent )
{
  gen = parent->gen + 1;
  t = parent->t;
  
  famtree = parent->famtree; // inhereits family tree
  fam_parent = parent->z1;
  famtree.push_back( fam_parent ); // add parent to family tree
  
//  if (parent->gen > 0)
//    hit_e = parent->hit_e;

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
}

ionBase* ionBase::spawnRecoil()
{
  ionBase *recoil = new ionBase;
  recoil->parent(this);
  return recoil;
}

void ionBase::assignType()
{
  if (z1 >=90 || z1 <= 10)
    type = LAT;
  else
    type = FG;
  
  reset();
}

double ionBase::RangeInFuel( std::string fueltype )
{
  double c[5];
  if( fueltype == "uc" )
  {
    if( z1 == 54 ) // if xenon
    {
      c[0]=-0.00567; c[1]=0.10018; c[2]=-0.55072; c[3]=1.57743; c[4]=-0.74034;
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

void ionBase::prep_FF()
{
  gen = 0.0;
  fam_fuel = 0;
  fam_fg = 0;
  tag = -1; // -1 = born in fuel
  type = FF;
  pot = simconf->pot_ff;
}

