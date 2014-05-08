#include <cmath>

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

void ionBase::prep_FF()
{
  gen = 0;
  fam_fuel = 0;
  fam_fg = 0;
  tag = -1; // -1 = born in fuel
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
    std::cout << "Warning: All fuel acts like uc for AddAndKill" << std::endl;
    if( z1 == 54 ) // if xenon
    {
      c[0]=-0.00462; c[1]=0.08158; c[2]=-0.4515; c[3]=1.47142; c[4]=-0.96243;
    }
    else if( z1 == 55 ) // if cesium
    {
      c[0]=-0.00567; c[1]=0.09466; c[2]=-0.48892; c[3]=1.48024; c[4]=-0.85407;
    }
    else if( z1 == 36 ) // if krypton
    {
      c[0]=-0.00567; c[1]=0.09466; c[2]=-0.48892; c[3]=1.48024; c[4]=-0.85407;
    }
  }
  else if( fueltype == "un" )
  {
    std::cout << "Warning: All fuel acts like uc for AddAndKill" << std::endl;
    if( z1 == 54 ) // if xenon
    {
      c[0]=-0.00462; c[1]=0.08158; c[2]=-0.4515; c[3]=1.47142; c[4]=-0.96243;
    }
    else if( z1 == 55 ) // if cesium
    {
      c[0]=-0.00567; c[1]=0.09466; c[2]=-0.48892; c[3]=1.48024; c[4]=-0.85407;
    }
    else if( z1 == 36 ) // if krypton
    {
      c[0]=-0.00567; c[1]=0.09466; c[2]=-0.48892; c[3]=1.48024; c[4]=-0.85407;
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
