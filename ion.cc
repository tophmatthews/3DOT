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
  md = 0;
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
