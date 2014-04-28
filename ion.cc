#include <cmath>

#include "ion.h"

ionBase::ionBase()
{ 
  ef = 5.0; // final energy
  t = 0.0; //clock
  pass = 0;
  punch = 0;
  escapee = false;
}

ionBase::ionBase( ionBase* prototype )
{ 
  ef = prototype->ef; // final energy
  t = prototype->t;   // clock
  z1 = prototype->z1; // atomic number
  m1 = prototype->m1; // mass
  e = prototype->e;   // energy
  pass = 0;
  punch = 0;
}

ionBase::ionBase( int _z1, double _m1, double _e ) : z1(_z1), m1(_m1), e(_e)
{
  ef = 5.0; // final energy
  t = 0.0; //clock
  escapee = false;
}

void ionBase::set_ef()
{
  ef = fmax( 5.0, 0.00001 * e ); // final energy -> either 5 or 0.00001 e
}

void ionBase::parent( ionBase *parent )
{
  ef = 5.0; // final energy

  gen = parent->gen + 1;
  t = parent->t;
  famtree = parent->famtree; // inhereits family tree
  escapee = false;

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
