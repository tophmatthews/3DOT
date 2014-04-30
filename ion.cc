#include <cmath>

#include "ion.h"

ionBase::ionBase()
{
  t = 0.0; //clock
  reset();
}

ionBase::ionBase( ionBase* prototype )
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
  Eend = 0;
}

void ionBase::parent( ionBase *parent )
{
  reset();
  gen = parent->gen + 1;
  t = parent->t;
  famtree = parent->famtree; // inhereits family tree

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
