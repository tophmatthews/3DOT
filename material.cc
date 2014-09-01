#include "material.h"
#include "simconf.h"

#include "functions.h"

#include <cmath>
#include <iostream>
using namespace std;

void materialBase::prepare()
{
  double tt = 0.0;

  // get total stoichiometry
  for( int i = 0; i < element.size(); i++ ) 
  {
    if( element[i]->t < 0.0 )
      element[i]->t = 0.0;
    tt += element[i]->t;
  }

  // turns stochiometry into fraction of total
  for( int i = 0; i < element.size(); i++ )
  {
    element[i]->t /= tt;
  }

  am = 0.0; // average mass of elements in material
  az = 0.0; // average atomic number of elements in material
  ar = 0.0; // average atomic radius of elements in material
  maxr = 0.0; // largest atomic radius of elements in material
  for( int i = 0; i < element.size(); i++ ) 
  {
    am += element[i]->m * element[i]->t;
    az += double( element[i]->z ) * element[i]->t;
    ar += simconf->scoef[element[i]->z-1].radius;
    maxr = (simconf->scoef[element[i]->z-1].radius > maxr) ? simconf->scoef[element[i]->z-1].radius : maxr;
  }
  
  // rho is bulk density in g/cc: [atoms/A] = [g/cc][atoms/mol][cc/A^3] / [avg g/mol]
  arho = rho * 0.6022 / am; // atomic density of material [atoms/Ang^3]
  
  double minls =cbrt( 1 / arho);
}

// make sure layers are prepare'd first!
void materialBase::average( const ionBase *pka )
{
  mu = pka->m1 / am; //

  double a0 = 0.5292; // [A] bohr radius in angstroms
  double e2 = 14.4; // [ev*A] squared charge of electron
  gamma = ( 4.0 * mu ) / sqr( 1.0 + mu ); // gamma using average material mass
  
  // universal stopping distance for material (old TRIM eq. 2-60)
  a = 0.8854 * a0 / ( pow( double(pka->z1), 0.23 ) + pow( az, 0.23 ) );

  // Center of mass average mass
  // Mc = pka->m1 * am / ( pka->m1 + am );
  // mean flight path
  // f = eps/E_0 (new TRIM eq 7-5): Ec = (E0 M_tar) / (Mpka + M_tar) (new trim 2-78)
  f = a * am / ( double(pka->z1) * az * e2 * ( pka->m1 + am ) );
  //eps0 = e0 * f;
  
  //setup for pmax calc: reduces to psi^2 / eps (newtrim eq 7-41
  epsdg = simconf->tmin * f / gamma; // reduces to epsilon*tmin/tmax

  // fd and kd determine how much recoil energy goes into el. loss and vaccancies
  fd = pow( 0.01 * az, -7.0 / 3.0 );
  kd = pow( 0.1334 * az, 2.0 / 3.0 ) / sqrtf( am ); 

  for( int i = 0; i < element.size(); i++ ) 
  {
    element[i]->my = pka->m1 / element[i]->m;
    element[i]->ec = 4.0 * element[i]->my / sqr( 1.0 + element[i]->my); //gamma
    
    // universal stopping distance
    element[i]->ai = a0 * .8854 / ( pow( double(pka->z1), 0.23 ) + pow( element[i]->z, 0.23 ) );
    
    // f = eps/E (new TRIM eq 7-5)
    element[i]->fi = element[i]->ai * element[i]->m / 
                     ( double(pka->z1) * double(element[i]->z) * e2 * ( pka->m1 + element[i]->m ) );
  }  
  dirty = false;
}

// make sure layers are prepare'd and averaged first!
double materialBase::getrstop( const ionBase *pka )
{
  double se = 0.0;
  for( int i = 0; i < element.size(); i++ ) 
    se += rstop( pka, element[i]->z ) * element[i]->t * arho;

  return se;
}

// proton electric stopping
double materialBase::rpstop( int z2p, double e )
{
  double pe, pe0, sl, sh, sp, velpwr;
  int z2 = z2p-1;
  // velocity proportional stopping below pe0
  pe0 = 25.0;
  pe = fmax( pe0, e );

  //for( int i=0;i<8;i++) fprintf( stderr, "%f ", pcoef[z2][i] );
  //fprintf( stderr, "\n" );

  // pcoef indices are one less than in the fortran version!
  sl = ( simconf->pcoef[z2][0] * pow( pe, simconf->pcoef[z2][1] ) ) +
       ( simconf->pcoef[z2][2] * pow( pe, simconf->pcoef[z2][3] ) );
  sh = simconf->pcoef[z2][4] / pow( pe, simconf->pcoef[z2][5] ) * 
       logf( simconf->pcoef[z2][6] / pe + simconf->pcoef[z2][7] * pe );
  sp = sl * sh / (sl + sh );
  if( e <= pe0 )
  {
    // velpwr is the power of velocity stopping below pe0
    if( z2p <= 6 )
      velpwr = 0.25;
    else
      velpwr = 0.45;
    sp *= pow( e/pe0, velpwr );
  }
  return sp;
}

//electronic stopping
double materialBase::rstop( const ionBase *ion, int z2 )
{
  double e, vrmin, yrmin, v, vr, yr, vmin, m1;
  double a, b, q, q1, l, l0, l1;
  double zeta;
  int z1 = ion->z1;
  double fz1 = double(z1), fz2 = double(z2);
  double eee, sp, power;
  double se;
  
  double v2_b = 50.0; // square of Bohr velocity [keV/u]
  
  // scoeff
  double lfctr = simconf->scoef[z1-1].lfctr;
  double mm1 = simconf->scoef[z1-1].mm1; // atomic mass # protons+neutrons
  double vfermi = simconf->scoef[z2-1].vfermi; // really v_f / v_0 -> fermi velocity over bohr velocity
  double atrho = simconf->scoef[z2-1].atrho; // atomic density
  //fprintf( stderr, "lfctr=%f mm1=%f vfermi=%f atrho=%e\n", lfctr, mm1, vfermi, atrho );

  if( ion->m1 == 0.0 ) 
    m1 = mm1;
  else
    m1 = ion->m1;

  e = 0.001 * ion->e / m1; // energy/mass [keV/amu]  ion->e = [eV], m1 = amu
  if( z1 == 1 ) 
  {
    fprintf( stderr, "proton stopping not yet implemented" );
  }
  else if( z1 == 2 )
  {
    fprintf( stderr, "alpha stopping not yet implemented" );
  }
  else
  {
    yrmin = 0.13;
    vrmin = 1.0;
    v = sqrtf( 2 * e / v2_b) / vfermi; // this gives v_ion/(vfermi * v_b) = v_ion/v_f (since vfermi = v/v_0

    if( v >= 1.0 )
      vr = v * vfermi * ( 1.0 + 1.0 / ( 5.0 * v*v ) ); // [SRIM eq 3-34] gives v_rel / v_0
    else
      vr = ( 3.0 * vfermi / 4.0 ) * ( 1.0 + ( 2.0 * v*v / 3.0 ) - pow( v, 4.0 ) / 15.0 ); //[SRIM eq 3-35] gives v_rel / v_0

    yr = fmax( yrmin, vr / pow(fz1,0.6667) ); // [SRIM eq 3-39] note vr = v_rel / v_0
    yr = fmax( yr, vrmin / pow(fz1,0.6667) );
    a = -0.803 * pow( yr, 0.3 ) + 1.3167 * pow( yr, 0.6 ) + 0.38157 * yr +  0.008983 * yr*yr; // [SRIM eq. 3-41] q = 1 + exp(-a)

    // ionization level of the ion at velocity yr
    q = fmin( 1.0, fmax( 0.0, 1.0 - exp( -fmin( a, 50.0 ) ) ) ); // [SRIM eq. 3-41] with protection for too low an a, and 0<q<1

    b = ( fmin( 0.43, fmax( 0.32, 0.12 + 0.025 * fz1 ) ) ) / pow( fz1, 0.3333 );
    l0 = ( 0.8 - q * fmin( 1.2, 0.6 + fz1 / 30.0) ) / pow( fz1, 0.3333 );
    if( q < 0.2 ) 
      l1 = 0.0;
    else if( q < fmax( 0.0, 0.9 - 0.025 * fz1 ) ) 
    {//210
      q1 = 0.2;
      l1 = b * ( q - 0.2 ) / fabs( fmax( 0.0, 0.9 - 0.025 * fz1 ) - 0.2000001 );
    }
    else if( q < fmax( 0.0, 1.0 - 0.025 * fmin( 16.0, fz1 ) ) )
      l1 = b;
    else
      l1 = b * ( 1.0 - q ) / ( 0.025 * fmin( 16.0, fz1 ) );

    l = fmax( l1, l0 * lfctr ); // l is capital Gamma kinda
    
    // zeta is gamma kinda, not sure what term in last paraentheses is [SRIM eq. 3-33]
    // note vfermi = v_f / v_0
    zeta = q + ( 1.0 / ( 2.0 * vfermi*vfermi ) ) * ( 1.0 - q ) * logf( 1.0 + sqr( 4.0 * l * vfermi / 1.919  ) ); //kf = 1.919 from (Brandt 1982)

    // add z1^3 effect
    a = -sqr( 7.6 - fmax( 0.0, logf( e ) ) ); // a is redefined here!!!
    zeta *= 1.0 + ( 1.0 / (fz1*fz1) ) * ( 0.18 + 0.0015 * fz2 ) * expf( a );

    if( yr <= fmax( yrmin, vrmin / pow( fz1, 0.6667 ) ) ) // missing third test of if yr<= z1^(2/3)
    {
      // calculate velocity stopping for  yr < yrmin
      vrmin = fmax( vrmin, yrmin * pow( fz1, 0.6667 ) ); // vrmin redefined here!!! 1 previously, now 1 or v_rel / v_0
      vmin = 0.5 * ( vrmin + sqrtf( fmax( 0.0, vrmin*vrmin - 0.8 * vfermi*vfermi ) ) );
      eee = v2_b / 2 * vmin*vmin; // kev/amu
      sp = rpstop( z2, eee );

      if( z2 == 6 || ( ( z2 == 14 || z2 == 32 ) && z1 <= 19 ) ) 
        power = 0.375;
      else
        power = 0.5;

      se = sp * sqr( zeta * fz1 ) * pow( e/eee, power );
      //printf("a: se[%d]=%f, %f %f %f %f %f %f\n", i, se[i], e, eee, power, zeta, fz1, sp );
    }
    else
    {
      sp = rpstop( z2, e );
      se = sp * sqr( zeta * fz1 ); //[SRIM eq 3.15]
    }
  } // END: heavy-ions

  return se * 10.0;
}
