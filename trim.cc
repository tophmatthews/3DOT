#include <math.h>
#include <stdio.h>

#include "trim.h"
#include "simconf.h"

#include "functions.h"
#include <iostream>
using namespace std;

//#define RANGECORRECT2

//
// all energies are in eV
//

// does a single ion cascade
void trimBase::trim( ionBase *pka_, queue<ionBase*> &recoils)
{
  // simconf should already be initialized
  pka = pka_;
  
  // make recoil queue available in overloadable functions
  recoil_queue_ptr = &recoils;

  double pl = 0.0;
  double max = 0.0;
  double e0kev = pka->e / 1000.0;
  int ic = 0;
  int nn, ie;
  double r1, r2, hh; // random numbers
  double eps, epsm; // reduced energy for (element, material), epsilon
  double eeg; // t1/2 minimum
  double ls;  // mean free path
  double p, b, r, see, dee;
  double s2, c2, ct, st;
  double rr, ex1, ex2, ex3, ex4, v ,v1;
  double fr, fr1, q, roc, sqe;
  double cc, aa, ff, co, delta;
  double mys2, myff, mydelta;
  double den;
  double rdir[3], perp[3], norm, psi;
  
  double testPos; // projected position. used for pathlength fix

  bool Rangefix; // switch for mfp fix
  double u[2];
  
  double p1, p2;
  double range;
  
  double wrap; //fix for periodic boundary conditions
  bool wrapped;

  // Constants for universal potentail, from [SRIM table 4-1]
  double u_c1 = 0.99229;
  double u_c2 = 0.011615;
  double u_c3 = 0.007122;
  double u_c5 = 9.3066;
  double u_c4 = 14.813;
  
  r1 = dr250();
  
  double bubcen[3];
  for( int i = 0; i < 3; i++ ) bubcen[i] = sample->w[i] / 2;
  
  terminate = false;
  edged = false;

  do // cycle for each collision
  {
    r2 = dr250();
    hh = dr250(); // selects element inside material to scatter from

    ic++; // number of collisions for each ion in each material
    material = sample->lookupMaterial( pka->pos );
    if( material == 0 ) break;

    v_norm( pka->dir ); // normalize direction vector

    // setup avg max. impact parameter
    epsm = pka->e * material->f; // reduced energy, epsilon
    eeg = sqrtf( epsm * material->epsdg ); // [TRI02450] reduces to eps * sqrt(tmin/tmax). t1/2 minimum? eq 4.62
    material->pmax = material->a / ( eeg + sqrtf( eeg ) + 0.125 * pow( eeg, 0.1 ) ); // pmax undefined before here

    ls = 1.0 / ( M_PI * sqr( material->pmax ) * material->arho );
    
    if( simconf->fullTraj )
      printf( "\n\nls: %f\n", ls);

    if( ic == 1 ) // if first recoil of atom
    {
      ls = ls * dr250();
      if( simconf->fullTraj ) printf( "  adjusted ls: %f\n", ls);
    }
    
    // Check to see if line intersects sphere
    // note: this only works for single sphere in center of sample space

    Rangefix = false;
    
    if( simconf->BoundaryFix )
    {
      double aaa = 0, bbb = 0, ccc = 0;
      for( int i = 0; i < 3; i++ )
      {
        testPos = pka->pos[i] + pka->dir[i] * ls;
        aaa += sqr( testPos - pka->pos[i] );
        bbb += 2 * ( testPos - pka->pos[i] ) * ( pka->pos[i] - sample->w[i] / 2 );
        ccc += sqr( sample->w[i]/2 ) + sqr( pka->pos[i] ) - sample->w[i] * pka->pos[i];
      }
      ccc -= sqr( simconf->bub_rad );

      double ddd = bbb*bbb - 4 * aaa * ccc;
      
      if( ddd > 0 )
      {
        u[0] = (- bbb - sqrtf( ddd )) / (2 * aaa);
        u[1] = (- bbb + sqrtf( ddd )) / (2 * aaa);
        if( u[0] >= 0 && u[0] <= 1 )
        {
          Rangefix = true;
          ls *=u[0];
        }
        else if( u[1] >=0 && u[1] <=1 )
        {
          Rangefix = true;
          ls *=u[1];
        }
      }
    
      if( simconf->fullTraj ) printf("u1: %f u2:%f\n", u[0], u[1]);
    }
    
    if(Rangefix) // boundary is crossed
    {
      ic = 0;       // reset collision count
      pka->punch++; // add to punch count
      
      ls += 0.1;                  // add a bit to ensure pka travels across boundary
      pka->e -= ls * material->getrstop( pka ); // electronic energy loss
      if( pka->e < 0.0)
        fprintf( stderr, " electronic energy loss stopped the ion. Broken recoil!!\n" );
      
      if( pka->tag >= 0 )
      {
        double escapeeng = -1 * pka->e * u[0]/abs(u[0]);
        pka->elist.push_back( escapeeng ); // set transition energy. + is leaving, - is entering
        if( escapeeng > 0 )
          pka->escapee = true;
      }

      pl += ls; // add travel length to total path length
      
      for( int i = 0; i < 3; i++ )
        pka->pos[i] += ls * pka->dir[i]; // update position
      
      if( simconf->fullTraj )
        printf( "CROSS at pos: %f %f %f e: %f ls*u: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, ls );
    }
    else // boundary is not crossed
    {
      // choose impact parameter
      p = material->pmax * sqrtf( r2 ); // random p, weighted towards pmax

      // determine which atom in the material will be hit
      for( nn = 0; nn < material->element.size(); nn++ )
      {
        hh -= material->element[nn]->t; // hh is random number
        if( hh <= 0 )
          break;
      }
      element = material->element[nn];

      eps = element->fi * pka->e; // epsilon of given element REDEFINED HERE!!
      b = p / element->ai;        // reduced impact parameter

      see = material->getrstop( pka ); // electronic stopping power
      dee = ls * see; // electronic energy loss: xsec * pathlength
      
      if( eps > 10.0 )
      {
        // use rutherford scattering. includes correction (1+b*(1+b))
        s2 = 1.0 / ( 1.0 + ( 1.0 + b * ( 1.0 + b ) ) * sqr( 2.0 * eps * b ) ); // sin^2(theta_c/2) [SRIM eq. 7-15]
        c2 = 1.0 - s2; // cos^2(theta_c/2)
      }
      
      else
      {
        // first guess at ion c.p.a.
        r = b;
        rr = -2.7 * logf( eps * b );
        if( rr >= b )
        {
          r = rr;
          rr = -2.7 * logf( eps * rr );
          if( rr >= b ) r = rr;
        }
        
        do
        { // Newton's method for r_0: distance of closest approach
          // universal
          // r = reduced radius x, or r/a_u
          ex1 = 0.18175 * exp( -3.1998 * r ); // [SRIM eq. 2-74]
          ex2 = 0.50986 * exp( -0.94229 * r );
          ex3 = 0.28022 * exp( -0.4029 * r );
          ex4 = 0.028171 * exp( -0.20162 * r );
          v = ( ex1 + ex2 + ex3 + ex4 ) / r; // universal screening function / reduced radius
          v1 = -( v + 3.1998 * ex1 + 0.94229 * ex2 + 0.4029 * ex3 + 0.20162 * ex4 ) / r; // -dv/dr

          fr = b*b / r + v * r / eps - r; // similar to [SRIM eq. 7-2] but really denomenator in [SRIM eq. 2-79]
          fr1 = - b*b / ( r*r ) + ( v + v1 * r ) / eps - 1.0; // d(fr)/dr
          q = fr / fr1; // f(x0)/f'(x0)
          r -= q; // Newton's method -> x1 = x0 - f(x0)/f'(x0)
        } while( fabs( q / r ) > 0.001 ); // Convergence criteria

        roc = -2.0 * ( eps - v ) / v1; // R_c [SRIM eq. 7-7], calculated from 7-4
        sqe = sqrtf( eps ); // square root of epsilon

        
        // b = p/a [SRIM eq. 7-7]
        // 5-parameter magic scattering calculation (universal pot.)
        cc = ( u_c2 + sqe ) / ( u_c3 + sqe ); // [SRIM eq. 7-12] beta
        aa = 2.0 * eps * ( 1.0 + ( u_c1 / sqe ) ) * pow( b, cc ); // [SRIM eq. 7-11] A
        
        ff = ( sqrtf( aa*aa + 1.0 ) - aa ) * ( ( u_c4 + eps ) / ( u_c5 + eps ) );
        //ff = 1 / ( sqrtf( aa*aa + 1.0 ) - aa ) * ( ( u_c4 + eps ) / ( u_c5 + eps ) ); // [SRIM eq. 7-11] G

        delta = ( r - b ) * aa * ff / ( ff + 1.0 );
        //delta =( r - b ) * aa / ( ff + 1.0 ); // [SRIM eq. 7-10]
        
        co = ( b + delta + roc ) / ( r + roc ); // MAGIC formula: cosine(theta/2) [SRIM eq. 7-8]
        c2 = co*co; //cosine(theta/2)^2
        s2 = 1.0 - c2; //sin(theta/2)^2 = 1-cos(theta/2)^2
      } // end MAGIC formulation
      
      ct = 2.0 * c2 - 1.0;        // cos(theta_c)
      st = sqrtf( 1.0 - ct*ct );  //sin(theta_c)
      
      // energy transferred to recoil atom
      den = element->ec * s2 * pka->e; // T=gamma*E*sin^2(theta_c/2). gamma*E = Tmax

      // advance clock pathlength/velocity
      pka->t += 10.1811859 * ls / sqrt( 2.0 * pka->e / pka->m1 );
      // time in fs! m in u, l in Ang, e in eV
      // 1000g/kg, 6.022e23/mol, 1.602e-19J/eV, 1e5m/s=1Ang/fs 1.0/0.09822038
      //printf( "se %d  %f [eV]  %f [keV/nm]  %f [nm]\n", pka->ionId, pka->e, see/100.0, pl/10.0 );

      pka->e -= dee; // subtract electronic energy loss
      
      if( pka->e < 0.0 && den > 100.0 ) 
        fprintf( stderr, " electronic energy loss stopped the ion. Broken recoil!!\n" );

      p1 = sqrtf( 2.0 * pka->m1 * pka->e ); // momentum before collision
      pka->e -= den;                        // subtract energy lost in collision
      if( pka->e < 0.0 ) pka->e = 0.0;
      p2 = sqrtf( 2.0 * pka->m1 * pka->e ); // momentum after collision 

      if( dee > max ) max = dee;

      pl += ls; // add travel length to total pka path length

      for( int i = 0; i < 3; i++ ) pka->pos[i] += pka->dir[i] * ls ;
      
      recoil = pka->spawnRecoil();
      
      for( int i = 0; i < 3; i++ ) // progress ion
      {
        // used to assign the new position to the recoil, but
        // we have to make sure the recoil starts in the appropriate material!
        //pka->pos[i] += pka->dir[i] * ( ls - simconf->tau ); // advance pka position
        recoil->dir[i] = pka->dir[i] * p1; // initial pka momentum vector
        
        // Fix for Periodic boundary condtions
        wrapped = false;
        if( sample->bc[i] == sampleBase::CUT && ( pka->pos[i] > sample->w[i] || pka->pos[i] < 0.0 ) )
          terminate = true;
        else if ( sample->bc[i] == sampleBase::PBC )
        {
          if (pka->pos[i] < 0) // fix for negative values
          {
            if( pka->gen != 0 )
            {
              if( simconf->AddAndKill )
                edged = true;// PBC for non-FF: move it and kill it
            }
            else
            {
              pka->pos[i] += sample->w[i];
              wrapped = true;
              if (pka->pos[i] < 0) // protection for too big a box, negative direction
              {
                printf("Failure in Trim.cc: simulation box. PKA wrapped around box in negative direction\n" );
                printf("dim: %i pos: %f w: %f\n", i, pka->pos[i], sample->w[i]);
                exit (EXIT_FAILURE);
              }
            }
          }
          
          if (pka->pos[i] > sample->w[i]) // fix for values past box wall
          {
            if( pka->gen != 0 )
            {
              if( simconf->AddAndKill )
                edged = true;// PBC for non-FF: move it and kill it
            }
            wrap = floor( pka->pos[i] / sample->w[i] );
            if (wrap > 1) // protection for too big a box
            {
              printf("Failure in Trim.cc: simulation box is too small\n" );
              printf("wrap: %f dim: %i pos: %f w: %f\n", wrap, i, pka->pos[i], sample->w[i]);
              exit (EXIT_FAILURE);
            }
            else
            {
              pka->pos[i] -= sample->w[i];
              wrapped = true;
            }
          }
          if( wrapped )
            pka->pass++;
        } // End PBC fix
      }
      
      recoil->e = den;              // assign recoil energy
      recoil->e -= element->Elbind; // subtract lattice binding energy from recoil
      recoil->m1 = element->m;      // assign mass
      recoil->z1 = element->z;      // assign atomic number
        
      // create a random vector perpendicular to pka.dir
      // there is a cleverer way by using the azimuthal angle of scatter...
      do
      {
        for( int i = 0; i < 3; i++ ) rdir[i] = dr250() - 0.5;
        v_cross( pka->dir, rdir, perp );
        norm = sqrtf( v_dot( perp, perp) );
      }
      while( norm == 0.0 );
      v_scale( perp, 1.0 / norm );

      psi = atan( st / ( ct + element->my ) ); // lab scattering angle [SRIM eq. 7.17]
      v_scale( pka->dir, cos( psi ) );

      // calculate new direction, subtract from old dir (stored in recoil)
      for( int i = 0; i < 3; i++ ) 
      {
        pka->dir[i] += perp[i] * sin( psi );
        recoil->dir[i] -= pka->dir[i] * p2;
      }
      
      if( edged )
      {
        if( simconf->fullTraj )
          printf( "EDGED at pos: %f %f %f e: %f ls*u: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, ls );
        
        terminate = true;
        ls = rangeInFuel( pka->e, simconf->fueltype );
        if( simconf->fullTraj )
          printf( "e: %f  ls: %f dir: %f %f %f\n", pka->e, ls, pka->dir[0], pka->dir[1], pka->dir[2]);
        pl += ls; // add travel length to total pka path length
        for( int i = 0; i < 3; i++ )
          pka->pos[i] += pka->dir[i] * ls;
        
        if( simconf->fullTraj )
          printf( "Final pos: %f %f %f\n\n", pka->pos[0], pka->pos[1], pka->pos[2]);
      }
      
      // put the recoil on the stack
      if( spawnRecoilLimit() && !terminate && simconf->makeRecoils )
      {
        v_norm( recoil->dir );
        recoil->tag = material->tag; // mark where new recoil was birthed
        
        if( recoil->tag >= 0 )
          recoil->elist.push_back( recoil->e );

        if( pka->md > 0 ) recoil->md = pka->md +1; // if pka is in range, mark recoil +1
        else recoil->md = 0;
        
        recoil->ionId = simconf->ionId++;           // set new id
        recoil->famtree.push_back( pka->z1 ); // add parent to family tree
        recoils.push( recoil );               // add recoil to stack
        
        if( simconf->fullTraj )
          printf( "SPAWN id: %li pos: %f %f %f e: %f\n", recoil->ionId, recoil->pos[0], recoil->pos[1], recoil->pos[2], recoil->e);
        
        if( simconf->fullTraj )
        {
          double asdf[3];
          for( int i=0; i<3; i++ ) asdf[i] = sample->w[i] / 2 - pka->pos[i];
          double RRR = sqrt( v_dot( asdf, asdf ));
          printf( "HIT pos: %f %f %f fromc: %f pkaz: %d pkae: %f recoilz: %d recoile:%f\n", pka->pos[0], pka->pos[1], pka->pos[2], RRR, pka->z1, pka->e, recoil->z1, recoil->e );
        }
      }
      else delete recoil;
      
    } // endif for Rangefix == false
  } while ( pka->e > pka->ef && !terminate );
  
  if( pka->tag >= 0 )
    pka->elist.push_back( pka->e );
  pka->travel = pl;
  if( simconf->fullTraj && pka->tag >= 0 ) printf( "\n" );
  if( simconf->fullTraj ) printf("particle killed \n\n");
}
