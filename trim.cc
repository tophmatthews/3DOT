#include <math.h>
#include <stdio.h>
#include <iostream>

#include "trim.h"

// does a single ion cascade
void trimBase::trim( ionBase *pka_, queue<ionBase*> &recoils)
{
  // simconf should already be initialized
  pka = pka_;
  
  // make recoil queue available in overloadable functions
  recoil_queue_ptr = &recoils;

  int ic = 0;
  double s2, c2, ct, st;
  double den;
  double rdir[3], perp[3], norm, psi;
  
  double p1, p2;
  double range;
  
  double r1 = dr250();
  
  double bubcen[3];
  for( int i = 0; i < 3; i++ )
    bubcen[i] = sample->w[i] / 2;
  
  terminate = false;
  edged = false;
  
  if (simconf->fullTraj)
  {
    printf( "\n+++=== ION START pos: %.2f %.2f %.2f ===+++\n", pka->pos[0], pka->pos[1], pka->pos[2]);
    printf( "     ion z: %i  e: %f \n", pka->z1, pka->e);
  }
  
  do // cycle for each collision
  {
    ++ic; // number of collisions for each ion in each material
    material = sample->lookupMaterial( pka->pos );
    oldMaterial = material;
    if (simconf->fullTraj)
      cout << "in material with avg z: " << material->az << endl;
    if (material == 0)
      break;

    v_norm( pka->dir ); // normalize direction vector
    
    setPmax( pka, material );
    
    // check changes
    
    double ls = 1.0 / ( M_PI * sqr( material->pmax ) * material->arho ); // (newtrim eq 7-28); // calculate path length
    
//    if (ls > 200 && material->az == 49)
//    {
//      cout << "SOMETHING BAD arho: " << material->arho << endl;
//      cout << "ls: " << ls << endl;
//      cout << "pmax: " << material->pmax << endl;
//      cout << "minls: " << material->minls << endl;
//    }
    
    ls = (ls > material->minls ? ls : material->minls);

    
    
    if (simconf->fullTraj)
      printf( "\nls: %f\n", ls);

    if ( ic == 1 ) // if first recoil of atom
    {
      ls = ls * dr250();
      if (simconf->fullTraj) printf( "  adjusted ls: %f\n", ls);
    }
    
    // Check to see if line intersects sphere
    // note: this only works for single sphere in center of sample space
    
    bool crossed = false;
    if (simconf->BoundaryFix)
      rangeFix( pka, sample, crossed, ls );
    
    if (crossed) // boundary is crossed
    {
      ic = 0;       // reset collision count
      pka->punch++; // add to punch count
      
      ls += 0.1;                  // add a bit to ensure pka travels across boundary
      
      if (simconf->calc_eloss)
        doELoss( pka, material, ls );

      if (pka->type >= FG)
      {
        pka->escapee = true;
        if (pka->Eout == 0)
          pka->Eout = pka->e; // set transition energy. + is leaving, - is entering
      }

      pka->travel += ls; // add travel length to total path length
      
      for (int i = 0; i < 3; ++i)
      {
        pka->pos[i] += ls * pka->dir[i]; // update position
      }
      
      if (simconf->fullTraj)
      {
        printf( "CROSS at pos: %.2f %.2f %.2f e: %f ls*u: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, ls );
        printf( "Leaving material with average Z %0.1f\n", material->az);
      }
    }
    else // boundary is not crossed
    {
      s2 = calcS2( pka, material );
      
      c2 = 1.0 - s2;
      ct = 2.0 * c2 - 1.0;        // cos(theta_c)
      st = sqrtf( 1.0 - ct*ct );  // sin(theta_c)
      
      // energy transferred to recoil atom
      den = element->ec * s2 * pka->e; // T=gamma*E*sin^2(theta_c/2). gamma*E = Tmax

      // advance clock pathlength/velocity
      pka->t += 10.1811859 * ls / sqrt( 2.0 * pka->e / pka->m1 );
      // time in fs! m in u, l in Ang, e in eV
      // 1000g/kg, 6.022e23/mol, 1.602e-19J/eV, 1e5m/s=1Ang/fs 1.0/0.09822038
      //printf( "se %d  %f [eV]  %f [keV/nm]  %f [nm]\n", pka->ionId, pka->e, see/100.0, pl/10.0 );

      if (simconf->calc_eloss)
        doELoss( pka, material, ls );
      
      pka->hit_e.push_back( floor(pka->e) );

      //if (pka->type == FG)
      //printf("%.0f, %.0f;\n", pka->e, den);
      
      p1 = sqrtf( 2.0 * pka->m1 * pka->e ); // momentum before collision
      pka->e -= den;                        // subtract energy lost in collision
      if (pka->e < 0.0)
        pka->e = 0.0;
      p2 = sqrtf( 2.0 * pka->m1 * pka->e ); // momentum after collision

      pka->travel += ls; // add travel length to total pka path length

      for (int i = 0; i < 3; ++i) // move ion to location of recoil
        pka->pos[i] += pka->dir[i] * ls ;
      
      if (simconf->fullTraj)
        printf( "HIT pos: %f %f %f \te: %f \tT: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, den);
      
      for( int i = 0; i < 3; ++i ) // ion location fix based on Boundary conditions
      {
        if (sample->bc[i] == sampleBase::CUT && ( pka->pos[i] > sample->w[i] || pka->pos[i] < 0.0 ))
          terminate = true;
        else if (sample->bc[i] == sampleBase::PBC)
        {
          bool wrapped = false;
          if (pka->pos[i] < 0) // fix for negative values
          {
            if (pka->type == FG && simconf->AddAndKill)
              edged = true;
            else
            {
              pka->pos[i] += sample->w[i];
              wrapped = true;
              if (pka->pos[i] < 0) // protection for too big a box, negative direction
              {
                printf("Failure in Trim.cc: simulation box too small. Ion wrapped around box in negative direction\n" );
                printf("ls: %f dim: %i pos: %f w: %f\n",ls,  i, pka->pos[i], sample->w[i]);
                printf("end pos: %f %f %f\n", pka->pos[0], pka->pos[1], pka->pos[2]);
                printf("original pos: %f %f %f\n", pka->posOld[0], pka->posOld[1], pka->posOld[2]);
                printf("old material az: %0.1f \t new material az: %0.1f\n", oldMaterial->az, material->az);
                exit (EXIT_FAILURE);
              }
            }
          }
          
          else if (pka->pos[i] > sample->w[i]) // fix for values past box wall
          {
            if (pka->type == FG && simconf->AddAndKill)
              edged = true;
            else
            {
              double wrap = floor( pka->pos[i] / sample->w[i] );
              if (wrap > 1) // protection for too big a box
              {
                printf("Failure in Trim.cc: simulation box is too small. Ion wrapped around box in positive direction\n" );
                printf("ls: %f wrap: %f dim: %i pos: %f w: %f\n", ls, wrap, i, pka->pos[i], sample->w[i]);
                printf("end pos: %f %f %f\n", pka->pos[0], pka->pos[1], pka->pos[2]);
                printf("original pos: %f %f %f\n", pka->posOld[0], pka->posOld[1], pka->posOld[2]);
                printf("old material az: %0.1f \t new material az: %0.1f\n", oldMaterial->az, material->az);
                exit (EXIT_FAILURE);
              }
              else
              {
                pka->pos[i] -= sample->w[i];
                wrapped = true;
              }
            }
          } // end positive boundary fix
          if (wrapped) pka->pass++;
        } // end PBC fix
      }
      
      recoil = pka->spawnRecoil();
      recoil->e = den;              // assign recoil energy
      recoil->e -= element->Elbind; // subtract lattice binding energy from recoil
      recoil->m1 = element->m;      // assign mass
      recoil->z1 = element->z;      // assign atomic number
      recoil->assignType();
      
      // create a random vector perpendicular to pka->dir
      // there is a cleverer way by using the azimuthal angle of scatter...
      do
      {
        for (int i = 0; i < 3; ++i)
          rdir[i] = dr250() - 0.5;
        v_cross( pka->dir, rdir, perp ); // takes 1st and 2nd, returns third
        norm = sqrtf( v_dot( perp, perp) );
      }
      while (norm == 0.0);
      v_scale( perp, 1.0 / norm );

      psi = atan( st / ( ct + element->my ) ); // lab scattering angle [SRIM eq. 7.17]
      v_scale( pka->dir, cos( psi ) );

      // calculate new direction, subtract from old dir (stored in recoil)
      for (int i = 0; i < 3; ++i)
      {
        recoil->dir[i] = pka->dir[i] * p1; // initial pka momentum vector
        pka->dir[i] += perp[i] * sin( psi );
        recoil->dir[i] -= pka->dir[i] * p2;
      }
      v_norm( recoil->dir );
      
      if (edged)
      {
        if (simconf->fullTraj)
          printf( "EDGED at pos: %f %f %f e: %f ls*u: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, ls );
        
        terminate = true;
        double est_range = pka->RangeInFuel( simconf->fueltype );
        if (simconf->fullTraj)
          printf( "  e: %f  ls: %f dir: %f %f %f\n", pka->e, ls, pka->dir[0], pka->dir[1], pka->dir[2]);
        
        // add travel length to total pka path length
        // will not work well for total travel length, since this is actually range. But good enough!
        pka->travel += est_range;
        for (int i = 0; i < 3; ++i)
          pka->pos[i] += pka->dir[i] * est_range;
        
        if (simconf->fullTraj)
          printf( "Final pos: %f %f %f\n\n", pka->pos[0], pka->pos[1], pka->pos[2]);
      }
      
      // put the recoil on the stack
      if (recoil->e > simconf->spawn_min_e && !terminate && simconf->makeRecoils)
      {
        recoil->tag = material->tag; // mark where new recoil was birthed
        
        if (recoil->type == FG)
          recoil->Ehit = recoil->e;
        
        recoil->ionId = simconf->ionId++;     // set new id then add
        recoils.push( recoil );               // add recoil to stack
        
        if (simconf->fullTraj)
        {
          double asdf[3];
          for (int i=0; i<3; ++i)
            asdf[i] = sample->w[i] / 2 - pka->pos[i];
          double RRR = sqrt( v_dot( asdf, asdf ));
          printf( "SPAWN fromcenter: %f pka->z: %d pka->e: %f recoil->z: %d recoil->e:%f\n", RRR, pka->z1, pka->e, recoil->z1, recoil->e );
        }
      }
      else delete recoil;
      
    } // endif for Rangefix == false
  } while (pka->e > pka->ef && !terminate);
  
  
  if (simconf->fullTraj && pka->tag >= 0 ) printf( "\n");
  if (simconf->fullTraj ) printf("particle killed \n\n");
}

void trimBase::doELoss( ionBase *pka, materialBase *material, double ls)
{
  pka->e -= ls * material->getrstop( pka );
  if (pka->e < 0.0)
    fprintf( stderr, " electronic energy loss stopped the ion. Broken recoil!!\n" );
}

void trimBase::rangeFix( ionBase *pka, sampleBase *sample, bool &rangefix_flag, double &ls )
{
  double a = 0;
  double b = 0;
  double c = 0;
  for (int i = 0; i < 3; i++)
  {
    if (pka->pos[i] < 0 || pka->pos[i] > simconf->length)
    {
      printf("in Trim.cc: ion wrapped around boundaries before intersection bubble");
      exit (EXIT_FAILURE);
    }

    double testPos = pka->pos[i] + pka->dir[i] * ls;
    a += sqr( testPos - pka->pos[i] );
    b += 2 * ( testPos - pka->pos[i] ) * ( pka->pos[i] - sample->w[i] / 2 );
    c += sqr( sample->w[i]/2 ) + sqr( pka->pos[i] ) - sample->w[i] * pka->pos[i];
  }
  c -= sqr( simconf->bub_rad );
  
  double d = b*b - 4 * a * c;
  
  double u[2];
  if (d > 0)
  {
    u[0] = (- b - sqrtf( d )) / (2 * a);
    u[1] = (- b + sqrtf( d )) / (2 * a);
    if (u[0] >= 0 && u[0] <= 1)
    {
      rangefix_flag = true;
      ls *=u[0];
    }
    else if (u[1] >=0 && u[1] <=1)
    {
      rangefix_flag = true;
      ls *=u[1];
    }
    if (simconf->fullTraj) printf("u1: %f u2:%f\n", u[0], u[1]);
  }
}

void trimBase::setPmax(ionBase *pka, materialBase *material)
{
  switch ( pka->pot )
  {
    case TRIM:
    {// setup avg max. impact parameter
      double epsm = pka->e * material->f; // reduced energy, epsilon
      double eeg = sqrtf( epsm * material->epsdg ); // chi (newtrim eq 7-41)
      material->pmax = material->a / ( eeg + sqrtf( eeg ) + 0.125 * pow( eeg, 0.1 ) ); // sets pmax (newtrim 7-40)
      break;
    }
    case RUTHERFORD:
    {
      double b_0 = material->a / material->f / pka->e; // reduced energy, epsilon
      double pmax_from_eng = b_0 / 2 * sqrt( material->gamma * pka->e / simconf->tmin - 1 );
      
      double sin2_min = sqr( sin( simconf->angmin * M_PI / 180) );
      double pmax_from_ang = b_0 / 2 * sqrt( 1 / sin2_min - 1);
      
      material->pmax = (pmax_from_eng > pmax_from_ang) ? pmax_from_eng : pmax_from_ang ;
      break;
    }
    case HARDSPHERE:
    {
      material->pmax = simconf->scoef[pka->z1-1].radius + material->maxr;
      break;
    }
    case NONE:
    {
      material->pmax = 100; // arbitrary
      break;
    }
  }
}

double trimBase::calcS2(ionBase *pka, materialBase *material)
{
  double r2 = dr250();
  double hh = dr250();
  
  // choose impact parameter: pmax set in trimBase::calcLs.
  double p = material->pmax * sqrtf( r2 ); // random p, weighted towards pmax
  
  // determine which atom in the material will be hit
  int nn;
  for (nn = 0; nn < material->element.size(); nn++)
  {
    hh -= material->element[nn]->t; // hh is random number
    if( hh <= 0 )
      break;
  }
  element = material->element[nn];
  
  double eps = element->fi * pka->e; // epsilon of given element
  double b = p / element->ai;        // reduced impact parameter
  
  double s2;
  switch (pka->pot)
  {
    case NONE:
      s2 = 0;
      break;
    case HARDSPHERE:
    {
      double R = simconf->scoef[pka->z1-1].radius + simconf->scoef[element->z-1].radius;
      if ( p > R)
        s2 = 0;
      else
        s2 = 1 - sqr( p / R );
      break;
    }
    case RUTHERFORD:
      s2 = 1.0 / ( 1.0 + sqr( 2.0 * eps * b ) ); // sin^2(theta_c/2) [SRIM eq. 7-15];
      break;
    case TRIM:
    {
      if (eps > 10.0)
      {
        // use rutherford scattering. includes correction (1+b*(1+b))
        s2 = 1.0 / ( 1.0 + ( 1.0 + b * ( 1.0 + b ) ) * sqr( 2.0 * eps * b ) ); // sin^2(theta_c/2) [SRIM eq. 7-15]
      }
      
      else
      {
        // first guess at ion c.p.a.
        double r = b;
        double rr = -2.7 * logf( eps * b );
        if( rr >= b )
        {
          r = rr;
          rr = -2.7 * logf( eps * rr );
          if( rr >= b ) r = rr;
        }
        
        double q, v, v1;
        
        // Constants for universal potentail, from [SRIM table 4-1]
        double u_c1 = 0.99229;
        double u_c2 = 0.011615;
        double u_c3 = 0.007122;
        double u_c5 = 9.3066;
        double u_c4 = 14.813;
        
        do
        { // Newton's method for r_0: distance of closest approach
          // universal
          // r = reduced radius x, or r/a_u
          double ex1 = 0.18175 * exp( -3.1998 * r ); // [SRIM eq. 2-74]
          double ex2 = 0.50986 * exp( -0.94229 * r );
          double ex3 = 0.28022 * exp( -0.4029 * r );
          double ex4 = 0.028171 * exp( -0.20162 * r );
          v = ( ex1 + ex2 + ex3 + ex4 ) / r; // universal screening function / reduced radius
          v1 = -( v + 3.1998 * ex1 + 0.94229 * ex2 + 0.4029 * ex3 + 0.20162 * ex4 ) / r; // -dv/dr
          
          double fr = b*b / r + v * r / eps - r; // similar to [SRIM eq. 7-2] but really denomenator in [SRIM eq. 2-79]
          double fr1 = - b*b / ( r*r ) + ( v + v1 * r ) / eps - 1.0; // d(fr)/dr
          q = fr / fr1; // f(x0)/f'(x0)
          r -= q; // Newton's method -> x1 = x0 - f(x0)/f'(x0)
        } while( fabs( q / r ) > 0.001 ); // Convergence criteria
        
        double roc = -2.0 * ( eps - v ) / v1; // R_c [SRIM eq. 7-7], calculated from 7-4
        double sqe = sqrtf( eps ); // square root of epsilon
        
        // b = p/a [SRIM eq. 7-7]
        // 5-parameter magic scattering calculation (universal pot.)
        double cc = ( u_c2 + sqe ) / ( u_c3 + sqe ); // [SRIM eq. 7-12] beta
        double aa = 2.0 * eps * ( 1.0 + ( u_c1 / sqe ) ) * pow( b, cc ); // [SRIM eq. 7-11] A
        
        double ff = ( sqrtf( aa*aa + 1.0 ) - aa ) * ( ( u_c4 + eps ) / ( u_c5 + eps ) );
        //ff = 1 / ( sqrtf( aa*aa + 1.0 ) - aa ) * ( ( u_c4 + eps ) / ( u_c5 + eps ) ); // [SRIM eq. 7-11] G
        
        double delta = ( r - b ) * aa * ff / ( ff + 1.0 );
        //delta =( r - b ) * aa / ( ff + 1.0 ); // [SRIM eq. 7-10]
        
        double co = ( b + delta + roc ) / ( r + roc ); // MAGIC formula: cosine(theta/2) [SRIM eq. 7-8]
        double c2 = co*co; //cosine(theta/2)^2
        s2 = 1.0 - c2; //sin(theta/2)^2 = 1-cos(theta/2)^2
      } // end MAGIC formulation
      //if ( pka->z1 == 92 && element->z == 92)
        //printf("%f %f\n", p, s2);
      //printf("p: %f \tpka->z: %i \ttarg->z: %i  \ts2: %f \n", p, pka->z1, element->z, s2);
      break;
    }
  }
  return s2;
}

void trimBase::calcTestS2(ionBase *pka, materialBase *material)
{
  double r2 = dr250();
  double hh = dr250();
  
  // choose impact parameter: pmax set in trimBase::calcLs.
  //double p = material->pmax * sqrtf( r2 ); // random p, weighted towards pmax
  double p = 1.25;

  // determine which atom in the material will be hit
  int nn;
//  for (nn = 0; nn < material->element.size(); nn++)
//  {
//    hh -= material->element[nn]->t; // hh is random number
//    if( hh <= 0 )
//      break;
//  }
  
  nn = 0;
  element = material->element[nn];
  cout << element->z << endl;
  
  double eps = element->fi * pka->e; // epsilon of given element
  double b = p / element->ai;        // reduced impact parameter
  
  double s2;
  cout << pka->pot << endl;
  switch (pka->pot)
  {
    case NONE:
      s2 = 0;
      break;
    case HARDSPHERE:
    {
      double R = simconf->scoef[pka->z1-1].radius + simconf->scoef[element->z-1].radius;
      cout << "p: " << p << " R: " << R << endl;
      if ( p > R)
        s2 = 0;
      else
        s2 = 1 - sqr( p / R );
      break;
    }
    case RUTHERFORD:
      s2 = 1.0 / ( 1.0 + sqr( 2.0 * eps * b ) ); // sin^2(theta_c/2) [SRIM eq. 7-15];
      break;
    case TRIM:
    {
      if (eps > 10.0)
      {
        // use rutherford scattering. includes correction (1+b*(1+b))
        s2 = 1.0 / ( 1.0 + ( 1.0 + b * ( 1.0 + b ) ) * sqr( 2.0 * eps * b ) ); // sin^2(theta_c/2) [SRIM eq. 7-15]
      }
      
      else
      {
        // first guess at ion c.p.a.
        double r = b;
        double rr = -2.7 * logf( eps * b );
        if( rr >= b )
        {
          r = rr;
          rr = -2.7 * logf( eps * rr );
          if( rr >= b ) r = rr;
        }
        
        double q, v, v1;
        
        // Constants for universal potentail, from [SRIM table 4-1]
        double u_c1 = 0.99229;
        double u_c2 = 0.011615;
        double u_c3 = 0.007122;
        double u_c5 = 9.3066;
        double u_c4 = 14.813;
        
        do
        { // Newton's method for r_0: distance of closest approach
          // universal
          // r = reduced radius x, or r/a_u
          double ex1 = 0.18175 * exp( -3.1998 * r ); // [SRIM eq. 2-74]
          double ex2 = 0.50986 * exp( -0.94229 * r );
          double ex3 = 0.28022 * exp( -0.4029 * r );
          double ex4 = 0.028171 * exp( -0.20162 * r );
          v = ( ex1 + ex2 + ex3 + ex4 ) / r; // universal screening function / reduced radius
          v1 = -( v + 3.1998 * ex1 + 0.94229 * ex2 + 0.4029 * ex3 + 0.20162 * ex4 ) / r; // -dv/dr
          
          double fr = b*b / r + v * r / eps - r; // similar to [SRIM eq. 7-2] but really denomenator in [SRIM eq. 2-79]
          double fr1 = - b*b / ( r*r ) + ( v + v1 * r ) / eps - 1.0; // d(fr)/dr
          q = fr / fr1; // f(x0)/f'(x0)
          r -= q; // Newton's method -> x1 = x0 - f(x0)/f'(x0)
        } while( fabs( q / r ) > 0.001 ); // Convergence criteria
        
        double roc = -2.0 * ( eps - v ) / v1; // R_c [SRIM eq. 7-7], calculated from 7-4
        double sqe = sqrtf( eps ); // square root of epsilon
        
        // b = p/a [SRIM eq. 7-7]
        // 5-parameter magic scattering calculation (universal pot.)
        double cc = ( u_c2 + sqe ) / ( u_c3 + sqe ); // [SRIM eq. 7-12] beta
        double aa = 2.0 * eps * ( 1.0 + ( u_c1 / sqe ) ) * pow( b, cc ); // [SRIM eq. 7-11] A
        
        double ff = ( sqrtf( aa*aa + 1.0 ) - aa ) * ( ( u_c4 + eps ) / ( u_c5 + eps ) );
        //ff = 1 / ( sqrtf( aa*aa + 1.0 ) - aa ) * ( ( u_c4 + eps ) / ( u_c5 + eps ) ); // [SRIM eq. 7-11] G
        
        double delta = ( r - b ) * aa * ff / ( ff + 1.0 );
        //delta =( r - b ) * aa / ( ff + 1.0 ); // [SRIM eq. 7-10]
        
        double co = ( b + delta + roc ) / ( r + roc ); // MAGIC formula: cosine(theta/2) [SRIM eq. 7-8]
        double c2 = co*co; //cosine(theta/2)^2
        s2 = 1.0 - c2; //sin(theta/2)^2 = 1-cos(theta/2)^2
      } // end MAGIC formulation
      //if ( pka->z1 == 92 && element->z == 92)
      //printf("%f %f\n", p, s2);
      //printf("p: %f \tpka->z: %i \ttarg->z: %i  \ts2: %f \n", p, pka->z1, element->z, s2);
      break;
    }
  }
  cout << p << ", " << pka->e << ", " << s2 << endl;
}

