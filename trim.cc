#include <math.h>
#include <stdio.h>
#include <iostream>

#include "trim.h"

// does a single ion cascade
void trimBase::trim( ionBase *pka_, queue<ionBase*> &recoils)
{
  pka = pka_;
  
  // make recoil queue available in overloadable functions
  recoil_queue_ptr = &recoils;

  pka->ic = 0;       // reset hit counter
  terminate = false; // reset kill flag
  
  //--- Initialize variables ---//
  double s2, c2, ct, st; // trig identities
  double den;            // energy transfered to recoil atom
  double rdir[3], perp[3], norm, psi;
  double p1, p2;         // momentum before and after collision
  double r1 = dr250();   // random number
  double oldls = 0.0;
  double newls = 0.0;
  
  if (simconf->fullTraj)
  {
    printf( "\n+++=== ION START pos: %f %f %f ===+++\n", pka->pos[0], pka->pos[1], pka->pos[2]);
    printf( "\tion z: %i  e: %f \n", pka->z1, pka->e);
    printf( "\tdir: %f %f %f\n", pka->dir[0], pka->dir[1], pka->dir[2]);
    printf( "\tpot: %i (0 is NONE, 1 is HARDSPHERE, 2 is RUTHERFORD, 3 is TRIM\n", pka->pot);
  }
  
  do // cycle for each collision
  {
    ++pka->ic; // number of collisions for each ion in each material
    material = sample->lookupMaterial( pka->pos );
    oldMaterial = material;
    if (simconf->fullTraj)
      printf( "\n\tnow in material with avg z: %0.1f", material->az);

    if (material == 0) break;

    v_norm( pka->dir ); // normalize direction vector
    
    //--- Calculate path length ls ---//
    setPmax( pka, material );
    double ls = 1.0 / ( M_PI * sqr( material->pmax ) * material->arho ); // (newtrim eq 7-28); // calculate path length
    
    if (simconf->monolayer)
      ls = material->minls;
    
    if (simconf->fullTraj)
      printf( "\n\tls: %f\tpos: %f %f %f\tdir: %f %f %f\n", ls,pka->pos[0], pka->pos[1], pka->pos[2], pka->dir[0], pka->dir[1], pka->dir[2]);
    
    if ( newls!=0 )
    {
      ls = newls;
      if (simconf->fullTraj) printf( "\t Continuined from previous oldls: %f newls: %f\n", oldls, newls);
    }
    else if ( pka->ic == 1 ) // if first recoil of atom
    {
      ls = ls * dr250();
      if (simconf->fullTraj) printf( "\tJust entered material. Adjusted ls: %f\n", ls);
    }
    
    if (ls < 0)
    {
      printf("Trim.cc: negative travel distance\n");
      exit (EXIT_FAILURE);
    }
    
    ls = (ls > material->minls ? ls : material->minls);
    
//    if ( pka->ic == 1)
//      ls -= 1.0; //subtract extra that was used to push it across before
    
    //printf("ls: %f material: %f\n", ls, material->az);
    
    newls = 0;
    if ( simconf->bub_rad != 0 )
    {
      //--- Check if pka crosses bubble surface ---//
      if ( bubbleCrossFix(pka, sample, ls) )
        continue; // if crosses bubble boundary, start over
    
      oldls = ls;
      //--- Check if pka crosses boundary surface ---//
      if ( boundsCrossFix(pka, sample, ls) )
      {
        newls = oldls - ls;
        continue; // if crosses sample boundary, start over
      }
    }
    
//--- TESTING BLOCK ---//
//    cout << "TESTING STILL" << endl;
//    double thisp = 1.1;
//    int thismat = 0;
//    pka->e = 1000;
//    s2 = calcTestS2(pka, material, thisp, thismat); //nn = 0 for u, nn = 1 for c
//    den = element->ec * s2 * pka->e;
//    printf("%f\n", s2);
//    break;
//
//    pka->e = 1000;
//    s2 = calcTestS2(pka, material, thisp, thismat); //nn = 0 for u, nn = 1 for c
//    den = element->ec * s2 * pka->e;
//    printf("%f\n", den);
//    
//    pka->e = 30000;
//    s2 = calcTestS2(pka, material, thisp, thismat); //nn = 0 for u, nn = 1 for c
//    den = element->ec * s2 * pka->e;
//    printf("%f\n", den);
//
//    printf("pka->e: %f\n",pka->e);
//    
//    for(double i=0.1;i<=2.5;i+=0.05)
//    {
//      s2 = calcTestS2(pka, material, i, 0); //nn = 0 for u, nn = 1 for c
//      printf("%f %f\n", i, s2);
//    }
//    
//    break;
//--- TESTING BLOCK ---//
    
    
//--- If doesn't hit bubble or boundary, then it hits an atom ---///
    s2 = calcS2( pka, material ); // sin^2(theta_c/2)
    c2 = 1.0 - s2;                // cos^2(theta_c/2)
    ct = 2.0 * c2 - 1.0;          // cos(theta_c)
    st = sqrtf( 1.0 - ct*ct );    // sin(theta_c)
    
    // energy transferred to recoil atom
    den = element->ec * s2 * pka->e; // T=gamma*E*sin^2(theta_c/2). gamma*E = Tmax

    // advance clock pathlength/velocity
    pka->t += 10.1811859 * ls / sqrt( 2.0 * pka->e / pka->m1 );
    // time in fs! m in u, l in Ang, e in eV
    // 1000g/kg, 6.022e23/mol, 1.602e-19J/eV, 1e5m/s=1Ang/fs 1.0/0.09822038
    //printf( "se %d  %f [eV]  %f [keV/nm]  %f [nm]\n", pka->ionId, pka->e, see/100.0, pl/10.0 );

    if ( simconf->calc_eloss )
    {
      doELoss( pka, material, ls );
      if ( pka->e < 0 )
        cout << "broken from hit" << endl;
    }
    
    pka->hit_e.push_back( floor(pka->e) );
    
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
    
    recoil = pka->spawnRecoil();
    recoil->e = den;              // assign recoil energy
    recoil->e -= element->Elbind; // subtract lattice binding energy from recoil
    recoil->m1 = element->m;      // assign mass
    recoil->z1 = element->z;      // assign atomic number
    recoil->assignType();
    
    //--- Create a random vector perpendicular to pka->dir ---//
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

    //--- Calculate new direction, subtract from old dir (stored in recoil) ---//
    for (int i = 0; i < 3; ++i)
    {
      recoil->dir[i] = pka->dir[i] * p1; // initial pka momentum vector
      pka->dir[i] += perp[i] * sin( psi );
      recoil->dir[i] -= pka->dir[i] * p2;
    }
    v_norm( recoil->dir );
    
    //--- Put the recoil on the stack ---//
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
          asdf[i] = sample->w[i] / 2.0 - pka->pos[i];
        double RRR = sqrt( v_dot( asdf, asdf ));
        printf( "SPAWN fromcenter: %f pka->z: %d pka->e: %f recoil->z: %d recoil->e:%f\n", RRR, pka->z1, pka->e, recoil->z1, recoil->e );
      }
    }
    else
      delete recoil;
    
  } while (pka->e > pka->ef && !terminate);
  
  
  if (simconf->fullTraj && pka->tag >= 0 ) printf( "\n");
  if (simconf->fullTraj ) printf("particle killed \n\n");
}

// ==================================================================================== //

void trimBase::doELoss( ionBase *pka, materialBase *material, double ls)
{
  double olde = pka->e;
  double rstop = material->getrstop( pka );
  pka->e -= ls * material->getrstop( pka );
  if (pka->e < 0.0)
  {
    fprintf( stderr, "\n electronic energy loss stopped the ion. Broken recoil!!\n" );
    fprintf( stderr, "pos: %f %f %f  olde: %f newe: %f\n", pka->pos[0], pka->pos[1], pka->pos[2], olde, pka->e);
    fprintf( stderr, "ls %f getrstop %f \n", ls, rstop );
    fprintf( stderr, "az: %f\n", material->az);
    exit (EXIT_FAILURE);
  }
}

// ==================================================================================== //

void trimBase::moveIonToBoundary( ionBase *pka, double ls )
{
  pka->ic = 0; // reset collision count
  
  if (simconf->calc_eloss) doELoss( pka, material, ls );
  
  pka->travel += ls; // add travel length to total path length
  
  for (int i = 0; i < 3; ++i)
    pka->pos[i] += ls * pka->dir[i]; // update position
}

// ==================================================================================== //

bool trimBase::bubbleCrossFix( ionBase *pka, sampleBase *sample, double &ls )
{
  // ---=== Check if Crossed first, and fix ls if true ===--- //
  bool crossed = false;
  double in_or_out = 0;;
  double a = 0;
  double b = 0;
  double c = 0;
  for (int i = 0; i < 3; i++)
  {
    double testPos = pka->pos[i] + pka->dir[i] * ls;
    a += sqr( testPos - pka->pos[i] );
    b += 2.0 * ( testPos - pka->pos[i] ) * ( pka->pos[i] - sample->w[i] / 2.0 );
    c += sqr( sample->w[i]/2.0 ) + sqr( pka->pos[i] ) - sample->w[i] * pka->pos[i];
  }
  c -= sqr( simconf->bub_rad );
  
  double d = b*b - 4.0 * a * c;
  
  if (d > 0)
  {
    double u[2];
    u[0] = (- b - sqrtf( d )) / (2.0 * a);
    u[1] = (- b + sqrtf( d )) / (2.0 * a);
    if (u[0] >= 0 && u[0] <= 1)
    {
      crossed = true;
      ls *= u[0];
      in_or_out = -1;
    }
    else if (u[1] >=0 && u[1] <=1)
    {
      crossed = true;
      ls *= u[1];
      in_or_out = 1;
    }
  }
  
  // ---=== If crossed then apply fixed ls =w==--- //
  if ( crossed )
  {
    pka->punch++; // add to punch count
    
    double oldmat = material->az;
    
    moveIonToBoundary ( pka, ls );
    
    if ( pka->e < 0 )
      cout << "broken from cross" << endl;
    
    double bub_norm[3];
    for (int i = 0; i < 3; ++i)
      bub_norm[i] = pka->pos[i] - simconf->length/2.0; // update position
    v_norm(bub_norm);
    

    for (int i = 0; i < 3; ++i)
      pka->pos[i] += in_or_out * 1.0 * bub_norm[i]; // update position
    
    materialBase *newmat;
    newmat = sample->lookupMaterial( pka->pos );
    if (newmat->az == material->az || newmat->am == material->am)
    {
      fprintf( stderr, "\n Crossed boundary but didn't register new material\n" );
      fprintf( stderr, "pos: %f %f %f e: %f\n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e);
      fprintf( stderr, "dir: %f %f %f\n", pka->dir[0], pka->dir[1], pka->dir[2]);
      fprintf( stderr, "direction: %f\n", in_or_out);
      exit (EXIT_FAILURE);
    }


    
    if ( pka->type == FG )
    {
      if ( !pka->escapee )
      {
        pka->escapee = true;
        pka->Eout = pka->e; // set transition energy
      }
    }
    
    if ( simconf->fullTraj )
    {
      printf( "CROSS at pos: %f %f %f e: %f ls*u: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, ls );
      printf( "\tLeaving material with average Z %0.1f\n", oldmat);
    }
  }
  return crossed;
}

// ==================================================================================== //

bool trimBase::boundsCrossFix( ionBase *pka, sampleBase *sample, double &ls )
{
  bool boundsCrossFix_test = false; // flag for testing output for this method
  
  if ( simconf->bounds != PBC )
  {
    printf("Non-periodic boundary conditions are not currently supported");
    exit (EXIT_FAILURE);
  }
  
  // ---=== First test if crossed the sample boundary ===--- //
  bool crossed = false;
  int crossed_dir[3] = {0};
  
  double testPos[3]; // test position
  for ( int i=0; i<3; ++i )
  {
    testPos[i] = pka->pos[i] + ls * pka->dir[i];
    
    if ( testPos[i] < 0 )
    {
      crossed = true;
      crossed_dir[i] = -1;
    }
    else if ( testPos[i] > sample->w[i] )
    {
      crossed = true;
      crossed_dir[i] = +1;
    }
    else
      crossed_dir[i] = 0;
  }
  if ( boundsCrossFix_test )
    printf ( "testPos: %f %f %f\n", testPos[0], testPos[1], testPos[2]);
  
  // ---=== Move particle to edge if crosses ===--- //
  
  if ( !crossed )
    return crossed;
  
  else
  { // --- Calculate distance to planes --- //
    if ( boundsCrossFix_test )
      cout << "hey crossed! " << crossed_dir[0] << " " << crossed_dir[1] << " " << crossed_dir[2] << endl;
    double d[3] = { ls }; //initially set all distances to path length. can only get smaller from here
    int whichone = -1;
    for ( int i=0; i<3; ++i )
    {
      if ( crossed_dir[i] != 0 )
      {
        double p_0 = ( sample->w[i] + crossed_dir[i] * sample->w[i] )/2.0;
        d[i] = ( p_0 - pka->pos[i] ) / pka->dir[i];
        if ( boundsCrossFix_test ) cout << "i: " << i << " d[i]: " << d[i] << endl;
        ls = ( ls < d[i] ? ls : d[i] );
        whichone = ( ls < d[i] ? whichone : i );
      }
    }
    if (boundsCrossFix_test)
      printf("the fix is for the %i direction. pos: %f %f %f", whichone, pka->pos[0], pka->pos[1], pka->pos[2]);
    if ( whichone < 0)
    {
      printf("Trim.cc: something is wrong with int whichone");
      exit (EXIT_FAILURE);
    }
    
    // --- Move ion to plane of first intersection --- //
    moveIonToBoundary ( pka, ls );
    
    if ( pka->e < 0 )
      cout << "broken from wrapped" << endl;
    
    //---< White boundary condition test >---\\
    
    //Choose coord 0,1,2, is +x, +y, +z
    int coord = floor(dr250() * 3);
    
    //choose which plane, 0 is +, 1 is -
    //double xxx = floor(dr250() * 2);
    double plane_sign = floor(dr250() * 2);
    
    //printf("coord: %i  plane_sign %f\n", coord, plane_sign);
    
    // Set position
    for (int i = 0; i < 3; ++i)
    {
      if (i != coord )
      {
        pka->pos[i] = simconf->length * dr250();
        pka->dir[i] = ( dr250() * 2.0 ) - 1.0;
      }
      else
      {
        pka->pos[i] = simconf->length - plane_sign * simconf->length; // if 1, then = 0, if 0 then on far surface
        pka->dir[i] = dr250() * floor( plane_sign * 2 - 0.5);
      }
    }
    
    v_norm(pka->dir);
    //printf("pos: %f %f %f\n",pka->pos[0], pka->pos[1], pka->pos[2]);
    //printf("dir: %f %f %f\n\n",pka->dir[0], pka->dir[1], pka->dir[2]);
  
    
    //pka->pos[whichone] += simconf->bit * pka->dir[whichone]; // update position
    
    ++pka->pass;
    
    if ( simconf->fullTraj )
      printf( "WRAPPED at pos: %f %f %f e: %f ls: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, ls );
    
    if ( pka->type != FG || !simconf->AddAndKill )
    {// --- Wrap particle and check fix only occurs once --- //
      for ( int i=0; i<3; ++i)
      { // if greater than or less than
        if ( pka->pos[i] < 0 || pka->pos[i] > sample->w[i] )
        {
          pka->pos[i] -= pka->pos[i] / abs(pka->pos[i]) * sample->w[i];
        }
      }
      
      if ( simconf->fullTraj )
        printf( "\tNow at pos: %f %f %f\n", pka->pos[0], pka->pos[1], pka->pos[2] );
      
      return crossed;
    }
    else
    {
      if (simconf->fullTraj)
        printf( "ADDANDKILLED at pos: %f %f %f e: %f ls: %f \n", pka->pos[0], pka->pos[1], pka->pos[2], pka->e, ls );
      
      terminate = true;
      double est_range = pka->RangeInFuel( simconf->fueltype );
      if (simconf->fullTraj)
        printf( "\tADDANDKILLED e: %f  range: %f dir: %f %f %f\n", pka->e, est_range, pka->dir[0], pka->dir[1], pka->dir[2]);
      
      // add travel length to total pka path length
      // will not work well for total travel length, since this is actually range. But good enough!
      pka->travel += est_range;
      for (int i = 0; i < 3; ++i)
        pka->pos[i] += pka->dir[i] * est_range;
      
      if (simconf->fullTraj)
        printf( "\n+++=== ION END pos: %f %f %f\n\n", pka->pos[0], pka->pos[1], pka->pos[2]);
      
      return crossed;
    }
  }
}

// ==================================================================================== //


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
      double pmax_from_eng = b_0 / 2.0 * sqrt( material->gamma * pka->e / simconf->tmin - 1.0 );
      
      double sin2_min = sqr( sin( simconf->angmin * M_PI / 180.0) );
      double pmax_from_ang = b_0 / 2.0 * sqrt( 1.0 / sin2_min - 1.0);
      
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
      material->pmax = 100.0; // arbitrary
      break;
    }
  }
}

// ==================================================================================== //

double trimBase::calcS2(ionBase *pka, materialBase *material)
{
  double r2 = dr250();
  double hh = dr250();
  
  // choose impact parameter: pmax set in trimBase::calcLs.
  if (simconf->monolayer)
    material->pmax = 2.17;
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
        s2 = 0.0;
      else
        s2 = 1.0 - sqr( p / R );
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

// ==================================================================================== //

double trimBase::calcTestS2(ionBase *pka, materialBase *material, double p, int nn)
{
  double r2 = dr250();
  double hh = dr250();
  
  // choose impact parameter: pmax set in trimBase::calcLs.
  //double p = material->pmax * sqrtf( r2 ); // random p, weighted towards pmax

  // determine which atom in the material will be hit
  //int nn;
//  for (nn = 0; nn < material->element.size(); nn++)
//  {
//    hh -= material->element[nn]->t; // hh is random number
//    if( hh <= 0 )
//      break;
//  }
  
  //nn = ;
  element = material->element[nn];
  //cout << element->z << endl;
  
  double eps = element->fi * pka->e; // epsilon of given element
  double b = p / element->ai;        // reduced impact parameter
  
  double s2;
  //cout << pka->pot << endl;
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
  //cout << p << ", " << pka->e << ", " << s2 << endl;
  
  
  
  return s2;
}

// ==================================================================================== //



