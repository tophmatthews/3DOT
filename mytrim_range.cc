/***************************************************************************
 *   Copyright (C) 2008 by Daniel Schwen   *
 *   daniel@schwen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <ctime>
#include <iostream>
#include <string>
#include <algorithm>

#include "simconf.h"
#include "element.h"
#include "material.h"
#include "sample_clusters.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include "sample.h"

#include "functions.h"


int main(int argc, char *argv[])
// argc = num of arguments at command line, including exe name,
// argv = pointer to supplied arguments
{
  time_t tstart, tend;
  tstart = time(0);
  
  string buildtime(__TIME__); // time of build
  string builddate(__DATE__);
  
  // Settings
  bool calc_bub_rho = false;
  bool runtrim = true;
  sampleBase::sampleBoundary bounds = sampleBase::CUT;
  
  double rangeavg;
  double rangemax;
  
  char fname[200];
  if( argc !=  7 ) // check if arguments are passed
  {
    //fprintf( stderr, "aasyntax:\n%s basename r Cbfactor\n\nCbfactor=1 => 7e-4 bubbles/nm^3\n", argv[0] );
    fprintf( stderr, "syntax: filename, length, fueltype, Z, M, E[keV], #");
    return 1;
  }
  
  //printf( "==+== %s.%s Started ==+==\n", argv[1], argv[2] );
  
  // Convert inputs to floats   
  
  double length = 1000000; // size of box in A.
  double zin = atof( argv[3] );    // z of ion
  double min = atof( argv[4] );    // mass of ion
  double ein = atof( argv[5] );    // energy of ion in keV
  double fissions = atof( argv[6] ); // number of ions to run
  
  
  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed
  
  // initialize global parameter structure and read data tables from file
  simconf = new simconfType;
  simconf->fueltype = argv[2];
  simconf->makeRecoils = false;
  simconf->BoundaryFix = false;
  simconf->AddAndKill  = false;
  
  // initialize sample structure. Passed values are xyz = w[3] = size of simulation
  sampleClusters *sample = new sampleClusters( length, length, length, bounds );
  
  // initialize trim engine for the sample
  trimBase *trim = new trimBase( sample );
  
  //fprintf( stderr, "sample built.\n" );
  
  materialBase *material;
  elementBase *element;
  
  if( simconf->fueltype == "uc" ) // uranium carbide
  {
    material = new materialBase( 12.3 ); // rho
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
    element = new elementBase;
    element->z = 6;
    element->m = 12.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else if( simconf->fueltype == "un" ) // uranium nitride
  {
    material = new materialBase( 12.9 ); // rho
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
    element = new elementBase;
    element->z = 7;
    element->m = 14.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else if( simconf->fueltype == "um" ) // uranium metal
  {
    material = new materialBase( 17.1 ); // rho
    element = new elementBase;
    element->z = 92;
    element->m = 238.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else if( simconf->fueltype == "xe" ) // uranium metal
  {
    material = new materialBase( 2.0 ); // rho
    element = new elementBase;
    element->z = 54;
    element->m = 132.0;
    element->t = 1.0;
    material->element.push_back( element );
  }
  else
  {
    fprintf( stderr, "Invalid fuel type specified");
  }
  
  material->prepare(); // all materials added
  sample->material.push_back( material ); // add material to sample

  
  // create a FIFO for recoils
  queue<ionBase*> recoils;
  
  vector<double> range;
  vector<double> crows;
  double norm;
  double pos1[3]; // initial position
  double dif[3];  // vector from center of bubble to location
  double crow; // range by crow's flight
  int oldionId;
  
  ionBase *ff1, *pka;
  
  double A1, Etot, E1;
  int Z1;
  
  snprintf( fname, 199, "output/temp/%s.%s.%.0f-%s.range" , argv[1], argv[2], zin, argv[5] );
  FILE *rangeFile = fopen( fname, "wt");
  
  snprintf( fname, 199, "output/temp/%s.%s.%.0f-%s.avgs" , argv[1], argv[2], zin, argv[5] );
  FILE *avgsFile = fopen( fname, "wt");
  
  // Start fissions
  for( int n = 1; n <= fissions; n++ )
  {
    oldionId = simconf->ionId;
    
    ff1 = new ionBase;
    ff1->gen = 0;  // generation (0 = FF)
    ff1->tag = -1; // -1 = born in fuel
    ff1->md = 0;

    // Assign FF data
    ff1->z1 = zin;
    ff1->m1 = min;
    ff1->e  = ein * 1e3; // eV

    // set direction
    ff1->dir[0] = 1;
    ff1->dir[1] = 0;
    ff1->dir[2] = 0;
    v_norm( ff1->dir );

    // set origin
    ff1->pos[0] = 0;
    ff1->pos[1] = length/2;
    ff1->pos[2] = length/2;

    ff1->ionId = simconf->ionId++;
    recoils.push( ff1 );
    
    if (runtrim == true)
    {
      while( !recoils.empty() )
      {
        pka = recoils.front();   // returns a reference to the front element
        recoils.pop();           // takes off first element
        sample->averages( pka ); // pre-calculations
        for( int i=0; i<3; i++ ) pos1[i] = pka->pos[i];
        
        trim->trim( pka, recoils );
      
        for( int i=0; i<3; i++ ) dif[i] = pos1[i] - pka->pos[i];
        crow = sqrt( v_dot( dif, dif ) );
        //fprintf( rangeFile, "%li \t%i \t%f \t%f \t%f \t%f \n", pka->ionId, pka->z1, pka->pos[0], pka->pos[1], pka->pos[2], crow);
        
        range.push_back(pka->pos[0]);
        crows.push_back(crow);
        
        if( (n % 1000) == 0 )
        {
          rangeavg = 0;
          for(std::vector<double>::iterator j=range.begin();j!=range.end();++j) rangeavg += *j;
          rangeavg /= range.size();
          rangemax = *max_element(range.begin(), range.end());
          cout << pka->ionId+1 << " " << rangeavg << " " << rangemax << endl;
        }
      
        delete pka;
      } // end of recoil list
    } // end of single fission
  } // End of all fissions
  
  // output infoFile data
  tend = time(0);
  
  
  
  double crowsavg;
  for(std::vector<double>::iterator j=crows.begin();j!=crows.end();++j) crowsavg += *j;
  crowsavg /= crows.size();
  double crowsmax = *max_element(crows.begin(), crows.end());
  //cout << cavg << endl;
  printf("fueltype, z, m, ions, E [keV], Range [A], Max [A], Crows [A], Max [A]\n");
  printf("%s %.0f \t%.2f \t%.0f \t%.3f \t%.2f \t%.2f \t%.2f \t%.2f\n", argv[2], zin, min, fissions, ein, rangeavg, rangemax, crowsavg, crowsmax);
  fprintf( avgsFile, " %.0f \t%.2f \t%.2f \t%.2f \t%.2f \t%.2f \t%.2f\n", zin, min, ein, rangeavg, rangemax, crowsavg, crowsmax);
  
  fclose( rangeFile );
  fclose( avgsFile );
  
  //printf( "==+== %s.%s-%s Finished ==+==\n", argv[1], argv[2], argv[3] );
  
  return EXIT_SUCCESS;
}