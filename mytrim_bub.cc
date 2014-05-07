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

#include "simconf.h"
#include "sample_clusters.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include "sample.h"
#include "functions.h"

int main(int argc, char *argv[])
{
  time_t tstart, tend;
  tstart = time(0);
  
  string buildtime(__TIME__); // time of build
  string builddate(__DATE__);
  
  // Settings
  bool runtrim = true;
  
  sampleBase::sampleBoundary bounds = sampleBase::PBC;
  
  bool save_infoFile  = false;
  bool save_fsnFile   = true;
  bool save_hitFile   = false;
  bool save_escFile   = true;
  bool save_rangeFile = false;
  
  if( argc != 6 ) // check if arguments are passed
  {
    fprintf( stderr, "syntax: filename bub_radius box_length fissions fueltype");
    return 1;
  }
  
  printf( "==+== %s.%s-%s Started ==+==\n", argv[1], argv[2], argv[3] );
  
  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed

  // initialize global parameter structure and read data tables from file
  simconf = new simconfType;
  simconf->bub_rad = atof( argv[2] );;
  simconf->fueltype = argv[5];
  
  // Convert inputs to floats
  double length = atof( argv[3] ); // size of box in A
  double fissions = atof( argv[4] ); // number of fissions to run

  // initialize sample structure. Passed values are xyz = w[3] = size of simulation
  sampleClusters *sample = new sampleClusters( length, length, length, bounds );
  
  // initialize trim engine for the sample
  trimBase *trim = new trimBase( sample );

  sample->initSpatialhash( int( sample->w[0] / simconf->bub_rad ) - 1,
                           int( sample->w[1] / simconf->bub_rad ) - 1,
                           int( sample->w[2] / simconf->bub_rad ) - 1 );

  sample->addCluster( length/2, length/2, length/2, simconf->bub_rad); // Add bubble in center of box
  sample->make_fuel( simconf->fueltype, sample, 0.9 );
  double bub_den = calc_rho(simconf->bub_rad, 132.0);
  sample->make_fg( simconf->fueltype, sample, bub_den );

  fprintf( stderr, "%s sample built.\n", simconf->fueltype.c_str() );
  std::cout << "Bubble density [g/cc]: " << bub_den << endl;
  
  // create a FIFO for recoils
  queue<ionBase*> recoils;

  double norm;                   // used to determine direction
  double pos1[3];                // initial position
  double dif[3];                 // vector from center of bubble to location
  double fromcenter[2];          // distance of (initial/final) ion location from center of bubble
  double path[2];                // path length of (ff1/ff2)
  int oldionId;                  // save for old ion id
  int pass[2];                   // # times (ff1/ff2) traverses length of problem
  int punch[2];                  // # times (ff1/ff2) hits bubble
  int hitNum, escNum, backInNum; // # times fg atoms is hit/escapes/knocked back in

  massInverter *m = new massInverter;
  energyInverter *e = new energyInverter;
  ionBase *ff1, *ff2, *pka;

  double A1, A2, Etot, E1, E2; // inputs for ff creation
  int Z1, Z2;

  // File creation
  
  FILE * hitFile = NULL;
  if( save_hitFile )
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%s-%s.hit" , argv[1], argv[2], argv[3] );
    hitFile = fopen( fname, "wt");
  }

  FILE * escFile = NULL;
  if( save_escFile )
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%s-%s.esc" , argv[1], argv[2], argv[3] );
    escFile = fopen( fname, "wt");
  }
  
  FILE * fsnFile = NULL;
  if( save_fsnFile )
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%s-%s.fsn" , argv[1], argv[2], argv[3] );
    fsnFile = fopen( fname, "wt");
  }
  
  FILE * rangeFile = NULL;
  if( save_rangeFile)
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%s-%s.range" , argv[1], argv[2], argv[3] );
    rangeFile = fopen( fname, "wt");
  }

  // Start fissions
  for( int n = 1; n <= fissions; ++n )
  {
    // reset counters
    escNum = 0;
    hitNum = 0;
    backInNum = 0;
    oldionId = simconf->ionId;
    
    // generate fission fragment data
    A1 = m->x( dr250() ); // Randomize first mass from double hump probability
    A2 = 235.0 - A1;
    e->setMass(A1);
    Etot = e->x( dr250() ); // This E is in MeV
    E1 = Etot * A2 / ( A1 + A2 );
    E2 = Etot - E1;
    Z1 = round( ( A1 * 92.0 ) / 235.0 );
    Z2 = 92 - Z1;
    
    // -- Spawn 1st fission fragment -- //
    ff1 = new ionBase;
    ff1->prep_FF();
    ff1->z1 = Z1;
    ff1->m1 = A1;
    ff1->e  = E1 * 1.0e6; // Change energy units to eV
    ff1->ionId = simconf->ionId++;
    
    // Random direction
    do
    { 
      for( int i = 0; i < 3; ++i ) ff1->dir[i] = dr250() - 0.5;
      norm = v_dot( ff1->dir, ff1->dir );
    } while( norm <= 0.0001 );
    v_scale( ff1->dir, 1.0 / sqrtf( norm ) );

    // random origin
    for( int i = 0; i < 3; ++i )
    {
      ff1->pos[i] = dr250() * sample->w[i];
      ff1->pos0[i] = ff1->pos[i]; // set orinal position
    }
    recoils.push( ff1 );

    // -- Spawn 2nd fission fragments -- //
    ff2 = new ionBase( *ff1 ); // copy constructor
    for( int i = 0; i < 3; ++i )
      ff2->dir[i] *= -1.0; // reverse direction
    ff2->z1 = Z2;
    ff2->m1 = A2;
    ff2->e  = E2 * 1.0e6;
    ff2->ionId = simconf->ionId++;
    recoils.push( ff2 );

    printf( "%s Fsn %i/%.0f Z1=%d (%.2f MeV)\t Z2=%d (%.2f MeV)\n", argv[1], n, fissions, Z1, E1, Z2, E2 );
  
    if (runtrim == true)
    {
      pass[0] = 0;
      while( !recoils.empty() )
      {
        pka = recoils.front();   // returns a reference to the front element
        recoils.pop();           // takes off first element
        sample->averages( pka ); // pre-calculations

        // -- pre-cascade ion analysis/processing -- //
        
        if( pka->tag >= 0) // if pka is xenon and not FF
        {
          // mark the first recoil that falls into the MD energy gap with 1 (child generations increase the number)
          if( pka->e > simconf->mdmin && pka->e < simconf->mdmax && pka->md == 0 )
            pka->md = 1;

          for( int i = 0; i < 3; i++ )
          {
            pos1[i] = pka->pos[i];
            dif[i] = sample->c[i][pka->tag] - pos1[i];
          }
          fromcenter[0] = sqrt( v_dot( dif, dif ) );
        }

        trim->trim( pka, recoils );
      
        // -- post-cascade ion analysis/processing -- //
        
        if( pka->gen == 0 ) // if ff
        {
          if( pass[0] == 0 ) // if first ff
          {
            pass[0]  = pka->pass;
            punch[0] = pka->punch;
            path[0]  = pka->travel;
          }
          else
          {
            pass[1]  = pka->pass;
            punch[1] = pka->punch;
            path[1]  = pka->travel;
          }
        }
        
        if( pka->tag >= 0 )
        {
          ++hitNum;

          for( int i = 0; i < 3; i++ )
            dif[i] = sample->c[i][pka->tag] - pka->pos[i];
          fromcenter[1] = sqrt( v_dot( dif, dif ) );
          
          // hitFile info
          if( save_hitFile )
          {
            // id, gen, pass, punch, r02cen [pm], rf2cen [pm], travel [pm],
            fprintf( hitFile, "%li\t%i\t%i\t%i\t", pka->ionId, pka->gen, pka->pass, pka->punch);
            fprintf( hitFile, "%f\t%f\t%f\t", fromcenter[0], fromcenter[1], pka->travel+0.5);
            fprintf( hitFile, "%f\t%f\t", pka->Ehit, pka->Eout);
            fprintf( hitFile, "%i\t%i\t%i", pka->fam_fuel, pka->fam_fg, pka->fam_parent);
            fprintf( hitFile, "\n");
          }
        
          if( pka->escapee )
          {
            ++escNum;
            if( fromcenter[1] < simconf->bub_rad )
              ++backInNum;
            if( save_escFile )
            {
              fprintf( escFile, "%li\t%i\t%i\t%i\t", pka->ionId, pka->gen, pka->pass, pka->punch);
              fprintf( escFile, "%f\t%f\t%f\t", fromcenter[0], fromcenter[1], pka->travel);
              fprintf( escFile, "%f\t%f\t", pka->Ehit, pka->Eout);
              fprintf( escFile, "%i\t%i\t%i", pka->fam_fuel, pka->fam_fg, pka->fam_parent);
              fprintf( escFile, "\n");
            }
          }// done with escaped fg data
        }// done with fg data
        delete pka;
      }// end of recoil list
    } // end of single fission
    
    if( save_fsnFile )
    {
      fprintf(fsnFile, "%i\t", n);
      fprintf(fsnFile, "%i\t%.3f\t%i\t%i\t%.1f\t", Z1, E1, pass[0], punch[0], path[0]);
      fprintf(fsnFile, "%i\t%.3f\t%i\t%i\t%.1f\t", Z2, E2, pass[1], punch[1], path[1]);
      fprintf(fsnFile, "%li\t%i\t%i\t%i\n", simconf->ionId - oldionId, hitNum, escNum, backInNum);
    }
  } // End of all fissions

  tend = time(0); // output infoFile data
  
  if( save_hitFile ) fclose( hitFile );
  if( save_escFile ) fclose( escFile );
  if( save_fsnFile ) fclose( fsnFile );
  if( save_rangeFile ) fclose( rangeFile );
  
  printf( "==+== %s.%s-%s Finished ==+==\n", argv[1], argv[2], argv[3] );
  printf( "Simulation time [s]:\t%.0f\n", difftime(tend,tstart));
  return EXIT_SUCCESS;
}