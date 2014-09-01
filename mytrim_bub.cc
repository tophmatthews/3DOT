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
#include "sample_single.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include "sample.h"
#include "functions.h"
#include "bubble.h"

int main(int argc, char *argv[])
{
  if( argc != 8 ) // check if arguments are passed
  {
    fprintf( stderr, "syntax: filename bub_radius box_length fissions fueltype legacy");
    return 1;
  }
  
  time_t tstart, tend;
  tstart = time(0);
  
  string buildtime(__TIME__); // time of build
  string builddate(__DATE__);
  
  // Settings
  bool runtrim = true;
  
  sampleBase::sampleBoundary bounds = sampleBase::PBC;
  
  bool save_fsnFile   = true;
  bool save_hitFile   = false;
  bool save_escFile   = true;
  bool save_rangeFile = false;
  

  //TODO: pce these all to simconf through initializaiton
  // initialize global parameter structure and read data tables from file
  simconf = new simconfType;
  simconf->run_name = argv[1];
  simconf->bub_rad = atof( argv[2] );
  simconf->length = atof( argv[3] );
  simconf->fissions = atof( argv[4] );
  simconf->fueltype = argv[5];
  
  if (atof(argv[6]) == 1)
  {
    printf( "Legacy calculation on\n");
    simconf->pot_ff = RUTHERFORD;
    simconf->pot_fg = HARDSPHERE;
    simconf->pot_lat = HARDSPHERE;
    simconf->calc_eloss = false;
  }

  simconf->bub_model = atof( argv[7] );

  printf( "==+== %s.%0.f-%0.f Started ==+==\n", simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
  
  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed

  // initialize sample structure. Passed values are xyz = w[3] = size of simulation
  sampleSingle *sample = new sampleSingle( simconf->length, simconf->length, simconf->length, bounds );

  sample->make_fuel( simconf->fueltype, sample, 0.9 );
  
  fprintf( stderr, "%s sample built.\n", simconf->fueltype.c_str() );
  
  // initialize trim engine for the sample
  trimBase *trim = new trimBase( sample );
  // initialze bubble structure.
  bubbleBase *bubble = new bubbleBase();
  
  //double bub_den = calc_rho(simconf->bub_rad, 132.0);
  sample->make_fg( sample, bubble->rho, false );

  std::cout << "Bubble density [g/cc]: " << bubble->rho << endl;
  
  // create a FIFO for recoils
  queue<ionBase*> recoils;

  //double norm;                   // used to determine direction
  double pos1[3];                // initial position
  double dif[3];                 // vector from center of bubble to location
  double fromcenter[2];          // distance of (initial/final) ion location from center of bubble
  double path[2];                // path length of (ff1/ff2)
  int oldionId;                  // save for old ion id
  int pass[2];                   // # times (ff1/ff2) traverses length of problem
  int punch[2];                  // # times (ff1/ff2) hits bubble
  int hitNum, escNum, backInNum; // # times fg atoms is hit/escapes/knocked back in

  ionBase *pka;


  // File creation
  
  FILE * hitFile = NULL;
  if( save_hitFile )
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%0.f-%0.f.hit" , simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
    hitFile = fopen( fname, "wt");
  }

  FILE * escFile = NULL;
  if( save_escFile )
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%0.f-%0.f.esc" , simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
    escFile = fopen( fname, "wt");
  }
  
  FILE * fsnFile = NULL;
  if( save_fsnFile )
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%0.f-%0.f.fsn" , simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
    fsnFile = fopen( fname, "wt");
  }
  
  FILE * rangeFile = NULL;
  if( save_rangeFile)
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%0.f-%0.f.range" , simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
    rangeFile = fopen( fname, "wt");
  }

  // Start fissions
  for( int n = 1; n <= simconf->fissions; ++n )
  {
    // reset counters
    escNum = 0;
    hitNum = 0;
    backInNum = 0;
    oldionId = simconf->ionId;
    
    ionBase *ionTrash;
    ionTrash->make_FF( recoils, n );
  
    if (runtrim == true)
    {
      pass[0] = 0;
      while( !recoils.empty() )
      {
        pka = recoils.front();   // returns a reference to the front element
        recoils.pop();           // takes off first element
        sample->averages( pka ); // pre-calculations

        // -- pre-cascade ion analysis/processing -- //
        
        if ( pka->type == FG )
        {
          for ( int i = 0; i < 3; ++i )
          {
            pos1[i] = pka->pos[i];
            dif[i] = sample->w[i] - pos1[i];
          }
          fromcenter[0] = sqrt( v_dot( dif, dif ) );
        }

        trim->trim( pka, recoils );
      
        // -- post-cascade ion analysis/processing -- //
        
        // hitFile info
        if( save_hitFile && pka->type == FG)
        {
          // type, escaped, hits [eV]
          //fprintf( hitFile, "%i, %i, ", pka->type, pka->escapee);
          for (std::vector<long>::iterator it = pka->hit_e.begin(); it != pka->hit_e.end(); ++it)
            fprintf( hitFile, "%li, ", *it);
          fprintf( hitFile, "\n");
        }
        
        if ( pka->type == FF ) // if ff
        {
          if ( pass[0] == 0 ) // if first ff
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
        
        if (pka->type == FG)
        {
          ++hitNum;

          for ( int i = 0; i < 3; ++i )
            dif[i] = sample->w[i] - pka->pos[i];
          fromcenter[1] = sqrt( v_dot( dif, dif ) );
        
          if (pka->escapee)
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
    
    if ( save_fsnFile )
    {
      fprintf(fsnFile, "%i\t", n);
      fprintf(fsnFile, "%i\t%.3f\t%i\t%i\t%.1f\t", 1, 1.0, pass[0], punch[0], path[0]);
      fprintf(fsnFile, "%i\t%.3f\t%i\t%i\t%.1f\t", 1, 1.0, pass[1], punch[1], path[1]);
      fprintf(fsnFile, "%li\t%i\t%i\t%i\n", simconf->ionId - oldionId, hitNum, escNum, backInNum);
    }
  } // End of all fissions

  tend = time(0); // output infoFile data
  
  if ( save_hitFile ) fclose( hitFile );
  if ( save_escFile ) fclose( escFile );
  if ( save_fsnFile ) fclose( fsnFile );
  if ( save_rangeFile ) fclose( rangeFile );
  
  printf( "==+== %s.%0.f-%0.f Finished ==+==\n", simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
  printf( "Simulation time [s]:\t%.0f\n", difftime(tend,tstart));
  return EXIT_SUCCESS;
}