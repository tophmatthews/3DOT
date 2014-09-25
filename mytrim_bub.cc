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
  time_t tstart, tend;
  tstart = time(0);

  /////////////////////////////////
  // ---=== Problem Setup ===--- //
  /////////////////////////////////
  
  //--- Basic Settings ---//
  bool runtrim = true;
  bool range_only = false;
  
  bool save_fsnFile   = true;
  bool save_hitFile   = false;
  bool save_escFile   = true;
  bool save_rangeFile = false;
  
  //--- initialize global parameter structure and read data tables from file ---//
  simconf = new simconfType;
  simconf->read_arg( argc, argv, range_only );
  printf( "==+== %s.%0.f-%0.f Started ==+==\n", simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
  
  //--- seed randomnumber generator from system entropy pool ---//
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed

  //--- initialize sample structure. Passed values are xyz = w[3] = size of simulation ---//
  sampleSingle *sample = new sampleSingle( simconf->length, simconf->length, simconf->length, simconf->bounds );
  sample->make_fuel( simconf->fueltype, sample, 1.0 );
  fprintf( stderr, "%s fuel built.\n", simconf->fueltype.c_str() );
  
  //--- initialize trim engine for the sample ---//
  trimBase *trim = new trimBase( sample );
  
  //--- initialze bubble structure. ---//
  bubbleBase *bubble = new bubbleBase();
  sample->make_fg( sample, bubble->rho, true );

  //--- Output file creation ---//
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
  if( save_rangeFile )
  {
    char fname[200];
    snprintf( fname, 199, "output/temp/%s.%0.f-%0.f.range" , simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
    rangeFile = fopen( fname, "wt");
  }
  
  //--- Create FIFO for recoils ---//
  queue<ionBase*> recoils;
  
  //--- Create ELoss vectors ---//
  vector<double> FF_eloss;
  vector<double> FG_eloss;
  vector<double> U_eloss;
  
  printf("Setup complete!\n");
    
  ///////////////////////////////
  // ---=== Run Problem ===--- //
  ///////////////////////////////
  printf("\n  --- Starting Fissions ---\n");
  for( int n = 1; n <= simconf->fissions; ++n )
  {
    if (simconf->fullTraj)
      printf("\t---< New Fission start >---\n");
    //--- reset counters ---//
    int escNum = 0;    // # escapees per fission
    int hitNum = 0;    // # bubble hits per fission
    int backInNum = 0; // # escapees are knocked back in per fission
    int oldionId = simconf->ionId;
    
    //--- initialize variables
    ionBase *pka;
    double dif[3] = {0};        // vector from center of bubble to location
    double fromcenter[2] = {0}; // distance of (initial/final) ion location from center of bubble
    double path[2] = {0};       // path length of (ff1/ff2)
    int pass[2] = {0};          // # times (ff1/ff2) traverses length of problem
    int punch[2] = {0};         // # times (ff1/ff2) hits bubble
    
    //--- Make fission fragments and run through collisions --///
    sample->make_FF( recoils, n );
    if (runtrim == true)
    {
      pass[0] = 0;
      while( !recoils.empty() )
      {
        pka = recoils.front();   // returns a reference to the front element
        recoils.pop();           // takes off first element
        sample->averages( pka ); // pre-calculations

        //--- pre-cascade ion analysis/processing ---//
        double starte = pka->e;
        
        // Distance from center to initial location
        if ( pka->type == FG )
          fromcenter[0] = sample->fromCenter( pka->pos );

        //--- Run pka TRIM engine ---///
        trim->trim( pka, recoils );
      
        //--- post-cascade ion analysis/processing ---//
        
        if ( simconf->save_eloss )
        {
          if ( pka->type == FG )
            FG_eloss.push_back( pka->elec_loss/starte );
          else if ( pka->type == FF )
            FF_eloss.push_back( pka->elec_loss/starte );
          else if ( pka->z1 > 90 )
            U_eloss.push_back( pka->elec_loss/starte );
        }
        
        // hitFile info
        if( save_hitFile && pka->type == FG)
        {
          // type, escaped, hits [eV]
          //fprintf( hitFile, "%i, %i, ", pka->type, pka->escapee);
          for (std::vector<long>::iterator it = pka->hit_e.begin(); it != pka->hit_e.end(); ++it)
            fprintf( hitFile, "%li, ", *it);
          fprintf( hitFile, "\n");
        }
        
        // FF statistics
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
        }// end FF statistics
        
        // FG statistics
        if (pka->type == FG)
        {
          ++hitNum;

          // Distance from center to final location
          fromcenter[1] = sample->fromCenter( pka->pos );
        
          if (pka->escapee)
          {
            ++escNum;
            if( fromcenter[1] < simconf->bub_rad ) // If end point is actually within bubble
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
  
  printf( "==+== %s.%0.f-%0.f Finished ==+==\n\n", simconf->run_name.c_str(), simconf->bub_rad, simconf->length );
  if ( simconf->save_eloss )
  {
    double FG_eloss_avg = v_avg( FG_eloss );
    double FF_eloss_avg = v_avg( FF_eloss );
    double U_eloss_avg  = v_avg( U_eloss );
    printf("Percent ELosses. FG: %.1f%%  FF: %.1f%%  U: %.1f%%\n", FG_eloss_avg*100, FF_eloss_avg*100, U_eloss_avg*100);
  }
  printf( "Simulation time [s]:\t%.0f\n", difftime(tend,tstart));
  return EXIT_SUCCESS;
}