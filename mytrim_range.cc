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
#include "sample_single.h"
#include "ion.h"
#include "trim.h"
#include "invert.h"
#include "sample.h"

#include "functions.h"


int main(int argc, char *argv[])
{
  time_t tstart, tend;
  tstart = time(0);
  
  // Settings
  bool runtrim = true;
  bool range_only = true;
  double length = 100000000.0; // size of box in A.
  
  double zin = atof( argv[3] );    // z of ion
  double min = atof( argv[4] );    // mass of ion
  double ein = atof( argv[5] );    // energy of ion in keV
  
  bool mono = true;
  if( zin == 0 || min == 0 || ein == 0 )
    mono = false;
  
  // seed randomnumber generator from system entropy pool
  FILE *urand = fopen( "/dev/random", "r" );
  int seed;
  fread( &seed, sizeof(int), 1, urand );
  fclose( urand );
  r250_init( seed<0 ? -seed : seed ); // random generator goes haywire with neg. seed
  
  // initialize global parameter structure and read data tables from file
  simconf = new simconfType;
  simconf->read_arg(argc, argv, range_only);
  
  
  // initialize sample structure. Passed values are xyz = w[3] = size of simulation
  sampleSingle *sample = new sampleSingle( length, length, length, simconf->bounds );
  sample->make_fuel( simconf->fueltype, sample, 1 );
  
  //cout << "making xe gas instead\n" << endl;
  //sample->make_fg( sample, 3, 1 );
  
  // initialize trim engine for the sample
  trimBase *trim = new trimBase( sample );
  
  // create a FIFO for recoils
  queue<ionBase*> recoils;
  
  // Declare variables
  vector<double> range; 
  vector<double> crows;
  vector<double> path;
  double norm;
  double pos1[3]; // initial position
  double dif[3];  // vector from center of bubble to location
  double crow; // range by crow's flight
  int oldionId;
  double rangeavg;
  double rangemax;
  double crowsavg;
  double crowsmax;
  double pathavg;
  double pathmax;
  
  ionBase *ff1, *pka;
  
  double A1, Etot, E1;
  int Z1;
  
  // prepare names
  char fname[200];
  snprintf( fname, 199, "output/temp/%s.%s.%.0f-%s.avgs" , argv[1], argv[2], zin, argv[5] );
  FILE *avgsFile = fopen( fname, "wt");
  
  massInverter *m = new massInverter;
  energyInverter *e = new energyInverter;
  
  printf("Fueltype: %s\n", simconf->fueltype.c_str());
  
  if (mono)
  {
    printf("Monoenergetic ion mode on\n");
    printf("Z: %0.1f \tM: %0.1f \t E[keV]: %0.1f\n", zin, min, ein);
  }
  else
    printf("FF ion mode on\n");
  
  cout << "# ion, rangeavg, crowsavg, pathavg" << endl << endl;
  // Start fissions
  for( int n = 1; n <= simconf->fissions; n++ )
  {
    oldionId = simconf->ionId;
    
    ff1 = new ionBase;
    ff1->prep_FF();
    
    if( mono )
    {
      // Assign FF data
      ff1->z1 = zin;
      ff1->m1 = min;
      ff1->e  = ein * 1.0e3; // e
    }
    else
    {
      double A1, A2, Etot, E1, E2; // inputs for ff creation
      int Z1;
      // generate fission fragment data
      A1 = m->x( dr250() ); // Randomize first mass from double hump probability
      A2 = 235.0 - A1;
      e->setMass(A1);
      Etot = e->x( dr250() ); // This E is in MeV
      E1 = Etot * A2 / ( A1 + A2 );
      Z1 = round( ( A1 * 92.0 ) / 235.0 );
      
      ff1->z1 = Z1;
      ff1->m1 = A1;
      ff1->e  = E1 * 1.0e6; // Change energy units to eV

    }
    // set direction
    ff1->dir[0] = 1.0;
    ff1->dir[1] = 0.0;
    ff1->dir[2] = 0.0;
    v_norm( ff1->dir );
    
    // set origin
    ff1->pos[0] = 0.0;
    ff1->pos[1] = 0.0;//length/2;
    ff1->pos[2] = 0.0;//length/2;
    
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
        
        pka->ef = simconf->fg_min_e;
        
        //cout << "manually setting potential!!!!" << endl;
        //pka->pot = HARDSPHERE;
        trim->trim( pka, recoils );
      
        for( int i=0; i<3; i++ ) dif[i] = pos1[i] - pka->pos[i];
        crow = sqrt( v_dot( dif, dif ) );
        //fprintf( rangeFile, "%li \t%i \t%f \t%f \t%f \t%f \n", pka->ionId, pka->z1, pka->pos[0], pka->pos[1], pka->pos[2], crow);
        
        if( pka->pos[0] >= length )
        {
          printf("box too small!");
          exit (EXIT_FAILURE);
        }
        
        range.push_back(pka->pos[0]);
        crows.push_back(crow);
        path.push_back(pka->travel);
        
        if (simconf->fullTraj)
          printf("Range: %f \tcrows: %f \tpath: %f\n",pka->pos[0], crow, pka->travel);
        
        if( (n % 10) == 0 )
        {
          rangeavg = 0;
          for(std::vector<double>::iterator j=range.begin();j!=range.end();++j) rangeavg += *j;
          rangeavg /= range.size();
          rangemax = *max_element(range.begin(), range.end());

          crowsavg = 0;
          for(std::vector<double>::iterator j=crows.begin();j!=crows.end();++j) crowsavg += *j;
          crowsavg /= crows.size();
          crowsmax = *max_element(crows.begin(), crows.end());
          
          pathavg = 0;
          for(std::vector<double>::iterator j=path.begin();j!=path.end();++j) pathavg += *j;
          pathavg /= path.size();
          pathmax = *max_element(path.begin(), path.end());
          
          cout << pka->ionId+1 << " " << rangeavg << " " << crowsavg << " " << pathavg << endl;
        }
      
        delete pka;
      } // end of recoil list
    } // end of single fission
  } // End of all fissions
  
  // output infoFile data
  
  
  rangeavg = 0;
  for(std::vector<double>::iterator j=range.begin();j!=range.end();++j) rangeavg += *j;
  rangeavg /= range.size();
  rangemax = *max_element(range.begin(), range.end());

  crowsavg = 0;
  for(std::vector<double>::iterator j=crows.begin();j!=crows.end();++j) crowsavg += *j;
  crowsavg /= crows.size();
  crowsmax = *max_element(crows.begin(), crows.end());
  
  pathavg = 0;
  for(std::vector<double>::iterator j=path.begin();j!=path.end();++j) pathavg += *j;
  pathavg /= path.size();
  pathmax = *max_element(path.begin(), path.end());
  
  tend = time(0);
  
  printf("\nfuel, \tz, \tmass, \tions, \tE [keV]\n");
  printf("%s \t\t%.0f \t%.2f \t%.0f \t%.3f \n", simconf->fueltype.c_str(), zin, min, simconf->fissions, ein);
  printf("Range [A], \tMax [A], \tCrows [A], \tMax [A], \tPath [A], \tMax [A]\n");
  printf("%.2f \t\t%.2f \t\t%.2f \t\t%.2f \t\t%.2f \t\t%.2f\n",rangeavg, rangemax, crowsavg, crowsmax, pathavg, pathmax);
  printf("\n\n%f %f\n",ein*1000.0, rangeavg);
  fprintf( avgsFile, " %.0f \t%.2f \t%.2f \t%.2f \t%.2f \t%.2f \t%.2f\n", zin, min, ein, rangeavg, rangemax, crowsavg, crowsmax);
  
  //fclose( rangeFile );
  fclose( avgsFile );
  
  cout <<"Time: " << tend - tstart << endl;
  
  //printf( "==+== %s.%s-%s Finished ==+==\n", argv[1], argv[2], argv[3] );
  
  return EXIT_SUCCESS;
}