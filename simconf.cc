#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "simconf.h"

simconfType *simconf;

simconfType::simconfType()
{
  ed = 25.0; // displacement energy
  tmin = 10.0; // Minimum transferred energy
  angmin = 10; // minimum angle
  da = 3.0; // angular grid for transmitted ions
  cw = 0.001; // channel width 1% of layer thickness
  mdmin = 10.0; // [eV] lower energy limit for ion tracking
  mdmax = 12000.0; // [eV] upper energy limit for ion tracking
  latp = 2.5; // [Ang] lattice parameter
  
  bub_rad = 0.0; // initiallization of bubble radius
  
  spawn_min_e = 300.0;
  fg_min_e = 5.0;
  ion_min_e = 300.0;
  
  fullTraj    = false; // output full trajectories
  makeRecoils = true;  // make recoils
  BoundaryFix = true;  //
  AddAndKill  = true;  // range estimation for gas atoms that reach box boundary
  calc_eloss  = true;
  monolayer   = false;
  
  // Declare potentials for atoms (NONE, HARDSPHERE, RUTHERFORD, TRIM)
  pot_ff  = TRIM;
  pot_fg  = TRIM;
  pot_lat = TRIM;
  
  // Declare boundary conditions
  bounds = PBC;
  ionId = 0; // set global ion id to zero (will be incremented for each new projectile)

  // initialize global statistics
  vacancies_created = 0;
  KP_vacancies = 0.0;
  
  // read data tables
  read_snuc();
  read_scoef();
}

void simconfType::read_arg( int argc, char *argv[], bool range_only )
{
  if ( !range_only )
  {// mytrim_bub
    if( argc != 8 ) // check if arguments are passed
    {
      fprintf( stderr, "syntax: filename bub_radius box_length fissions fueltype legacy bub_model");
      exit (EXIT_FAILURE);
    }
    run_name = argv[1];
    bub_rad = atof( argv[2] );
    length = atof( argv[3] );
    fissions = atof( argv[4] );
    fueltype = argv[5];
    
    if ( bub_rad * 2 > length )
    {
      fprintf( stderr, "simconf.cc: box size too small");
      exit (EXIT_FAILURE);
    }
    
    
    if (atof(argv[6]) == 1)
    {
      printf( "Legacy calculation on\n");
      pot_ff = RUTHERFORD;
      pot_fg = HARDSPHERE;
      pot_lat = HARDSPHERE;
      calc_eloss = false;
    }
    simconf->bub_model = atof( argv[7] );
  }
  else
  {//mytrim_range
    if ( argc != 8 )
    {
      fprintf( stderr, "syntax: filename, fueltype, Z, M, E[keV], #");
      exit (EXIT_FAILURE);
    }
    fueltype = argv[2];
    fissions = atof( argv[6] );

    if (atof(argv[7]) == 1)
    {
      printf( "Legacy calculation on\n");
      pot_ff = RUTHERFORD;
      pot_fg = HARDSPHERE;
      pot_lat = HARDSPHERE;
      calc_eloss = false;
    }
    
    // Other range settings
    makeRecoils = false;
    BoundaryFix = false;
    AddAndKill  = false;
  }
}

void simconfType::read_snuc()
{
  char buf[200], fname[500];

  //snprintf( fname, 500, "%s/SNUC03.dat", getenv("HOME") );
  snprintf( fname, 500, "%s/SNUC03.dat", DATA_DIR );
  FILE *sf = fopen( fname, "rt" );
  if( sf == 0 )
  {
    fprintf( stderr, "unable to open %s\n", fname );
    exit(1);
  }
  for( int i = 0; i < 92; i++ )
    for( int j = i; j < 92; j++ )
    {
      fscanf( sf, "%*d %*d %lf %lf %lf %lf\n",
        &snuc[j][i][0], &snuc[j][i][1], &snuc[j][i][2], &snuc[j][i][3] );
      for( int n = 0; n < 4; n++ )
        snuc[i][j][n] = snuc[j][i][n];
    }
  fclose( sf );
}

void simconfType::read_scoef()
{
  char buf[2001], fname[500];
  FILE *sf;

  //snprintf( fname, 500, "%s/SCOEF.95A", getenv("HOME") );
  snprintf( fname, 500, "%s/SCOEF.95A", DATA_DIR );
  sf = fopen( fname, "rt" );
  fgets( buf, 2000, sf ); // header
  fgets( buf, 2000, sf ); // header
  for ( int i = 0; i < 92; ++i )
  {
    fscanf( sf, "%*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &scoef[i].mm1, &scoef[i].m1, &scoef[i].mnat,
      &scoef[i].rho, &scoef[i].atrho, &scoef[i].vfermi, &scoef[i].heat,
      &pcoef[i][0], &pcoef[i][1], &pcoef[i][2], &pcoef[i][3],
      &pcoef[i][4], &pcoef[i][5], &pcoef[i][6], &pcoef[i][7] );
    //printf( "%f %f\n", scoef[i].mm1, pcoef[i][7] );
  }
  fclose( sf );

  //snprintf( fname, 500, "%s/SLFCTR.dat", getenv("HOME") );
  snprintf( fname, 500, "%s/SLFCTR.dat", DATA_DIR );
  sf = fopen( fname, "rt" );
  fgets( buf, 2000, sf ); // header
  for ( int i = 0; i < 92; ++i )
    fscanf( sf, "%*d %lf\n", &scoef[i].lfctr );
  fclose( sf );

  snprintf( fname, 500, "%s/ELRAD.dat", DATA_DIR );
  sf = fopen( fname, "rt" );
  int z;
  for ( int i = 0; i < 92; ++i )
  {
    fscanf( sf, "%*d %s %s %lf\n", scoef[i].sym, scoef[i].name, &scoef[i].radius );
    scoef[i].radius /= 100;
  }
  fclose( sf );
}

