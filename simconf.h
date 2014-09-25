#ifndef SIMCONF_H
#define SIMCONF_H 1

#include <string>
#include <iostream>

// ZBL coefficients a,b,c,d for all element pairs from Z=1..92
struct scoefLine {
  char sym[3], name[30];
  double mm1, m1, mnat, rho, atrho, vfermi, heat, lfctr;
  double radius;
};

enum potentials { NONE, HARDSPHERE, RUTHERFORD, TRIM };
enum boundaryCondition { PBC, INF, CUT }; // periodic, infinitly large, cut off cascades

struct simconfType {
  double ed, alfa, alpha, tmin, angmin, da, cw;
  double mdmin, mdmax; // energy ranges for ion selection
  long int ionId;
  double latp;
  
  double spawn_min_e;
  double ion_min_e;
  double fg_min_e;

  // input from commandline
  std::string run_name, fueltype;
  double bub_rad, length, fissions;
  double bub_model;
  
  // tables from files
  scoefLine scoef[92];
  double pcoef[92][8];
  double snuc[92][92][4];

  bool fullTraj;
  bool makeRecoils;
  bool BoundaryFix;
  bool AddAndKill;
  bool calc_eloss;
  bool save_eloss;
  bool monolayer;
  
  potentials pot_ff, pot_fg, pot_lat;
  boundaryCondition bounds;

  // statistics of the simulation run
  int vacancies_created;
  double KP_vacancies;
  
  void read_arg( int argc, char *argv[], bool range_only );
  
  simconfType();
private:

  void read_scoef();
  void read_snuc();
};

extern simconfType *simconf;

#endif