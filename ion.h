#ifndef ION_H
#define ION_H 1

#include <vector>
#include <string>
#include <queue>

#include "invert.h"
#include "r250.h"
#include "functions.h"
#include "simconf.h"
//#include "sample_single.h"

using namespace std;

enum ion_type {FF, LAT, FG};

struct ionBase {
  int z1;
  int ic; // Hit counter for trim engine
  double m1, e;
  double dir[3], pos[3], pos0[3], posOld[3]; // normalized velocity vector, and position, and birthed position
  double t; // internal clock
  double travel; // pathlength
  int pass;  // # times through box
  int punch; // # times through bubble
  bool escapee;
  double Ehit, Eout;
  int fam_fuel, fam_fg, fam_parent;
  potentials pot;
  
  std::vector<int> famtree; // family tree of ion
  
  std::vector<long> hit_e; // energies of each hit

  int tag, gen;
  long ionId;
  ion_type type;

  double ef;

  ionBase();
  ionBase( ionBase* prototype );
  ionBase( int _z1, double _m1, double _e );
  virtual ~ionBase() {};

  virtual void parent( ionBase* parent );
  virtual ionBase* spawnRecoil();
  
  virtual void assignType();
  
  void reset();
  double RangeInFuel( std::string fueltype );
  virtual void prep_FF();
  virtual void calcLs();
};

#endif
