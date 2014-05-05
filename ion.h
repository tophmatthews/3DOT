#ifndef ION_H
#define ION_H 1

#include <vector>

struct ionBase {
  int z1;
  double m1, e;
  double dir[3], pos[3], pos0[3]; // normalized velocity vector, and position, and birthed position
  double t; // internal clock
  double travel; // pathlength
  int pass;  // # times through box
  int punch; // # times through bubble
  bool escapee;
  double Ehit, Eout;
  int fam_fuel, fam_fg, fam_parent;
  
  std::vector<int> famtree; // family tree of ion

  int tag, gen;
  long ionId;
  int md; // generation after first ion falling into the MD energy gap ( 200eV - 12000eV )

  double ef;

  ionBase();
  ionBase( ionBase* prototype );
  ionBase( int _z1, double _m1, double _e );
  virtual ~ionBase() {};

  virtual void parent( ionBase* parent );
  virtual ionBase* spawnRecoil();

  void set_ef();
  void reset();
  void prep_FF();
};

#endif
