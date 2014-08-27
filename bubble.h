#ifndef BUBBLE_H
#define BUBBLE_H 1

#include "interpolation.h"
#include "simconf.h"

using namespace std;

struct bubbleBase {
  double rho;  // mass density in [g/cc]
  double arho; // atomic density [atoms/Ang]?
  double _rad;  // bubble radius [Ang]
  double _temp; // temperature [K]
  
  enum eos {VDW, RONCHI};
  
  bubbleBase(double temp, eos model);
  
  eos _model;
  
  void calcArho();
};


#endif