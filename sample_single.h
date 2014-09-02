#ifndef SAMPLE_SINGLE_H
#define SAMPLE_SINGLE_H 1

#include "sample.h"
#include "functions.h"

struct sampleSingle : sampleBase {

  sampleSingle( double x = 10000.0, double y = 10000.0, double z = 10000.0, boundaryCondition b = PBC );

  virtual materialBase* lookupMaterial( double* pos );
};

#endif