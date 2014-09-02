#ifndef SAMPLE_CLUSTERS_H
#define SAMPLE_CLUSTERS_H 1

#include "sample.h"

struct sampleClusters : sampleBase {

  double sd, kd[3]; // half the spatial diagonal of a hash block, hash block size
  int *sh, kn[3]; // spatial hash and its dimensions

  int *cl, cn, cnm; // cluster linklist, actual number of clusters (incl. ghosts) and number reserved
  double *c[4]; // three arrays for cluster x,y,z,r^2 coordinates
  double cmr; // maximum cluster radius in the sample

  sampleClusters( double x = 10000.0, double y = 10000.0, double z = 10000.0, boundaryCondition b = PBC );

  virtual materialBase* lookupMaterial( double* pos );

  int lookupCluster( double* pos, double dr = 0.0 );
  void initSpatialhash( int x, int y, int z );
  void clearSpatialHash();
  void clearClusters();
  void addCluster( double x, double y, double z, double r );
  void addRandomClusters( int n, double r, double dr = 0.0 );
protected:
  void reallocClusters( int n );
};

#endif