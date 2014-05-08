#ifndef TRIM_H
#define TRIM_H 1

#include <vector>
#include <queue>

#include "material.h"
#include "sample.h"

using namespace std;

class trimBase {
public:
  void trim( ionBase *pka, queue<ionBase*> &recoils );
  trimBase( sampleBase *sample_ ) : sample( sample_ ) {};

protected:
  sampleBase *sample;
  ionBase *pka, *recoil;
  materialBase *material;
  materialBase *testMaterial;
  elementBase *element;
  queue<ionBase*> *recoil_queue_ptr;
  bool terminate;
  bool edged;

  // by default only follow recoils with E > 12eV
  virtual bool spawnRecoilLimit() { return recoil->e > 12.0; };
  virtual void vacancyCreation() {};
};

#endif
