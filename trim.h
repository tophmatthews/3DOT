#ifndef TRIM_H
#define TRIM_H 1

#include <vector>
#include <queue>

#include "ion.h"
#include "material.h"
#include "sample.h"
#include "simconf.h"
#include "functions.h"

using namespace std;

class trimBase {
public:
  void trim( ionBase *pka, queue<ionBase*> &recoils );
  trimBase( sampleBase *sample_ ) : sample( sample_ ) {};

protected:
  sampleBase *sample;
  ionBase *pka, *recoil;
  materialBase *material, *oldMaterial;
  elementBase *element;
  queue<ionBase*> *recoil_queue_ptr;
  bool terminate;
  bool edged;

  virtual void doELoss( ionBase *pka, materialBase *material, double ls);
  virtual void rangeFix( ionBase *pka, sampleBase *sample, bool &rangefix_flag, double &ls );

  virtual void setPmax(ionBase *pka, materialBase *material);
  virtual double calcS2(ionBase *pka, materialBase *material);
  virtual void calcTestS2(ionBase *pka, materialBase *material);
};

#endif
