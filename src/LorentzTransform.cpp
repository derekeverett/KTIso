#pragma once
#include <math.h>
#include "Parameter.h"

class fourVector {
private:

public:
  fourVector();
  ~fourVector();

  float p0, p1, p2, p3;

  fourVector boost_by_u(fourVector p, fourVector u);


  };

  fourVector::fourVector() {

  }

  fourVector::~fourVector() {
  }

  //use this function to initialize energy density within JETSCAPE
  fourVector fourVector::boost_by_u(fourVector p, fourVector u)
  {

  }
