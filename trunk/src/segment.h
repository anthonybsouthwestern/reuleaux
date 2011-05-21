#ifndef __SEGMENT__
#define __SEGMENT__

#include "curve.h" 
#include "psc.h"

// see curve.h for inherited class
class Segment : public Curve {

 public:

  Segment(PSC * c, int a, int b);
  
  Point MidPoint(int l);

  void PreSplit();
  void PreSplit(double tv);

  void SplitAtDist(int l, double r, vector<Point> & pts);
  
};

#endif
