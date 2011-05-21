#ifndef __CURVE__
#define __CURVE__

#include "psc.h"

enum CurveType
  {CURVE_GENERIC,CURVE_SEGMENT, CURVE_CIRCLE, CURVE_ARC, CURVE_BEZIER, CURVE_CUBBEZ, CURVE_NBEZIER};

class Curve {

 public:

  CurveType type;

  string name;
  int id;

  PSC * c;

  Curve(PSC * c);
  //Curve(PSC * c, int a, int b);
  
  // insert vertex i in subsegment (l,right[l])
  void InsertVertex(int i, int l); 

  // virtual functions must be implemented for each type of curve
  // considered...
  
  // return the midpoint of the arc (l,right[l])
  virtual Point MidPoint(int l) {Point p; return p;};


  // add vertices to a curve so that no subcurve has total variation 
  // in orientation of longer than totalVariation
  // if no argument is given, is a dfault value of pi/4
  virtual void PreSplit() {};
  virtual void PreSplit(double totalVariation) {};
  // example of presplit: given an input circle containing no
  //  input verties, PreSplit() should add 4 vertices to the circle.


  // find points that are distance r away from (intput) vertex l
  // this is an important operation for adding vertices when handling
  // acute input angles
  //
  // this should not insert the vertex into the mesh or curve, only
  // return the appropriate point(s) to be handled by the main program
  virtual void SplitAtDist(int l, double r, vector <Point> & pts) {};

  void SetName(string s) {name = s;}
  
  // curve iterator
  class iterator {
  public:
    iterator();
    iterator(int, Curve *);
    void operator++();
    void operator++(int);
    iterator& operator=(const iterator & it);
    bool operator==(const iterator & it);
    bool operator!=(const iterator & it);
    int operator*();
     

  private:
    int current;
    Curve * curve;

  };
  
  iterator begin();
  iterator end();
  

  int leftEndpoint,rightEndpoint;

  map<int,int> left;
  map<int,int> right;
  
 private:


  
};

class Subsegment {
 public:
  int leftEndpoint;
  Curve * inputCurve;
 private:
};

bool operator<(const Subsegment & x, const Subsegment & y);

#endif
