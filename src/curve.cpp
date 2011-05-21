#include "allh.h"

// note: circles should replace the constructor
//         to get the circular connections working correctly...
Curve::Curve(PSC * c) {
  this->c = c;
  type = CURVE_GENERIC;
}
/*// this code moved to the segment class...
Curve::Curve(PSC * c, int a, int b) {
  this->c = c;
  left[b] = a;
  right[b] = -1;
  left[a] = -1;
  right[a] = b;

  leftEndpoint = a;
  rightEndpoint = b;
}
*/
/*
void Segment::InsertVertex(int i) {
  InsertVertex(i,left.begin()->first);
}
*/

// insert a vertex if the subsegment it is 
// going to break is known...
void Curve::InsertVertex(int i, int l) {
  int r = right[l];

  if (l == -1 || r == -1) {
    cout << "ERROR: Segment::InsertVertex() failed" << endl;
    cout << "l="<< l <<", r=" << r << endl;
    exit(0);
  }

  right[l] = i;
  right[i] = r;
  left[i] = l;
  left[r] = i;

  c->pts->l2s[i].insert(this);

}

bool operator<(const Subsegment & a, const Subsegment & b)
{


  if (a.leftEndpoint < b.leftEndpoint) {
    return true;
  } else if (a.leftEndpoint == b.leftEndpoint && 
	     a.inputCurve < b.inputCurve) {
    return true;
  }
  return false;
}



// subsegment iterator

Curve::iterator::iterator() {

}

Curve::iterator::iterator(int i, Curve * s) {
  current = i;
  curve = s;
}

void Curve::iterator::operator++() {
  if (current == -1) {
    cout << "ERROR: Problem incrementing Segment::iterator" 
	 << endl;
  }
  current = curve->right[current];
}

void Curve::iterator::operator++(int a) {
  ++(*this);
}

Curve::iterator & Curve::iterator::operator=(const Curve::iterator & it) {
  current = it.current;
  curve = it.curve;
  return *this;
}


bool Curve::iterator::operator==(const Curve::iterator & it) {
  if (curve == it.curve && current == it.current) {
    return true;
  }
  return false;
}

bool Curve::iterator::operator!=(const Curve::iterator & it) {
    return !(*this == it);
}

int Curve::iterator::operator*() {
  return current;
}

Curve::iterator Curve::begin() {
  return iterator(leftEndpoint,this);
}

Curve::iterator Curve::end() {
  return iterator(rightEndpoint,this);
}
