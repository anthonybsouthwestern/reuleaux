#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

// 
// Triangle data structure
// 
// Each triangle stores a list of its vertex IDs
// as well as a list of IDs of the neighboring
// triangles
//
// Also, the radius and circumcenter of each triangle
// are stored since this data is used frequently
//

class Triangle {
 public:
  // list of vertices in the simplex
  int v[3];
  // list of neighbors
  Triangle *n[3];
  // IMPORTANT: n[0] is through the 12 edge
  // IMPORTANT: n[1] is through the 02 edge
  // IMPORTANT: n[2] is through the 01 edge
  
  double rsq;
  Point cc;

  int visited;

  Triangle();
  Triangle(Triangle * s);
  
  Triangle(int a, int b, int c);
  
  
  void set(int a, int b, int c);
  void link(Triangle * s1, Triangle * s2, Triangle * s3);
  int isNeighbor(Triangle * s);
  
  void SetVertex(int i, int myv);
  void SetNeighbor(int i, Triangle * myn);
  
  bool contains(int i);

  bool ncontains(Triangle * i);
  
  // useful utilities...
  friend ostream & operator<<(ostream &, Triangle *);
  
 protected:
  
 private:
};

bool operator<(Triangle & x, Triangle & y);

bool operator<(const Face & x, const Face & y);

#endif
