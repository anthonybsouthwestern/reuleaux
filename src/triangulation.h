#ifndef __TRIANGULATION__
#define __TRIANGULATION__

//
// triangulation.h
// 
// class for managing the Delaunay triangulation
//
// Long-term TODO List:
//  * Implement RemovePoint routines
//  * Improve output functions to handle smooth features
//
class Triangulation {
public:
  
  PSC * c;
  
  double bbox;         // size of the bounding box...
  double upperbound[2];
  double lowerbound[2];
  
  map < int, Triangle * >p2t;   // map from point to some simplex containing it...
  
  int np;          // number of points
  int nt;          // number of triangles

  Point center;
  
  // constructor...
  Triangulation(PSC * c) {this->c = c; np = 0; nt = 0;}
  
  // high level operations...
  void Neighbors(int i, set < int >&nbs);
  void Neighbors(int i, vector < int >&nbs);
  void Triangles(int p, set < Triangle * >&star);
  void Triangles(int p, vector < Triangle * >&star);
  
  // find the triangle that contains p in its circumcircle.
  Triangle *Locate(Point p);
  Triangle *Locate(int i);
  // hint is a starting point...
  Triangle *Locate(Point p, int hint);
  Triangle *Locate(int i, int hint);
  // s is a starting simplex...
  Triangle *Locate(Point p, Triangle * s);
  Triangle *Locate(int i, Triangle * s);
  
  // clear the cavity and add new triangles
// s is the starting simplex for the search for i
// cavity returns all new simplices added...
   void InsertPoint(int i);
   void InsertPoint(int i, int hint);
   void InsertPoint(int i, Triangle * s);
   void InsertPoint(int i, set < Triangle * >&oldcavity,
                    set < Triangle * >&newcavity);
   void InsertPoint(int i, int hint, set < Triangle * >&oldcavity,
                    set < Triangle * >&newcavity);
   void InsertPoint(int i, Triangle * s, set < Triangle * >&oldcavity,
                    set < Triangle * >&newcavity);
   void InsertPoint(int i, vector < Triangle * >&oldcavity,
                    vector < Triangle * >&newcavity);
   void InsertPoint(int i, int hint, vector < Triangle * >&oldcavity,
                    vector < Triangle * >&newcavity);
   void InsertPoint(int i, Triangle * s, vector < Triangle * >&oldcavity,
                    vector < Triangle * >&newcavity);
   
   // Utility Function...
   // Insert Point into triangulation but do not update p2s...
   // make this private? (2-18-08)

   // no longer used... remove
   //void InsertPointNoP2T(int i, Triangle * s, set < Triangle * >&oldcavity,
   //                     set < Triangle * >&newcavity);


   // FIXME: this has not been implemented...
   void RemovePoint(int i, set < Triangle * >&cavity);
   void RemovePoint(int i, set < Triangle * >&oldcavity,
                    set < Triangle * >&cavity);
   
   // return the cavity and the boundary...
   void Cavity(int p, Triangle * t, set < Triangle * >&cavity,
               set < Face > &bndry);
   void Cavity(Point & p, Triangle * t, set < Triangle * >&cavity,
               set < Face > &bndry);
   //void Cavity(Point & p, set < int >&nbs);
   void Cavity(Point & p, Triangle * t, set < int >&nbs);

   void Cavity(int p, Triangle * t, vector < Triangle * >&cavity,
               vector < Face > &bndry);
   void Cavity(Point & p, Triangle * t, vector < Triangle * >&cavity,
               vector < Face > &bndry);
   // come back to these...
   //void Cavity(Point & p, vector < int >&nbs);
   void Cavity(Point & p, Triangle * t, vector < int >&nbs);

   
   // add the new tets with the new point
   void Stitch(int i, set < Face > &bndry, set < Triangle * >&newcavity);
   void Stitch(int i, vector < Face > &bndry, vector < Triangle * >&newcavity);
   
   void SetRadiusAndCC(Triangle * s);
      
   // set up the bounding box...
   void Initialize();
   
   void SetP2t(Triangle * t);
   
   
   // triangle iterator
   class titerator {
   public:
     titerator();
     titerator(Triangle * s);
     void operator++();
     void operator++(int);
     titerator& operator=(const titerator & it);
     bool operator==(const titerator & it);
     bool operator!=(const titerator & it);
     Triangle * operator*();
     
   private:
     set < Triangle * >v;
     deque < Triangle * >q;
   };
   
   titerator tbegin();
   titerator tend();
   
   // vertex iterator
   class viterator {
   public:
     viterator();
     viterator(map< int, Triangle *>::iterator m);
     void operator++();
     void operator++(int);
     viterator& operator=(const viterator & it);
     bool operator==(const viterator & it);
     bool operator!=(const viterator & it);
     int operator*();
     
   private:
     map < int, Triangle * >::iterator mit;
     
   };
   
   viterator vbegin();
   viterator vend();
   
 private:
   
   
};

#endif
