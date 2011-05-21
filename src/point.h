#ifndef __Point_H__
#define __Point_H__

using namespace std;

/////////////////////////////////////////////////////
//
//  Point Class
//
/////////////////////////////////////////////////////

class Vector;

class Point {
   friend class Vector;

public:
   double xy[2];

   Point() {
      xy[0] = 0;
      xy[1] = 0;
   } 
   Point(int a, int b) {
      xy[0] = a;
      xy[1] = b;
   }
   Point(double a, double b) {
      xy[0] = a;
      xy[1] = b;
   }
   ~Point() {
   };

   friend istream & operator>>(istream &, Point &);
   friend ostream & operator<<(ostream &, Point);

   Point & operator =(const Point & p);
   int operator==(Point);
   int operator!=(Point);

   Vector operator-(Point);
   Point operator+(Vector);
   Point operator-(Vector);
   Point & operator+=(Vector);
   Point & operator-=(Vector);

// Scalar Multiplication
   friend Point operator*(int, Point);
   friend Point operator*(double, Point);
   friend Point operator*(Point, int);
   friend Point operator*(Point, double);
// Scalar Division
   friend Point operator/(Point, int);
   friend Point operator/(Point, double);


// Point Addition
   friend Point operator+(Point, Point);

// Point Relations
   friend double d(Point, Point);       // Distance
   friend double d2(Point, Point);      // Distance Squared
};

/////////////////////////////////////////////////////
//
//  Vector Class
//
/////////////////////////////////////////////////////


class Vector : public Point {
public:
   Vector() : Point() {
   };
   Vector(int a, int b) : Point(a, b) {
   };
   Vector(double a, double b) : Point(a, b) {
   };
   ~Vector() {
   };

   friend ostream & operator<<(ostream &, Vector);

   Vector operator-();

// Scalar Multiplication
   friend Vector operator*(int, Vector);
   friend Vector operator*(double, Vector);
   friend Vector operator*(Vector, int);
   friend Vector operator*(Vector, double);
// Scalar Division
   friend Vector operator/(Vector, int);
   friend Vector operator/(Vector, double);


// Vector Arithmetic Operations
   Vector operator+(Vector);
   Vector operator-(Vector);
   double operator*(Vector);    // dot product
   double operator*(Point);     // dot product
   Vector operator^(Vector);    // cross product

   Vector & operator*=(double);
   Vector & operator/=(double);
   Vector & operator+=(Vector);
   Vector & operator-=(Vector);
   Vector & operator^=(Vector);

   double len() {               // vector length
      return sqrt(xy[0] *xy[0] + xy[1] *xy[1]);
   }
   double len2() {              // vector length squared
      return (xy[0] *xy[0] + xy[1] *xy[1]);
   }

   void normalize();

};

#endif
