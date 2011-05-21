Reuleaux Code Base
5-21-2011

1.  Building the Code

In many cases simply type the following command in the root directory:

make

If there are build errors, it is likely related to some configuration
options that can be set in file: src/makefile

A successful code build should create the executable bin/reuleaux 
  (or bin/reuleaux.exe in windows/cygwin).


2.  Running the Code

There are two ways to run the code.  In both cases, commands are given
  to run the code from the top level Reuleaux directory.

A) Directly call the executable followed by the filename of the 
    data file:

./bin/reuleaux data/ReuleauxExample/ReuleauxExample

This loads both data/ReuleauxExample/ReuleauxExample.psc for input complex 
and data/su/ReuleauxExample.cfg for the configuration options.

B) Use the run target in the makefile:

make run f=ReuleauxExample

This performs the exact same call as in method A.  


3.  Source Code Overview

Here is a brief description of the main source files.

reuleaux.cpp: The main program is located here.
psc.cpp: The main data structure for the "mesh", i.e. this keeps
         pointers to all the important data structures and has 
         routines to perform main funcitons (Ruppert's algorithm)
triangulation.cpp: Data structure to hold a traingulation.  The 
         triangulation is stored as a list of triangles.
triangle.cpp: The triangle data structure.  A triangle contains
         a list of IDs of its vertices and a list of pointers to
         the neighboring triangles.
trianglepool.cpp: Does memory management for triangulation.cpp.
points.cpp: Data structure to keep a list of all the vertices
         in the mesh.  The IDs in each triangle point to this
         structure.
point.cpp: Simple classes for storing points and vectors.
geometry.cpp: Routines for performing geometric computations, i.e.
         computing circumcenters and incircle tests.
curve.cpp: The generic curve class.  
segment.cpp: Inherits from the curve class for a special type of
	 curves: straight line segments.
stats.cpp: Place to keep track of stats during the mesh generation.
         Currently stats doesn't do very much.
configuration.cpp: Loads and stores configuration info from .cfg.
iniparser.cpp/dictionary.cpp: An open source code for parsing input
         files.  Used by configuration.cpp.  More info at:
         http://ndevilla.free.fr/iniparser
predicates.c: Code for performing exact insphere checks.  From:
         http://www.cs.cmu.edu/~quake/robust.html  