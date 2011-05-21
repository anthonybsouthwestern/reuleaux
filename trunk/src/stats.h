#ifndef __STATS_H__
#define __STATS_H__

#include <iostream>
#include <time.h>

//
// Stats data structure
//
// This data structure is used for keeping track of 
// information about the mesh generation procedure
//



using namespace std;

class Stats {
public:


  // some of these are currently unused...
  double nloc;
  double nlocsteps;
  double worstloc;
  void AddLocate(int numSteps);
  
  double nreject;
  
  int starttime;
  int endtime;
  int logtime;
  void StartTime();
  void EndTime();
  
  int proce;
  int procf;
  int procs;
  
  int splitnotet;
  int splitedge;
  int splittri;
  void AddSplit(int size);
  
  
  Stats();
  
  void PrintStatsReport();
  
 private:
};

#endif                          // !defined(__STATS_H__)
