#include "stats.h"

Stats::Stats()
{
   nloc = 0;
   nlocsteps = 0;
   worstloc = 0;
   nreject = 0;
   starttime = -1;
   endtime = -1;
   proce = 0;
   procf = 0;
   procs = 0;
   splitnotet = 0;
   splitedge = 0;
   splittri = 0;

}

void Stats::AddLocate(int numSteps)
{
   nloc++;
   nlocsteps += (double) numSteps;
   if (numSteps > worstloc) {
      worstloc = numSteps;
   }
}

void Stats::StartTime()
{
   starttime = time(NULL);
}

void Stats::EndTime()
{
   endtime = time(NULL);
}


void Stats::AddSplit(int size)
{
   if (size == 2) {
      splitedge++;
   }
   if (size == 3) {
      splittri++;
   }
}


void Stats::PrintStatsReport()
{


   cout << "Queue Items Processed: " << proce << "/" << procf << "/" <<
   procs << endl;
   cout << "nloc " << nloc << " nlocsteps " << nlocsteps << " worstloc "
        << worstloc << endl;

   cout << "Split Called on Nonexistent simplex: " << splitnotet << endl;
   cout << "Edge Splits: " << splitedge << ", Triangle Splits: " 
	<< splittri << endl;

   cout << "Number of circumcenter rejections: " << nreject << endl;

   int mytime = endtime - starttime;
   cout << "------------------------------------------------" << endl;
   cout << "Run Time: " << mytime / 60 << " min, " << mytime %
   60 << " sec" << endl;
   cout << "------------------------------------------------" << endl;

}
