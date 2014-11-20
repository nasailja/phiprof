/*
This file is part of the phiprof library

Copyright 2011, 2012 Finnish Meteorological Institute

Phiprof is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library. If not, see <http://www.gnu.org/licenses/>.


Modified in 2014-11 from https://github.com/fmihpc/phiprof see the
version at https://github.com/nasailja/phiprof for more info
*/


#include "mpi.h"
#include <iostream>
#include <iomanip> 
#include <string>    
#include <vector>     
#include "phiprof.hpp"

using namespace std;

double compute(double seconds){
   const double overhead=0.000002; //estimated
   double a,b,c,d;
   int  i,iterations=1000;
   double t1,t2;

   a=1.00004;
   b=1e-10;
   c=1.00011;
   d=0;

   t2=t1=MPI_Wtime();
   while(t2-t1<seconds-overhead){
      for(i=0;i<iterations;i++){
         d+=a*b+c;
      }
      t2=MPI_Wtime();
   }
   return d;
}


int main(int argc,char **argv){

   int rank;
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   const int nIterations=1000000;
   if(rank==0)
      cout << "Measuring performance of start-stop calls" <<endl;

   phiprof::Phiprof profiler;
   profiler.start("Benchmarking phiprofr"); 
   profiler.start("Initalized timers using ID");

   if(rank==0)
      cout << "  1/3" <<endl;
   int id_a=profiler.initializeTimer("a","A with ID");
   for(int i=0;i<nIterations;i++){
      profiler.start(id_a);
      profiler.stop(id_a);
   }
   profiler.stop("Initalized timers using ID",nIterations,"start-stop");

   if(rank==0)
      cout << "  2/3" <<endl;
   profiler.start("Re-initialized timers using ID");
   for(int i=0;i<nIterations;i++){
      id_a=profiler.initializeTimer("a","A with ID");
      profiler.start(id_a);
      profiler.stop(id_a);
   }
   profiler.stop("Re-initialized timers using ID",nIterations,"start-stop");

   if(rank==0)
      cout << "  3/3" <<endl;
   profiler.start("Timers using labels");
   for(int i=0;i<nIterations;i++){
      profiler.start("a");
      profiler.stop("a");
   }
   profiler.stop("Timers using labels",nIterations*2,"start-stop");
   profiler.stop("Benchmarking phiprofr"); 

   MPI_Barrier(MPI_COMM_WORLD);

   if(rank==0)
      cout << "Measuring accuracy" <<endl;
   profiler.start("Test accuracy");
   if(rank==0)
      cout << "  1/2" <<endl;
   profiler.start("100x0.1s computations"); 
   for(int i=0;i<100;i++){
      profiler.start("compute");
      compute(0.1);
      profiler.stop("compute");
   }
   profiler.stop("100x0.1s computations");

   if(rank==0)
      cout << "  2/2" <<endl;
   MPI_Barrier(MPI_COMM_WORLD);
   profiler.start("100x0.1s computations + logprofile"); 
   for(int i=0;i<100;i++){
      profiler.start("compute");
      compute(0.1);
      profiler.printLogProfile(MPI_COMM_WORLD,i);
      profiler.printLogProfile(MPI_COMM_WORLD,i,"profile_log_maxlev1"," ",1);
      profiler.stop("compute");
   }
   profiler.stop("100x0.1s computations + logprofile"); 

   profiler.stop("Test accuracy");
   
   MPI_Barrier(MPI_COMM_WORLD);
   double t1=MPI_Wtime();
   profiler.print(MPI_COMM_WORLD);
   profiler.print(MPI_COMM_WORLD,"profile_minfrac0.01",0.01);
   if(rank==0) cout<< "Print time is "<<MPI_Wtime()-t1<<endl;

   MPI_Finalize();
}
