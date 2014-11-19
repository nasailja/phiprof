/*
This file is part of the phiprof library

Copyright 2011, 2012 Finnish Meteorological Institute

Phiprof is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Phiprof is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#include "iostream"

#include "mpi.h"
#include "phiprof.hpp"

int main(int argc,char **argv){

   int rank;
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   phiprof::Phiprof profiler;
   profiler.start("Greeting");
   std::cout << "Hello world from rank " << rank << std::endl;
   profiler.stop("Greeting", 1, "greetings");

   MPI_Barrier(MPI_COMM_WORLD);
   profiler.print(MPI_COMM_WORLD);
   if (rank == 0) {
	   std::cout << "Profiling results were printed into profile_*.txt" << std::endl;
   }

   MPI_Finalize();
   return 0;
}
