/*
Phiprof - Parallel Hierarchical Profiler

Copyright 2010, 2011, 2012 Finnish Meteorological Institute
Copyright 2014 Ilja Honkonen

Phiprof is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Phiprof is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Modified in 2014-11-19 from the version at
https://github.com/fmihpc/phiprof see
https://github.com/nasailja/phiprof for details
*/

#ifndef PHIPROF_HPP
#define PHIPROF_HPP


#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>


#ifdef _OPENMP
#include "omp.h"
#endif

#include "mpi.h"

#ifdef CRAYPAT
#include "pat_api.h"
#endif


namespace phiprof
{

class Phiprof {
public:

     struct TimerData{
	 std::string label;          //print label 
	 std::string workUnitLabel;   //unit for the counter workUnitCount
	 double workUnits;        // how many units of work have we done. If -1 it is not counted, or printed
	 std::vector<std::string> groups; // What user-defined groups does this timer belong to, e.g., "MPI", "IO", etc..
	
	 int id; // unique id identifying this timer (index for timers)
	 int parentId;  //key of parent (id)
	 int threads; //threads active when this timer was called
  	 std::vector<int> childIds; //children of this timer
	
	 int level;  //what hierarchy level
	 int count; //how many times have this been accumulated
	 double time; // total time accumulated
	 double startTime; //Starting time of previous start() call

         bool active;
      };
      //used with MPI reduction operator
      struct doubleRankPair {
	 double val;
	 int rank;
      };

      // Updated in collectStats, only valid on root rank
      struct TimerStatistics {
         std::vector<int> id; //id of the timer at this index in the statistics vectors
         std::vector<int> level;
         std::vector<double> timeSum;
         std::vector<doubleRankPair> timeMax;
         std::vector<doubleRankPair> timeMin;
         std::vector<double> timeTotalFraction;
         std::vector<double> timeParentFraction;
         std::vector<bool> hasWorkUnits;
         std::vector<double> workUnitsSum;
         std::vector<int> countSum;
     	 std::vector<int> threadsSum;
      };

      struct GroupStatistics {
         std::vector<std::string> name; 
         std::vector<double> timeSum;
         std::vector<doubleRankPair> timeMax;
         std::vector<doubleRankPair> timeMin;
         std::vector<double> timeTotalFraction;
      };

     //vector with timers, cumulative and for log print
      std::vector<TimerData> _cumulativeTimers;
      std::vector<TimerData> _logTimers;
      
      //current position in timer hierarchy 
      int _currentId=-1;

      //is this profiler initialized
      bool _initialized=false;

      //used to store time when we start printing. This is used to correctly asses active timer time
      double _printStartTime;

      //defines print-area widths for print() output
      const int _indentWidth=2; //how many spaces each level is indented
      const int _floatWidth=10; //width of float fields;
      const int _intWidth=6;   //width of int fields;
      const int _unitWidth=4;  //width of workunit label
      const int _levelWidth=5; //width of level label

      //this function returns the time in seconds 
      double getTime() {
	 return MPI_Wtime();
      }
      //this function returns the accuracy of the timer     
      double getTick() {
	 return MPI_Wtick();
      }

   //get id number of a timer, return -1 if it does not exist
   int getId(const std::string &label){
      //find child with this id
      int childId=-1;
#pragma omp master
      {
         for(unsigned int i=0;i<_cumulativeTimers[_currentId].childIds.size();i++)
            if (_cumulativeTimers[_cumulativeTimers[_currentId].childIds[i]].label==label){
               childId=_cumulativeTimers[_currentId].childIds[i];
               break;
            }
      }
      return childId;
   }
   
//construct a new timer for the current level.
      int constructTimer(const std::string &label,int parentId,const std::vector<std::string> groups){
	 TimerData timerData;
	 timerData.label=label;
	 timerData.groups=groups;
	 timerData.id=_cumulativeTimers.size(); //will be added last to timers vector
	 timerData.parentId=parentId;
	 
	 timerData.workUnits=-1;
	 timerData.workUnitLabel="";
	 
	 if(parentId!=-1) 
	    timerData.level=_cumulativeTimers[parentId].level+1;
	 else //this is the special case when one adds the root timer
	    timerData.level=0;
	 timerData.time=0;
	 timerData.startTime=-1;
	 timerData.count=0;
         timerData.active=false;
	 timerData.threads=1;
#ifdef _OPENMP
	 timerData.threads=omp_get_num_threads();
#endif
	 //timerData.work  UnitCount initialized in stop
	 //add timer, to both vectors   
	 _cumulativeTimers.push_back(timerData);
         _logTimers.push_back(timerData);
	 //add timer to tree, both in _cumulativeTimers, and _logTimers
	 if(parentId!=-1){
	    _cumulativeTimers[parentId].childIds.push_back(timerData.id);
            _logTimers[parentId].childIds.push_back(timerData.id);
         }
         return timerData.id;
      }
      
         
   //start timer, with id
   bool start(int id){
      bool success=true;
#pragma omp master
      {
         if(_currentId!=_cumulativeTimers[id].parentId){
            std::cerr << "PHIPROF-ERROR: Starting timer that is not a child of the current profiling region" <<std::endl;
            success= false;
         }
         else{
            _currentId=id;      
            //start timer
            _cumulativeTimers[_currentId].startTime=getTime();
            _cumulativeTimers[_currentId].active=true;

            //start log timer
            _logTimers[_currentId].startTime=_cumulativeTimers[_currentId].startTime;
            _logTimers[_currentId].active=true;

         
#ifdef CRAYPAT
            PAT_region_begin(_currentId+1,getFullLabel(_cumulativeTimers,_currentId,true).c_str());
#endif
         }
      }
      return success;        
   }
  
      //initialize profiler, called by first start/initializeTimer. This adds the root timer
      bool init(){
	 if(!_initialized){
	    _initialized=true;
	    std::vector<std::string> group;
	    group.push_back("Total");
	    //no timer yet        
	    _currentId=-1;
	    //mainId will be 0, parent is -1 (does not exist)
	    int id=constructTimer("total",-1,group);
	    //start root timer, is stopped in print.      
	    start(id);
	 }
	 return true;
      }

   //stop a timer defined by id
   bool stop (int id,
              const double workUnits = -1.0,
              const std::string &workUnitLabel = ""){
      bool success=true;
#pragma omp master
      {
         double stopTime=getTime();
         if(id != _currentId ){
            std::cerr << "PHIPROF-ERROR: id missmatch in profile::stop Stopping "<< id <<" at level " << _cumulativeTimers[_currentId].level << std::endl;
            success=false;
         }

         else {
         
#ifdef CRAYPAT
            PAT_region_end(_currentId+1);
#endif  

         
            //handle workUnits for _cumulativeTimers               
            if(_cumulativeTimers[_currentId].count!=0){
               //if this, or a previous, stop did not include work units then do not add them
               //work units have to be defined for all stops with a certain (full)label
               if(workUnits<0 || _cumulativeTimers[_currentId].workUnits<0){
                  _cumulativeTimers[_currentId].workUnits=-1;
                  _logTimers[_currentId].workUnits=-1;
               }
               else{
                  _cumulativeTimers[_currentId].workUnits+=workUnits;
               }
            }
            else{
               //firsttime, initialize workUnit stuff here
               if(workUnits>=0.0 ){
                  //we have workUnits for this counter
                  _cumulativeTimers[_currentId].workUnits=workUnits;
                  _cumulativeTimers[_currentId].workUnitLabel=workUnitLabel;
               }
               else{
                  // no workUnits for this counter
                  _cumulativeTimers[_currentId].workUnits=-1.0;
               }
            }
            
            //handle workUnits for _logTimers
            if(_logTimers[_currentId].count!=0){
               //if this, or a previous, stop did not include work units then do not add t hem
               //work units have to be defined for all stops with a certain (full)label
               if(workUnits<0 || _logTimers[_currentId].workUnits<0){
                  _logTimers[_currentId].workUnits=-1;
               }
               else{
                  _logTimers[_currentId].workUnits+=workUnits;
               }
            }
            else{
               //firsttime, initialize workUnit stuff here
               if(workUnits>=0.0 ){
                  //we  have workUnits for this counter
                  _logTimers[_currentId].workUnits=workUnits;
                  _logTimers[_currentId].workUnitLabel=workUnitLabel;
               }
               else{
                  //  no workUnits for this counter
                  _logTimers[_currentId].workUnits=-1.0;
               }
            }
            
         //stop _cumulativeTimers & _logTimers timer              
            _cumulativeTimers[_currentId].time+=(stopTime-_cumulativeTimers[_currentId].startTime);
            _cumulativeTimers[_currentId].count++;
            _cumulativeTimers[_currentId].active=false;
            _logTimers[_currentId].time+=(stopTime-_logTimers[_currentId].startTime);
            _logTimers[_currentId].count++;
            _logTimers[_currentId].active=false;
            
      
            //go down in hierarchy    
            _currentId=_cumulativeTimers[_currentId].parentId;
         }
      }
      return success;
   }

    

  /**
   * Initialize a timer, with a particular label   
   *
   * Initialize a timer. This enables one to define groups, and to use
   * the return id value for more efficient starts/stops in tight
   * loops. If this function is called for an existing timer it will
   * simply just return the id.
   *
   *
   * @param label
   *   Name for the timer to be created. This will not yet start the timer, that has to be done separately with a start. 
   * @param groups
   *   The groups to which this timer belongs. Groups can be used combine times for different timers to logical groups, e.g., MPI, io, compute...
   * @return
   *   The id of the timer
   */
   int initializeTimer(
      const std::string &label,
      const std::vector<std::string> &groups
   ) {
      //check if the global profiler is initialized
      if(!_initialized)
         init();
      int id=getId(label); //check if it exists
      if(id>=0)
         //do nothing if it exists      
         return id; 
      else
        //create new timer if it did not exist
        return constructTimer(label,_currentId,groups);
   }


   //initialize a timer, with a particular label belonging to a group
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const std::string &label,const std::string &group){
      //check if the global profiler is initialized
      std::vector<std::string> groups;
      groups.push_back(group);
      return initializeTimer(label,groups);
   }

   //initialize a timer, with a particular label belonging to two groups
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const std::string &label,
		     const std::string &group1,
		     const std::string &group2
		     ){
      std::vector<std::string> groups;
      groups.push_back(group1);
      groups.push_back(group2);
      return initializeTimer(label,groups);
   }

   //initialize a timer, with a particular label belonging to three groups
   //returns id of new  timer. If timer exists, then that id is returned.
   int initializeTimer(const std::string &label,
		     const std::string &group1,
		     const std::string &group2,
		     const std::string &group3
		     ){
      std::vector<std::string> groups;
      groups.push_back(group1);
      groups.push_back(group2);
      groups.push_back(group3);
      return initializeTimer(label,groups);
   }


   //initialize a timer, with a particular label belonging to no group
   //returns id of new timer. If timer exists, then that id is returned.
   int initializeTimer(const std::string &label){
      std::vector<std::string> groups; //empty vector
      return initializeTimer(label,groups);
   }


  /**
   * Start a profiling timer.
   *
   * This function starts a timer with a certain label (name). If
   * timer does not exist, then it is created. The timer is
   * automatically started in the current active location in the tree
   * of timers. Thus the same start command in the code can start
   * different timers, if the current active timer is different.
   *
   * @param label
   *   Name for the timer to be start. 
   * @return
   *   Returns true if timer started successfully.
   */
   bool start(const std::string &label){
#pragma omp master
      {
         //If the timer exists, then initializeTimer just returns its id, otherwise it is constructed.
         //Make the timer the current one
         _currentId=initializeTimer(label);
         //start timer
         _cumulativeTimers[_currentId].startTime=getTime();
         _cumulativeTimers[_currentId].active=true;

         //start log timer   
         _logTimers[_currentId].startTime=_cumulativeTimers[_currentId].startTime;
         _logTimers[_currentId].active=true;
      
#ifdef CRAYPAT
         PAT_region_begin(_currentId+1,getFullLabel(_cumulativeTimers,_currentId,true).c_str());
#endif
      }
      return true;        
   }
   
  /**
   * Stop a profiling timer.
   *
   * This function stops a timer with a certain label (name). The
   * label has to match the currently last opened timer. One can also
   * (optionally) report how many workunits was done during this
   * start-stop timed segment, e.g. GB for IO routines, Cells for
   * grid-based solvers. Note, all stops for a particular timer has to
   * report workunits, otherwise the workunits will not be reported.
   *
   * @param label 
   *   Name for the timer to be stopped.     
   * @param workunits 
   *   (optional) Default is for no workunits to be
   *   collected.Amount of workunits that was done during this timer
   *   segment. If value is negative, then no workunit statistics will
   *   be collected.
   * @param workUnitLabel
   *   (optional) Name describing the unit of the workunits, e.g. "GB", "Flop", "Cells",...
   * @return
   *   Returns true if timer stopped successfully.
   */
   bool stop (const std::string &label,
              const double workUnits = -1.0,
              const std::string &workUnitLabel = ""){
      bool success=true;
#pragma omp master
      {
         if(label != _cumulativeTimers[_currentId].label ){
            std::cerr << "PHIPROF-ERROR: label missmatch in profile::stop  when stopping "<< label <<
               ". The started timer is "<< _cumulativeTimers[_currentId].label<< " at level " << _cumulativeTimers[_currentId].level << std::endl;
            success=false;
         }
         if(success)
            success=stop(_currentId,workUnits,workUnitLabel);
      }
      return success;
   }

   
////-------------------------------------------------------------------------
///  Hash functions
////-------------------------------------------------------------------------      
      
      //djb2 hash function copied from
      //http://www.cse.yorku.ca/~oz/hash.html
      unsigned long hash(const char *str)
      {
	 unsigned long hash = 5381;
	 int c;
	 while ( (c = *str++) )
	       hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
	 return hash;
      }

         
      //Hash value identifying all labels, groups and workunitlabels.
      //If any strings differ, hash should differ. Computed recursively in the same way as prints
      int getTimersHash(const std::vector<TimerData> &timers,int id=0){
         unsigned long hashValue;
         //add hash values from label, workunitlabel and groups. Everything has to match.
         hashValue=hash(timers[id].label.c_str());
         hashValue+=hash(timers[id].workUnitLabel.c_str());
         for (std::vector<std::string>::const_iterator g = timers[id].groups.begin();g != timers[id].groups.end(); ++g ) {
            hashValue+=hash((*g).c_str());
         }
         
         for(unsigned int i=0;i<timers[id].childIds.size();i++){
            hashValue+=getTimersHash(timers,timers[id].childIds[i]);
         }

         
         // MPI_Comm_split needs a non-zero value
         if (hashValue == 0) {
            return 1;
         } else {
            //we return an integer value
            return hashValue%std::numeric_limits<int>::max();
         }
      }

      //get full hierarchical name for a timer
      //can have either the timer-label first (reverse order), or last.
      std::string getFullLabel(const std::vector<TimerData> &timers,int id,bool reverse=false){
         //create a label with all hierarchical levels     
         std::vector<std::string> labels;
         while(id>0){
            labels.push_back(_cumulativeTimers[id].label);
            id=_cumulativeTimers[id].parentId;
         }
         
         std::string fullLabel;
         if(reverse){
            std::vector<std::string>::iterator it;
            for ( it=labels.begin() ; it != labels.end(); ++it ){
               fullLabel+=*it;
               fullLabel+="\\";   
            }
         }
         else{
            std::vector<std::string>::reverse_iterator rit;
            for ( rit=labels.rbegin() ; rit != labels.rend(); ++rit ){
               fullLabel+="/";
               fullLabel+=*rit;
            }
         }

         
         return fullLabel;
      }
         
////-------------------------------------------------------------------------
///  Timer handling functions functions
////-------------------------------------------------------------------------            
      


   

      
////-------------------------------------------------------------------------
///  Collect statistics functions
////-------------------------------------------------------------------------            


      double getTimerTime(int id,const std::vector<TimerData> &timers){
         double time;
         time=timers[id].time;
         //add time uptill print for active timers
         if(timers[id].active)
            time+=_printStartTime-timers[id].startTime;
         return time;
      }
   
     
      
      double getGroupInfo(std::string group,int id,const std::vector<TimerData> &timers , MPI_Comm comm){
         double groupTime=0.0;
	 for (std::vector<std::string>::const_iterator g = timers[id].groups.begin();
              g != timers[id].groups.end(); ++g ) {
	    if(group==*g){
               groupTime=getTimerTime(id,timers);
               return groupTime; // do not collect for children when this is already in group.Avoid double counting
	    }
	 }
      
         //recursively collect time data   
         for(unsigned int i=0;i<timers[id].childIds.size();i++){
            groupTime+=getGroupInfo(group,timers[id].childIds[i],timers,comm);
	 }
	 return groupTime;
      }
      

      
      void collectGroupStats(GroupStatistics &stats,const std::vector<TimerData> &timers,MPI_Comm &comm){
	 int rank,nProcesses;
         //per process info. updated in collectStats
         std::vector<double> time;
         std::vector<doubleRankPair> timeRank;
         std::map<std::string,std::vector<int> > groups;
         doubleRankPair in;
         int totalIndex=0; //where we store the group for total time (called Total, in timer id=0)
         
	 MPI_Comm_rank(comm,&rank);
	 MPI_Comm_size(comm,&nProcesses);

         
         //construct map from groups to timers in group 
	 for(unsigned int id=0;id<timers.size();id++){
	    for (std::vector<std::string>::const_iterator group = timers[id].groups.begin();
                 group != timers[id].groups.end(); ++group ) {
	       groups[*group].push_back(id);
	    }
	 }
         int nGroups=groups.size();
         stats.name.clear(); // we will use push_back to add names to this vector
         
         //collect data for groups
         for(std::map<std::string, std::vector<int> >::iterator group=groups.begin();
	    group!=groups.end();++group){
	    double groupTime=0.0;
            if(group->first=="Total")
               totalIndex=stats.name.size();

	    groupTime=getGroupInfo(group->first,0,timers,comm);
            stats.name.push_back(group->first);
            time.push_back(groupTime);
            in.val=groupTime;
            in.rank=rank;
            timeRank.push_back(in);
         }

         //Compute statistics using reduce operations
         if(rank==0){
            //reserve space for reduce operations output
            stats.timeSum.resize(nGroups);
            stats.timeMax.resize(nGroups);
            stats.timeMin.resize(nGroups);
            stats.timeTotalFraction.resize(nGroups);

            MPI_Reduce(&(time[0]),&(stats.timeSum[0]),nGroups,MPI_DOUBLE,MPI_SUM,0,comm);
            MPI_Reduce(&(timeRank[0]),&(stats.timeMax[0]),nGroups,MPI_DOUBLE_INT,MPI_MAXLOC,0,comm);
            MPI_Reduce(&(timeRank[0]),&(stats.timeMin[0]),nGroups,MPI_DOUBLE_INT,MPI_MINLOC,0,comm);

            for(int i=0;i<nGroups;i++){
               if(stats.timeSum[totalIndex]>0)
                  stats.timeTotalFraction[i]=stats.timeSum[i]/stats.timeSum[totalIndex];
               else
                  stats.timeTotalFraction[i]=0.0;
            }
         }
         else{
            //not masterank, we do not resize and use stats vectors
            MPI_Reduce(&(time[0]),NULL,nGroups,MPI_DOUBLE,MPI_SUM,0,comm);
            MPI_Reduce(&(timeRank[0]),NULL,nGroups,MPI_DOUBLE_INT,MPI_MAXLOC,0,comm);
            MPI_Reduce(&(timeRank[0]),NULL,nGroups,MPI_DOUBLE_INT,MPI_MINLOC,0,comm);            
         }

      }


//collect timer stats, call children recursively. In original code this should be called for the first id=0
// reportRank is the rank to be used in the report, not the rank in the comm communicator
      void collectTimerStats(TimerStatistics &stats,const std::vector<TimerData> &timers,MPI_Comm &comm,int reportRank,int id=0,int parentIndex=0){
         int rank;
         //per process info. updated in collectStats
         static std::vector<double> time;
         static std::vector<doubleRankPair> timeRank;
         static std::vector<double> workUnits;
         static std::vector<int> count;
         static std::vector<int> threads;
         static std::vector<int> parentIndices;
         int currentIndex;
         doubleRankPair in;

         MPI_Comm_rank(comm,&rank); 
         //first time we call  this function
         if(id==0){
            time.clear();
            timeRank.clear();
            count.clear();
            threads.clear();
	    workUnits.clear();
            parentIndices.clear();
            stats.id.clear();
            stats.level.clear();
         }
         
         //collect statistics

         currentIndex=stats.id.size();
         
         double currentTime=getTimerTime(id,timers);
         
         stats.id.push_back(id);
         stats.level.push_back(timers[id].level);
         time.push_back(currentTime);
         in.val=currentTime;
         in.rank=reportRank;
         timeRank.push_back(in);
         count.push_back(timers[id].count);
         threads.push_back(timers[id].threads);
         workUnits.push_back(timers[id].workUnits);
         parentIndices.push_back(parentIndex);
         
         double childTime=0;
         //collect data for children. Also compute total time spent in children
         for(unsigned int i=0;i<timers[id].childIds.size();i++){
            childTime+=getTimerTime(timers[id].childIds[i],timers);
            collectTimerStats(stats,timers,comm,reportRank,timers[id].childIds[i],currentIndex);
         }
         
         if(timers[id].childIds.size()>0){
            //Added timings for other time. These are assigned id=-1
            stats.id.push_back(-1);
            stats.level.push_back(timers[timers[id].childIds[0]].level); //same level as children
            time.push_back(currentTime-childTime);
            in.val=currentTime-childTime;
            in.rank=reportRank;
            timeRank.push_back(in);
            count.push_back(timers[id].count);
	    threads.push_back(timers[id].threads);
            workUnits.push_back(-1);
            parentIndices.push_back(currentIndex);
         }

         
         //End of function for id=0, we have now collected all timer data.
         //compute statistics now
         if(id==0){
            int nTimers=time.size();
            if(rank==0){
               stats.timeSum.resize(nTimers);
               stats.timeMax.resize(nTimers);
               stats.timeMin.resize(nTimers);
               stats.workUnitsSum.resize(nTimers);
               stats.hasWorkUnits.resize(nTimers);
               stats.countSum.resize(nTimers);
               stats.threadsSum.resize(nTimers);

               stats.timeTotalFraction.resize(nTimers);
               stats.timeParentFraction.resize(nTimers);
               std::vector<double> workUnitsMin;
               workUnitsMin.resize(nTimers);
               
               MPI_Reduce(&(time[0]),&(stats.timeSum[0]),nTimers,MPI_DOUBLE,MPI_SUM,0,comm);
               MPI_Reduce(&(timeRank[0]),&(stats.timeMax[0]),nTimers,MPI_DOUBLE_INT,MPI_MAXLOC,0,comm);
               MPI_Reduce(&(timeRank[0]),&(stats.timeMin[0]),nTimers,MPI_DOUBLE_INT,MPI_MINLOC,0,comm);
               
               MPI_Reduce(&(workUnits[0]),&(stats.workUnitsSum[0]),nTimers,MPI_DOUBLE,MPI_SUM,0,comm);
               MPI_Reduce(&(workUnits[0]),&(workUnitsMin[0]),nTimers,MPI_DOUBLE,MPI_MIN,0,comm);
               MPI_Reduce(&(count[0]),&(stats.countSum[0]),nTimers,MPI_INT,MPI_SUM,0,comm);
               MPI_Reduce(&(threads[0]),&(stats.threadsSum[0]),nTimers,MPI_INT,MPI_SUM,0,comm);
               
               for(int i=0;i<nTimers;i++){
                  if(workUnitsMin[i]<0)
                     stats.hasWorkUnits[i]=false;
                  else
                     stats.hasWorkUnits[i]=true;
                  
                  if(stats.timeSum[0]>0)
                     stats.timeTotalFraction[i]=stats.timeSum[i]/stats.timeSum[0];
                  else
                     stats.timeTotalFraction[i]=0.0;
                  
                  if(stats.timeSum[parentIndices[i]]>0)
                     stats.timeParentFraction[i]=stats.timeSum[i]/stats.timeSum[parentIndices[i]];
                  else
                     stats.timeParentFraction[i]=0.0;
               }
            }
            else{
               //not masterank, we do not resize and use stats vectors
               MPI_Reduce(&(time[0]),NULL,nTimers,MPI_DOUBLE,MPI_SUM,0,comm);
               MPI_Reduce(&(timeRank[0]),NULL,nTimers,MPI_DOUBLE_INT,MPI_MAXLOC,0,comm);
               MPI_Reduce(&(timeRank[0]),NULL,nTimers,MPI_DOUBLE_INT,MPI_MINLOC,0,comm);
               
               MPI_Reduce(&(workUnits[0]),NULL,nTimers,MPI_DOUBLE,MPI_SUM,0,comm);
               MPI_Reduce(&(workUnits[0]),NULL,nTimers,MPI_DOUBLE,MPI_MIN,0,comm);
               MPI_Reduce(&(count[0]),NULL,nTimers,MPI_INT,MPI_SUM,0,comm);               
	       MPI_Reduce(&(threads[0]),NULL,nTimers,MPI_INT,MPI_SUM,0,comm);
            }
            //clear temporary data structures
            time.clear();
            timeRank.clear();
            count.clear();
            threads.clear();
            workUnits.clear();
            parentIndices.clear();
         }
      }
      

//remove print time from timings by pushing forward start_time
      void removePrintTime(double endPrintTime,std::vector<TimerData> &timers,int id=0){
         if(timers[id].active){
            //push start time so that ew push it forward
            timers[id].startTime+=endPrintTime-_printStartTime;
            for(unsigned int i=0;i<timers[id].childIds.size();i++){
               removePrintTime(endPrintTime,timers,timers[id].childIds[i]);
            }
         }

      }

      //reset logtimes in timers to zero.

      void resetTime(double endPrintTime,std::vector<TimerData> &timers,int id=0){
         timers[id].time=0;
         timers[id].count=0;
         timers[id].workUnits=0;
                     
         if(timers[id].active){
            timers[id].startTime=endPrintTime;
         }

         for(unsigned int i=0;i<timers[id].childIds.size();i++){
            resetTime(endPrintTime,timers,timers[id].childIds[i]);
         }
      }            
      



////-------------------------------------------------------------------------
///  PRINT functions

////-------------------------------------------------------------------------      

      // Creating map of group names to group one-letter ID
      // We assume same timers exists in all timer vectors, so can use just one here
      void getGroupIds(std::map<std::string, std::string> &groupIds,size_t &groupWidth,const std::vector<TimerData> &timers){
         groupIds.clear();
         groupWidth=6;
         //add groups to map
         for(unsigned int id=0;id<timers.size();id++) {
            size_t width = timers[id].groups.size();
            groupWidth = std::max(width, groupWidth);               
            for (std::vector<std::string>::const_iterator group = timers[id].groups.begin();
                 group != timers[id].groups.end(); ++group ) {
               groupIds[*group] = *group;
            }
         }
         //assing letters
         int character = 65; // ASCII A, TODO skip after Z to a, see ascii table
         for(std::map<std::string, std::string>::const_iterator group = groupIds.begin();
             group != groupIds.end(); ++group) {
            groupIds[group->first] = character++;
         }
      
      }


      //print all timers in stats
      bool printTreeTimerStatistics(const TimerStatistics &stats,const std::vector<TimerData> &timers,double minFraction,size_t labelWidth,size_t groupWidth,int totalWidth,
                                    const std::map<std::string,std::string> &groupIds,int nProcesses,std::fstream &output){

         for(int i=0;i<totalWidth/2-5;i++) output <<"-";
         output << " Profile ";
         for(int i=0;i<totalWidth/2-5;i++) output <<"-";
         output<<std::endl;

         for(unsigned int i=1;i<stats.id.size();i++){
            int id=stats.id[i];
            if(stats.timeTotalFraction[i]>=minFraction){
               //print timer if enough time is spent in it
               bool hasNoGroups = true;
               int indent=(stats.level[i]-1)*_indentWidth;
               output << std::setw(_levelWidth+1) << stats.level[i];
               
               if(id!=-1){
		 //other label has no groups
                  for (std::vector<std::string>::const_iterator group = timers[id].groups.begin();
                       group != timers[id].groups.end(); ++group) {
		    std::string groupId=groupIds.count(*group) ? groupIds.find(*group)->second : std::string();
		    output << std::setw(1) << groupId;
		    hasNoGroups = false;
                  }
               }
               
               if(hasNoGroups) output << std::setw(groupWidth+1) << "";
               else output << std::setw(groupWidth-timers[id].groups.size()+1) << "";
	       
               output << std::setw(indent) << "";
               if(id!=-1){
                  output << std::setw(labelWidth+1-indent) << std::setiosflags(std::ios::left) << timers[id].label;
               }
               else{
                  output << std::setw(labelWidth+1-indent) << std::setiosflags(std::ios::left) << "Other";
               }

               output << std::setw(_floatWidth+1) << stats.threadsSum[i]/nProcesses;        	       
               output << std::setw(_floatWidth) << stats.timeSum[i]/nProcesses;
               output << std::setw(_floatWidth) << 100.0*stats.timeParentFraction[i];
               output << std::setw(_floatWidth) << stats.timeMax[i].val;
               output << std::setw(_intWidth)   << stats.timeMax[i].rank;
               output << std::setw(_floatWidth) << stats.timeMin[i].val;
               output << std::setw(_intWidth)   << stats.timeMin[i].rank;
               output << std::setw(_floatWidth) << stats.countSum[i]/nProcesses;

               if(stats.hasWorkUnits[i]){
                     //print if units defined for all processes
                     //note how the total rate is computed. This is to avoid one process with little data to     
                     //skew results one way or the other                        
                  if(stats.timeSum[i]>0){
                     output << std::setw(_floatWidth) << nProcesses*(stats.workUnitsSum[i]/stats.timeSum[i]);
                     output << std::setw(_floatWidth) << stats.workUnitsSum[i]/stats.timeSum[i];
                  }
                  else if (stats.workUnitsSum[i]>0){
                     //time is zero
                     output << std::setw(_floatWidth) << "inf";
                     output << std::setw(_floatWidth) << "inf";
                  }
                  else {
                     //zero time zero units
                     output << std::setw(_floatWidth) << 0;
                     output << std::setw(_floatWidth) << 0;
                  }
                  output << timers[id].workUnitLabel<<"/s";                     
               }
               output<<std::endl;
            }
            
	 }
	 return true;
      }
      
      
      bool printTreeFooter(int totalWidth,std::fstream &output){
         for(int i=0;i<totalWidth;i++) output <<"-";
         output<<std::endl;
	 return true;
      }
      
      bool printTreeHeader(double minFraction,size_t labelWidth,size_t groupWidth,int totalWidth,int nProcs,std::fstream &output){
         for(int i=0;i<totalWidth;i++) output <<"-";
         output<<std::endl;
         output << "Phiprof results with time fraction of total time larger than " << minFraction;
         output<<std::endl;
         output << "Processes in set of timers " << nProcs;
#ifdef _OPENMP
	 output << " with (up to) " << omp_get_max_threads() << " threads ";
#endif 
         output<<std::endl;
         for(int i=0;i<totalWidth;i++) output <<"-";
         output<<std::endl;
         output<<std::setw(_levelWidth+1+groupWidth+1+labelWidth+1)<< std::setiosflags(std::ios::left) << "";
         output<<std::setw(_floatWidth)<< "Threads";
         output<<std::setw(4*_floatWidth+2*_intWidth) <<"Time(s)";
         output<<std::setw(_floatWidth)<<"Calls";
         output<<std::setw(2*_floatWidth)<<"Workunit-rate";
         output<<std::endl;
         output<<std::setw(_levelWidth+1)<< "Level";	    
         output<<std::setw(groupWidth+1)<< "Groups";
//         output << std::setw(1) << "|";
         output<<std::setw(labelWidth+1)<< "Label";
//         output << std::setw(1) << "|";
	    //  time
         output<<std::setw(_floatWidth) <<"Average";
         output<<std::setw(_floatWidth) <<"Average";
         output<<std::setw(_floatWidth) <<"parent %";
         output<<std::setw(_floatWidth) <<"Maximum";
         output<<std::setw(_intWidth) << "Rank";
         output<<std::setw(_floatWidth)<< "Minimum";
         output<<std::setw(_intWidth) << "Rank";
         //call count
         output<<std::setw(_floatWidth) << "Average";
         // workunit rate    
         output<<std::setw(_floatWidth) << "Average";
         output<<std::setw(_floatWidth) << "Per process";
         output<<std::endl;
      
	 return true;
      }

            //print groups
      bool printTreeGroupStatistics(const GroupStatistics &groupStats,
                                double minFraction,
                                size_t labelWidth,
                                size_t groupWidth,
                                int totalWidth,
                                const std::map<std::string, std::string> &groupIds,
                                int nProcesses,
                                std::fstream &output){

         for(int i=0;i<totalWidth/2 -4;i++) output <<"-";
         output <<" Groups ";
         for(int i=0;i<totalWidth/2 -3;i++) output <<"-";
         output<<std::endl;

         for(unsigned int i=0;i<groupStats.name.size();i++){            
            if(minFraction<=groupStats.timeTotalFraction[i]){
	      std::string groupId= groupIds.count(groupStats.name[i]) ? groupIds.find(groupStats.name[i])->second : std::string();
	      output << std::setw(_levelWidth+1) << " ";
	      output << std::setw(groupWidth+1) << groupId;
	      output << std::setw(labelWidth+1) << groupStats.name[i];
	      output << std::setw(_floatWidth) << " ";
	      output << std::setw(_floatWidth) << groupStats.timeSum[i]/nProcesses;
	      output << std::setw(_floatWidth) << 100.0*groupStats.timeTotalFraction[i];
	      output << std::setw(_floatWidth) << groupStats.timeMax[i].val;
	      output << std::setw(_intWidth)   << groupStats.timeMax[i].rank;
	      output << std::setw(_floatWidth) << groupStats.timeMin[i].val;
	      output << std::setw(_intWidth)   << groupStats.timeMin[i].rank;
	      output << std::endl;
	    }
         }

         return true;
      }
         
      
      
      //print out global timers
      //If any labels differ, then this print will deadlock. Only call it with a communicator that is guaranteed to be consistent on all processes.
      bool printTree(const TimerStatistics &stats,const GroupStatistics &groupStats,const std::vector<TimerData> &timers,double minFraction,std::string fileName, MPI_Comm comm){
         int rank,nProcesses;
         std::map<std::string, std::string> groupIds;
         size_t labelWidth=0;    //width of column with timer labels
         size_t groupWidth=0;    //width of column with group letters
         int totalWidth=6;       //total width of the table
      
         MPI_Comm_rank(comm,&rank);
         MPI_Comm_size(comm,&nProcesses);
         
         //compute labelWidth
         if(rank==0){
            std::fstream output;
            output.open(fileName.c_str(), std::fstream::out);
            if (output.good() == false)
               return false;

            
            for (unsigned int i=0;i<timers.size();i++){
               size_t width=timers[i].label.length()+(timers[i].level-1)*_indentWidth;
               labelWidth=std::max(labelWidth,width);
            }

            getGroupIds(groupIds,groupWidth,timers);
	    
            //make sure we use default floats, and not fixed or other format
            output <<std::resetiosflags( std::ios::floatfield );
            //set       float p  rec  ision
            output <<std::setprecision(_floatWidth-6); //6 is needed for ".", "e+xx" and a space
            
            unsigned int labelDelimeterWidth=2;
            totalWidth=_levelWidth+1+groupWidth+1+labelWidth+1+labelDelimeterWidth+_floatWidth*7+_intWidth*2+_unitWidth;
         
            //print header 
            printTreeHeader(minFraction,labelWidth,groupWidth,totalWidth,nProcesses,output);
            //print out all labels recursively
            printTreeTimerStatistics(stats,timers,minFraction,labelWidth,groupWidth,totalWidth,groupIds,nProcesses,output);
            //print groups
            printTreeGroupStatistics(groupStats,minFraction,labelWidth,groupWidth,totalWidth,groupIds,nProcesses,output);
            //print footer  
            printTreeFooter(totalWidth,output);
            // start root timer again in case we continue and call print several times
            output.close();
            
         }
         return true;
      }
      
      
      //print out global timers and groups in log format

      bool printLog(const TimerStatistics &stats,const GroupStatistics &groupStats,const std::vector<TimerData> &timers,std::string fileName, std::string separator,double simulationTime,int maxLevel,MPI_Comm comm){
         int rank,nProcesses;
         MPI_Comm_rank(comm,&rank);
         MPI_Comm_size(comm,&nProcesses);
         
         //compute labelWidth
         if(rank==0){
            std::fstream output;
            bool fileExists;
            int lineWidth=80;
            //test if file exists
            output.open(fileName.c_str(), std::fstream::in);
            fileExists=output.good();
            output.close();

            //FIXME, we should instead check if this profiler haswritten out anything and overwrite old files
            
            if(fileExists)
               output.open(fileName.c_str(), std::fstream::app|std::fstream::out);
            else
               output.open(fileName.c_str(), std::fstream::out);

            if(!fileExists){
               int column=1;
               size_t groupWidth;
               std::map<std::string, std::string> groupIds;               
               getGroupIds(groupIds,groupWidth,timers);
               
               //write header
               output << "#";
               for(int i=0;i<lineWidth;i++) output << "-";
               output << std::endl;
               output << "# Profile with "<< nProcesses <<" number of processes"<<std::endl;
               output << "# In first column of header, the first column for that timer/group is printed"<<std::endl;
               output << "# Each timer has data for count, time and workunits, 9 columns in total."<<std::endl;
               output << "#       +0  Count average"<<std::endl;
               output << "#       +1  Time  average"<<std::endl;
               output << "#       +2        percent of parent time"<<std::endl;
               output << "#       +3        max value"<<std::endl;
               output << "#       +4        max rank"<<std::endl;
               output << "#       +5        min value"<<std::endl;
               output << "#       +6        min rank"<<std::endl;
               output << "#       +7  Workunits (/s) average total"<<std::endl;
               output << "#       +8                 average per process"<<std::endl;
               output << "# Each group has data for time, 6 columns in total"<<std::endl;
               output << "#       +0  Time  average"<<std::endl;
               output << "#       +1        percent of total time"<<std::endl;
               output << "#       +2        max value"<<std::endl;
               output << "#       +3        max rank"<<std::endl;
               output << "#       +4        min value"<<std::endl;
               output << "#       +5        min rank"<<std::endl;
               for(int i=0;i<lineWidth/2 -4;i++) output <<"-";
               output <<" User data ";
               for(int i=0;i<lineWidth/2 -3;i++) output <<"-";
               output<<std::endl;
               output <<"#";
               output << std::setw(6) << column;
               output << std::setw(groupWidth+1) << "";
               output << "Simulation time";
               output << std::endl;
               column+=1;    

               
               output << "#";               
               for(int i=0;i<lineWidth/2 -4;i++) output <<"-";
               output <<" Groups ";
               for(int i=0;i<lineWidth/2 -3;i++) output <<"-";
               output<<std::endl;
               
               for(unsigned int i=0;i<groupStats.name.size();i++){            
		 std::string groupId=groupIds.count(groupStats.name[i]) ? groupIds.find(groupStats.name[i])->second : std::string();
                  output <<"#";
                  output << std::setw(6) << column;
                  output << std::setw(2) << groupId<<" ";
                  output << groupStats.name[i];
                  output << std::endl;
                  column+=6;
               }

               output << "#";               
               for(int i=0;i<lineWidth/2 -4;i++) output <<"-";
               output <<" Timers ";
               for(int i=0;i<lineWidth/2 -3;i++) output <<"-";
               output<<std::endl;

               
               for(unsigned int i=1;i<stats.id.size();i++){
                  int id=stats.id[i];
                  bool hasNoGroups = true;
                  int indent=(stats.level[i]-1)*_indentWidth;
                  if(stats.level[i]<=maxLevel || maxLevel<=0){
                     output <<"#";
                     output << std::setw(6) << column <<" ";
                     if(id!=-1){
                        //other label has no groups
                        for (std::vector<std::string>::const_iterator group = timers[id].groups.begin();
                             group != timers[id].groups.end(); ++group) {
			  std::string groupId=groupIds.count(*group) ? groupIds.find(*group)->second : std::string();
			  output << std::setw(2) << groupId;
			  hasNoGroups = false;
                        }
                     }
                     
                     if(hasNoGroups) output << std::setw(groupWidth+1) << "";
                     else output << std::setw(groupWidth-timers[id].groups.size()+1) << "";
                     
                     output << std::setw(indent) << "";
                     if(id!=-1){
                        output   << timers[id].label;
                     }
                     else{
                        output  << "Other";
                     }
                     output<<std::endl;
                     column+=9;
                  }
               }
               
               output << "#";
               for(int i=0;i<lineWidth;i++) output <<"-";
               output << std::endl;

               //end header write
            }

            output << simulationTime<<separator;
            
            for(unsigned int i=0;i<groupStats.name.size();i++){            
               output << groupStats.timeSum[i]/nProcesses<<separator;
               output << 100.0*groupStats.timeTotalFraction[i]<<separator;
               output << groupStats.timeMax[i].val<<separator;
               output << groupStats.timeMax[i].rank<<separator;
               output << groupStats.timeMin[i].val<<separator;
               output << groupStats.timeMin[i].rank<<separator;
             }

            for(unsigned int i=1;i<stats.id.size();i++){
               if(stats.level[i]<=maxLevel || maxLevel<=0){
                  output <<  stats.countSum[i]/nProcesses<<separator;
                  output <<  stats.timeSum[i]/nProcesses<<separator;
                  output <<  100.0*stats.timeParentFraction[i]<<separator;
                  output <<  stats.timeMax[i].val<<separator;
                  output <<  stats.timeMax[i].rank<<separator;
                  output <<  stats.timeMin[i].val<<separator;
                  output <<  stats.timeMin[i].rank<<separator;
                  
                  
                  if(stats.hasWorkUnits[i]){
                     //print if units defined for all processes
                     //note how the total rate is computed. This is to avoid one process with little data to     
                     //skew results one way or the other                        
                     if(stats.timeSum[i]>0){
                        output <<nProcesses*(stats.workUnitsSum[i]/stats.timeSum[i])<<separator;
                        output << stats.workUnitsSum[i]/stats.timeSum[i]<<separator;
                     }
                     
                     else if (stats.workUnitsSum[i]>0){
                        //time is zero
                        output << "inf"<<separator;
                        output << "inf"<<separator;
                     }
                     else{
                        //zero time zero units, or no time un
                        output << 0<<separator;
                        output << 0<<separator;
                     }
                  }
                  else{
                     //no workunits;
                     output << "nan" <<separator<< "nan" <<separator;
                  }
               } //if level<=maxlevel
            }
            output<<std::endl;
            output.close();

         }
            
         return true;
      }

      
      bool getPrintCommunicator(int &printIndex,int &timersHash,MPI_Comm &printComm, MPI_Comm comm){
         int mySuccess=1;
         int success;
         int myRank;
         MPI_Comm_rank(comm,&myRank);

//hash should be the same for _cumulativeTimers, _logTimers

         timersHash=getTimersHash(_cumulativeTimers);

         
         int result = MPI_Comm_split(comm, timersHash, 0, &printComm);
         
         if (result != MPI_SUCCESS) {
            int error_string_len = MPI_MAX_ERROR_STRING;
            char error_string[MPI_MAX_ERROR_STRING + 1];
            MPI_Error_string(result, error_string, &error_string_len);
            error_string[MPI_MAX_ERROR_STRING] = 0;
            if(myRank==0)
               std::cerr << "PHIPROF-ERROR: Error splitting communicator for printing: " << error_string << std::endl;
            mySuccess=0;
         }
         
         //Now compute the id for the print comms. First communicatr is given index 0, next 1 and so on.
         int printCommRank;
         int printCommSize;
         //get rank in printComm
         MPI_Comm_rank(printComm,&printCommRank);
         MPI_Comm_size(printComm,&printCommSize);
         //communicator with printComm masters(rank=0), this will be used to number the printComm's
         MPI_Comm printCommMasters;
         MPI_Comm_split(comm,printCommRank==0,-printCommSize,&printCommMasters);
         MPI_Comm_rank(printCommMasters,&printIndex);
         MPI_Comm_free(&printCommMasters);
         MPI_Bcast(&printIndex,1,MPI_INT,0,printComm);

         //check that the hashes at least have the same number of timers, just to be sure and to avoid crashes...

         int nTimers=_cumulativeTimers.size();
         int maxTimers;
         int minTimers;
         
         MPI_Reduce(&nTimers,&maxTimers,1,MPI_INT,MPI_MAX,0,printComm);
         MPI_Reduce(&nTimers,&minTimers,1,MPI_INT,MPI_MIN,0,printComm);

         if(printCommRank==0) { 
            if(minTimers !=  maxTimers) {
               std::cerr << "PHIPROF-ERROR: Missmatch in number of timers, hash conflict?  maxTimers = " << maxTimers << " minTimers = " << minTimers << std::endl;
               mySuccess=0;
            }
         }
         
         //if any process failed, the whole routine failed
         MPI_Allreduce(&mySuccess,&success,1,MPI_INT,MPI_MIN,comm);
         
         return success;

      }

   
  /**
   * Print the  current timer state in a human readable file
   *
   * This function will print the timer statistics in a text based
   * hierarchical form into file(s), each unique set of hierarchical
   * profiles (labels, hierarchy, workunits) will be written out to a
   * separate file. This function will print the times since the
   * ininitalization of phiprof in the first start call. It can be
   * called multiple times, and will not close currently active
   * timers. The time spent in active timers uptill the print call is
   * taken into account, and the time spent in the print function will
   * be corrected for in them.
   *
   *
   * @param comm
   *   Communicator for processes that print their profile.
   * @param fileNamePrefix
   *   (optional) Default value is "profile"
   *   The first part of the filename where the profile is printed. Each
   *   unique set of timers (label, hierarchy, workunits) will be
   *   assigned a unique hash number and the profile will be written
   *   out into a file called fileprefix_hash.txt
   * @param minFraction
   *   (optional) Default value is to print all timers
   *   minFraction can be used to filter the timers being printed so
   *   that only the ones with a meaningfull amount of time are
   *   printed. Only timers with (timer time)/(total time)>=minFraction
   *   are printed. If minfraction is <=0.0 then all timers are printed.
   * @return
   *   Returns true if pofile printed successfully.
   */

   bool print(
      MPI_Comm comm,
      std::string fileNamePrefix = "profile",
      double minFraction = 0.0
   ) {
#pragma omp master
      {
         int timersHash,printIndex;
         int rank,nProcesses;
         MPI_Comm printComm;
         TimerStatistics stats;
         GroupStatistics groupStats;
      
         //_printStartTime defined in namespace, used to correct timings for open timers
         _printStartTime=getTime();
         MPI_Barrier(comm);
         
         MPI_Comm_rank(comm,&rank);
         MPI_Comm_size(comm,&nProcesses);
         
         //get hash value of timers and the print communicator
         if(getPrintCommunicator(printIndex,timersHash,printComm,comm)) {
            //generate file name
            std::stringstream fname;
            fname << fileNamePrefix << "_" << printIndex << ".txt";
            collectTimerStats(stats,_cumulativeTimers,printComm,rank); //rank is in respect to comm given to print!
            collectGroupStats(groupStats,_cumulativeTimers,printComm);
            printTree(stats,groupStats,_cumulativeTimers,minFraction,fname.str(),printComm);
         }
         
         MPI_Comm_free(&printComm);
         MPI_Barrier(comm);
         
         double endPrintTime=getTime();
         removePrintTime(endPrintTime,_cumulativeTimers);
         removePrintTime(endPrintTime,_logTimers);
      }
      return true;
   }
   
    /**
   * Print the current timer state in a easily parsable format
   *
   * This function will print the timer statistics in a text based
   * parsable format, each unique set of hierarchical profiles
   * (labels, hierarchy, workunits) will be written out to a separate
   * file. This function will print the times since the last call to
   * printLogProfile, or since initialization for the first call to
   * printLogProfile. It is meant to be called multiple times, to
   * track the develpoment of performance metrics. 
   *
   *
   * @param comm
   *   Communicator for processes that print their profile.
   * @param fileNamePrefix
   *   (optional) Default value is "profile_log"
   *   The first part of the filename where the profile is printed. Each
   *   unique set of timers (label, hierarchy, workunits) will be
   *   assigned a unique hash number and the profile will be written
   *   out into a file called fileprefix_hash.txt
   * @param separator
   *   (optional) Default value is " "
   *   The separator between fields in the file.
   * @param maxLevel
   *   (optional) Default value is to print all timers
   *   Maxlevel can be used to limit the number of timers that are
   *   printed out. Only timers with a level in the
   *   hierarchy<=maxLevel are printed.
   * @return
   *   Returns true if pofile printed successfully.
   */
   bool printLogProfile(
      MPI_Comm comm,
      double simulationTime,
      std::string fileNamePrefix = "profile_log",
      std::string separator = " ",
      int maxLevel = 0
   ) {
#pragma omp master
      {
         int timersHash,printIndex;
         int rank,nProcesses;
         MPI_Comm printComm;
         TimerStatistics stats;
         GroupStatistics groupStats;
         
         //_printStartTime defined in namespace, used to correct timings for open timers
         _printStartTime=getTime();
         MPI_Barrier(comm);
         
         MPI_Comm_rank(comm,&rank);
         MPI_Comm_size(comm,&nProcesses);
         
         //get hash value of timers and the print communicator
         if(getPrintCommunicator(printIndex,timersHash,printComm,comm)) {
         
            //generate file name, here we use timers hash to avoid overwriting old data!
            std::stringstream fname;
            fname << fileNamePrefix << "_" << timersHash << ".txt";
            
            //collect statistics
            collectTimerStats(stats,_logTimers,printComm,rank);
            collectGroupStats(groupStats,_logTimers,printComm);
            //print log
            printLog(stats,groupStats,_logTimers,fname.str(),separator,simulationTime,maxLevel,printComm);
         }
         MPI_Comm_free(&printComm);
         MPI_Barrier(comm);
         
         double endPrintTime=getTime();
         removePrintTime(endPrintTime,_cumulativeTimers);
         resetTime(endPrintTime,_logTimers);
      }
      return true;
   }

   /**
   * Assert function
   *
   * It prints out the error message, the suplied line and file
   * information. It also prints out the current position in the timer
   * stack in phiprof. Finally it terminates the process, and thus the
   * whole parallel program. Can also be used through the
   * phiprof_assert macro
   *
   * @param condition
   *   If false, then the error is raised.
   * @param error_message
   *   The error message
   * @param file
   *   The source  file where this is called, typically supplied by  __FILE__
   * @param line     
   *   The line in the source  file where this is called, typically supplied by  __LINE__ 
   */
   void phiprofAssert(
      bool condition,
      const std::string& error_message,
      const std::string& file = "",
      int line = 0
   ) {
#ifndef NDEBUG
      if(!condition) {
#pragma omp critical
         {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            std::cerr << "ASSERT ERROR on process "<< rank << ": " << error_message      
                 << ", File: " << file << ", Line: " << line
                 << ", phiproftimer:" << getFullLabel(_cumulativeTimers,_currentId,false)
                 << std::endl;       
            exit(1);
         }
      }
#endif
      return;
   }

}; // class Phiprof
} // namespace phiprof


#endif
