/*===========================================================================*
 * This file is part of a solver for the Vehicle Routing Problem             *
 * developed using the BiCePS Linear Integer Solver (BLIS).                  *
 *                                                                           *
 * This solver is distributed under the Common Public License as part of     * 
 * the COIN-OR repository (http://www.coin-or.org).                          *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2007 Yan Xu and Ted Ralphs.                                 *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef VrpModel_h_
#define VrpModel_h_

//#############################################################################

#include <vector>

#include "BlisModel.h"
#include "VrpVariable.h"
#include "VrpCommonTypes.h"

//#############################################################################

/** Model class for VRP. */
class VrpModel : public BlisModel 
{    
private:

   char name_[100];
   int vertnum_;
   int edgenum_;
   int numroutes_;
   int depot_;
   int capacity_;
   int wtype_;
   int *demand_;
   int *posx_;
   int *posy_;
   double *coordx_;
   double *coordy_;
   double *coordz_;
   std::vector<VrpVariable *> edges_;
   std::vector<best_tours> curTour_;
   
public:

    /** Default construtor. */
    VrpModel() : vertnum_(0), edgenum_(0), numroutes_(0), depot_(0),
       capacity_(0), wtype_(0){
       demand_ = 0;
       posx_ = 0;
       posy_ = 0;
       coordx_ = 0;
       coordy_ = 0;
       coordz_ = 0;
    }

    VrpModel(const char* dataFile){
       readInstance(dataFile);
    }
    
    /** Destructor. */
    virtual ~VrpModel() {}

    /** Read in the instance data. */
    virtual void readInstance(const char* dateFile);

    /** User's criteria for a feasible solution. Return true (feasible)
     *	or false (infeasible) 
     */
    virtual bool userFeasibleSolution() { return true; }

    int index (int v0, int v1) {
       return(v0 < v1 ? v1*(v1 - 1)/2 + v0 : v0*(v0 - 1)/2 + v1);
    }

    int compute_cost(int v0, int v1); 
   
    /** 1) Load problem to lp solver. Assume lp solver is ready. 
     *  2) Set objective sense
     *  3) Set integer
     */
    void loadProblem(int numVars,
		     int numCons,
		     std::vector<VrpVariable *> vars);
 
};

//#############################################################################

#endif
