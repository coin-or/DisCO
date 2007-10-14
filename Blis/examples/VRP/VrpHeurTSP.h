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

#ifndef VrpHeurTSP_h_
#define VrpHeurTSP_h_

//#############################################################################

#include <vector>

#include "CoinPackedVector.hpp"

#include "BlisHeuristic.h"
#include "VrpModel.h"

//#############################################################################
#if 0
class VrpAdjList 
{
private:
    std::vector<CoinPackedVector *> list_;
    VrpModel *model_;
public:
    VrpAdjList() {}
    VrpAdjList(VrpModel *m) { 
	model_ = m; 
	createAdjList(m);
    }
    
    void createAdjList(VrpModel *model);
};
#endif

//#############################################################################

class VrpHeurTSP : public BlisHeuristic {
private:
    /** Illegal Assignment operator. */ 
    VrpHeurTSP & operator=(const VrpHeurTSP& rhs);
    
protected:
    /* Stored the predetermined next vertex to visit for vertex k if 
       the value determined_[k] greater than zero. */
    //int *determined_;
    
    /* Adjacent list of all vertices. */
    std::vector<CoinPackedVector *> adjList_;

    /** Create adjacent list for each vertex. */
    void createAdjList(VrpModel *model);

public:
    /** Default Constructor. */
    VrpHeurTSP() {}
    
    /** Constructor with model. */
    VrpHeurTSP(VrpModel * model, const char *name,
               BlisHeurStrategy strategy, int freq)
        :
        BlisHeuristic(model, name, strategy, freq) 
    { 
        createAdjList(model);
    }

    /** Destructor. */
    ~VrpHeurTSP() 
    { 
        //delete [] determined_; 
    }
    
    /** Returns 0 if no solution, 1 if valid solution. newSolution stores
        the solution in dense format. */
    virtual bool searchSolution(double & objectiveValue, double * newSolution);
};
#endif

//#############################################################################
