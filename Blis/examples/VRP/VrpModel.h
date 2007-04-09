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

#include "BlisModel.h"

//#############################################################################

/** Model class for VRP. */
class VrpModel : public BlisModel 
{    
private:
    
public:

    /** Default construtor. */
    VrpModel() {}

    /** Destructor. */
    virtual ~VrpModel() {}

    /** Read in the instance data. */
    virtual void readInstance(const char* dateFile);

    /** User's criteria for a feasible solution. Return true (feasible)
	or false (infeasible) 
    */
    virtual bool userFeasibleSolution() { return true; }
    
};

//#############################################################################

#endif
