/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
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
