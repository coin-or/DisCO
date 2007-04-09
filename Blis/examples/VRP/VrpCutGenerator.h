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

#ifndef VrpCutGenerator_h_
#define VrpCutGenerator_h_

//#############################################################################

#include "BlisConGenerator.h"

//#############################################################################

class VrpCutGenerator : public BlisConGenerator  
{
private:
    
public:
    
    /** Default construtor. */
    VrpCutGenerator() {}

    /** Destructor. */
    virtual ~VrpCutGenerator() {}
    
    /** Generate cons for the client model.

	Evaluate the state of the client model and decide whether to 
	generate cons. The generated cons are inserted into and returned 
	in the collection of cons \p cs.
	
	If \p fullScan is true, the generator is obliged to call the CGL
	\c generateCuts routine.  Otherwise, it is free to make a local 
	decision. The current implementation uses \c strategy_ to decide.
        
	The routine returns true if reoptimisation is needed (because the 
	state of the solver interface has been modified).
    */
    virtual bool generateCons(OsiCuts &cs, bool fullScan);
};

//#############################################################################

#endif
