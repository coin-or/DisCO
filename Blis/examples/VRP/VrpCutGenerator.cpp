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

#include "BlisModel.h"
#include "VrpCutGenerator.h"

//#############################################################################

bool 
VrpCutGenerator::generateCons(OsiCuts &cs, bool fullScan)
{
    // Generator cuts for VRP

    // Can get LP solution information from solver;
    OsiSolverInterface * solver = model_->solver();

    
    return false;
}

//#############################################################################
