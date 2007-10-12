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

#include "VrpHeurTSP.h"
#include "VrpModel.h"

//#############################################################################

// Given current LP solution. Rouding number closed to 1 to 1
// Let those edge be in the route, then use nearest neighbor idea
// to construct a tour.
bool
VrpHeurTSP::searchSolution(double & objectiveValue, double * newSolution)
{
    bool feasible = false;
    


    
    
}

//#############################################################################

void VrpHeurTSP::findNearest()
{
    


    
}

//#############################################################################
