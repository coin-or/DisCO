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

#include "VrpModel.h"

//#############################################################################
void
VrpModel::readInstance(const char* dateFile)
{
    
    
}

//#############################################################################

void 
VrpModel::loadProblem(int numVars,
		      int numCons,
		      BcpsVariable **variable,
		      double *conLower,
		      double *conUpper)
{
    VrpVariable * var = NULL;
    int k;

    double* collb = new double [numVars];
    double* colub = new double [numVars];
    double* obj = new double [numVars];

    CoinBigIndex * start = new CoinBigIndex [numVars+1];
    int* index = NULL;
    double* value = NULL;
    
    // Get collb, colub, obj, and matrix from variables
    for (k = 0; k < numVars; ++k) {
	var = dynamic_cast<VrpVariable *>(variable[k]);
    }

    // Load to lp solver
    origLpSolver_->loadProblem(numVars, numCons,
			       start, index, value,
			       collb, colub, obj,
			       conLower, conUpper);
    
    origLpSolver_->setObjSense(objSense_);

    origLpSolver_->setInteger(intColIndices_, numIntObjects_);    
    
}

//#############################################################################
