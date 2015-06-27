/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "DcoBranchStrategyMaxInf.hpp"
#include "DcoObjectInt.hpp"

//#############################################################################

// Copy constructor 
DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf (
    const DcoBranchStrategyMaxInf & rhs)
    : BcpsBranchStrategy()
{
    bestChangeUp_ = rhs.bestChangeUp_;
    bestNumberUp_ = rhs.bestNumberUp_;
    bestChangeDown_ = rhs.bestChangeDown_;
    bestNumberDown_ = rhs.bestNumberDown_;
    type_ = rhs.type_;
}

//#############################################################################

/** Create a set of candidate branching objects. */
int 
DcoBranchStrategyMaxInf::createCandBranchObjects(int numPassesLeft,
						  double ub)
{

    int numInfs = 0;
    
    int i, col, preferDir, maxInfDir = 0, maxScoreDir = 0;
    
    double score, maxScore = 0.0;
    double infeasibility, maxInf = 0.0;
    
    DcoModel *model = dynamic_cast<DcoModel *>(model_);
    
    DcoObjectInt * intObject = 0;
    DcoObjectInt * maxInfIntObject = 0;
    DcoObjectInt * maxScoreIntObject = 0;
    
    int numObjects = model->numObjects();
    
    double *objCoef = model->getObjCoef();
    
    for (i = 0; i < numObjects; ++i) {
	
        // TODO: currently all integer object.
        intObject = dynamic_cast<DcoObjectInt *>(model->objects(i));
        infeasibility = intObject->infeasibility(model, preferDir);
        
        if (infeasibility) {
            ++numInfs;
            
            if (infeasibility > maxInf) {
                maxInfIntObject = intObject;
                maxInfDir = preferDir;
                maxInf = infeasibility;
            }
            
            col = intObject->columnIndex();
            score = ALPS_FABS(objCoef[col] * infeasibility);
            
            if (score > maxScore) {
                maxScoreIntObject = intObject;
                maxScoreDir = preferDir;
                maxScore = score;
            }
        }
    }

    assert(numInfs > 0);
    
    if (maxScoreIntObject) {
        maxInfIntObject = maxInfIntObject;
        maxInfDir = maxScoreDir;
    }
    
    numBranchObjects_ = 1;
    branchObjects_ = new BcpsBranchObject* [1];
    branchObjects_[0] = maxInfIntObject->createBranchObject(model,
                                                            maxInfDir);
    
    return 0;
}

//#############################################################################

/** Compare branching object thisOne to bestSoFar. If thisOne is better 
    than bestObject, return branching direction(1 or -1), otherwise
    return 0. 
    If bestSorFar is NULL, then always return branching direction(1 or -1).
*/
int
DcoBranchStrategyMaxInf::betterBranchObject(BcpsBranchObject * thisOne,
					     BcpsBranchObject * bestSoFar)
{
    return thisOne->getDirection();
}

//#############################################################################
