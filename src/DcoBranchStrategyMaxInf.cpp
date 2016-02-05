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
    : BcpsBranchStrategy() {
  bestChangeUp_ = rhs.bestChangeUp_;
  bestNumberUp_ = rhs.bestNumberUp_;
  bestChangeDown_ = rhs.bestChangeDown_;
  bestNumberDown_ = rhs.bestNumberDown_;
  type_ = rhs.type_;
}

//#############################################################################

/** Create a set of candidate branching objects. */
int DcoBranchStrategyMaxInf::createCandBranchObjects(int numPassesLeft,
						     double ub) {
  int numInfs = 0;
  // max inf variables
  double maxInf = 0.0;
  int maxInfDir = 0;
  DcoObjectInt * maxInfIntObject = 0;
  // max score variables
  double maxScore = 0.0;
  int maxScoreDir = 0;
  DcoObjectInt * maxScoreIntObject = 0;
  // get model
  DcoModel * model = dynamic_cast<DcoModel*>(model_);
  // get objective coef
  double const * objCoef = model->getObjCoef();
  // for all objects
  for (int i=0; i<model->numObjects(); ++i) {
    // TODO: currently all integer object.
    int preferDir;
    DcoObjectInt * intObject = dynamic_cast<DcoObjectInt*>(model->objects(i));
    double infeasibility = intObject->infeasibility(model, preferDir);
    if (infeasibility) {
      ++numInfs;
      if (infeasibility > maxInf) {
	maxInfIntObject = intObject;
	maxInfDir = preferDir;
	maxInf = infeasibility;
      }
      int col = intObject->columnIndex();
      double score = ALPS_FABS(objCoef[col] * infeasibility);
      if (score > maxScore) {
	maxScoreIntObject = intObject;
	maxScoreDir = preferDir;
	maxScore = score;
      }
    }
  }
  if (numInfs==0) {
    std::cerr << "All variables are feasible." << std::endl;
    throw std::exception();
  }
  if (maxScoreIntObject) {
    maxInfIntObject = maxScoreIntObject;
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
int DcoBranchStrategyMaxInf::betterBranchObject(BcpsBranchObject * thisOne,
						BcpsBranchObject * bestSoFar) {
  return thisOne->getDirection();
}

//#############################################################################
