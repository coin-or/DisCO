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

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinWarmStartBasis.hpp"

#include "Alps.h"

#include "DcoBranchStrategyRel.hpp"
#include "DcoHelp.hpp"
#include "DcoObjectInt.hpp"

//#############################################################################

struct DcoPseuoGreater
{
    bool operator()(double x, double y) const {
	return (x > y);
    }
};

//#############################################################################

// Copy constructor
DcoBranchStrategyRel::DcoBranchStrategyRel (
    const DcoBranchStrategyRel & rhs
    )
    :
    BcpsBranchStrategy()
{
    bestChangeUp_ = rhs.bestChangeUp_;
    bestNumberUp_ = rhs.bestNumberUp_;
    bestChangeDown_ = rhs.bestChangeDown_;
    bestNumberDown_ = rhs.bestNumberDown_;
    reliability_ = rhs.reliability_;
    type_ = rhs.type_;
}

//#############################################################################

/** Compare object thisOne to the bestSorFar. The compare is based on
    pseudocost. */
int
DcoBranchStrategyRel::betterBranchObject(BcpsBranchObject * thisOne,
					     BcpsBranchObject * bestSoFar)
{
    int betterDirection = 0;
    double bestChange;

    if (!bestSoFar) {
	bestChange = -1.0;
    }
    else {
	bestChange = bestChangeUp_;
    }

    double upCost = thisOne->getUpScore();

    if (upCost > bestChange) {
	betterDirection = thisOne->getDirection();
	bestChangeUp_ = upCost;
    }

    return betterDirection;
}

//#############################################################################

/** Create a set of candidate branching objects. */
int DcoBranchStrategyRel::createCandBranchObjects(int numPassesLeft,
						  double ub) {
  int bStatus = 0;
  int pass, colInd;
  int preferDir, saveLimit;
  int numFirsts  = 0;
  int numInfs = 0;
  int minCount = 0;
  int numLowerTightens = 0;
  int numUpperTightens = 0;
  double lpX, score, infeasibility, downDeg, upDeg, sumDeg = 0.0;
  bool roundAgain, downKeep, downGood, upKeep, upGood;
  int * lbInd = NULL;
  int * ubInd = NULL;
  double * newLB = NULL;
  double * newUB = NULL;
  double * origSolution = NULL;
  DcoModel * model = dynamic_cast<DcoModel *>(model_);
  int MAX_ITER = 2;
#if defined(__OA__)
  OsiSolverInterface * solver = model->solver();
#else
  OsiConicSolverInterface * solver = model->solver();
#endif
  int numCols = model->getNumCols();
  int numObjects = model->numObjects();
  //------------------------------------------------------
  // Check if max time is reached or no pass is left.
  //------------------------------------------------------
  double timeLimit = model->AlpsPar()->entry(AlpsParams::timeLimit);
  AlpsKnowledgeBroker * broker = model->getKnowledgeBroker();
  bool maxTimeReached = (broker->timer().getTime() > timeLimit);
  bool selectNow = false;
  if (maxTimeReached || !numPassesLeft) {
    selectNow = true;
#ifdef DISCO_DEBUG
    printf("REL: CREATE: maxTimeReached %d, numPassesLeft %d\n",
	   maxTimeReached, numPassesLeft);
#endif
  }
  // Store first time objects.
  std::vector<DcoObjectInt*> firstObjects;
  // Store infeasible objects.
  std::vector<DcoObjectInt*> infObjects;
  // TODO: check if sorting is expensive.
  std::multimap<double, DcoObjectInt*, DcoPseuoGreater> sortedObjects;
  double objValue = solver->getObjSense() * solver->getObjValue();
  // get current bounds
  double const * colLB = solver->getColLower();
  double const * colUB = solver->getColUpper();
  // save original bounds, solver bounds will get updated in strong branching
  double * origColLB = new double[numCols];
  double * origColUB = new double[numCols];
  std::copy(colLB, colLB+numCols, origColLB);
  std::copy(colUB, colUB+numCols, origColUB);
  // get look ahead parameter
  int lookAhead = dynamic_cast<DcoParams*>
    (model->DcoPar())->entry(DcoParams::lookAhead);
  //------------------------------------------------------
  // Backup solver status and mark hot start.
  //-----------------------------------------------------
  origSolution = new double[numCols];
  std::copy(solver->getColSolution(), solver->getColSolution()+numCols,
	    origSolution);
  //------------------------------------------------------
  // Find the infeasible objects.
  // NOTE: we might go round this loop twice if we are feed in
  //       a "feasible" solution.
  //------------------------------------------------------
  for (pass = 0; pass < 2; ++pass) {
    numInfs = 0;
    BcpsObject * object = NULL;
    infObjects.clear();
    firstObjects.clear();
    // iterate over objects and store infeasible integer objects in infObjects
    // populate firstObjects
    for (int i=0; i<numObjects; ++i) {
      object = model->objects(i);
      infeasibility = object->infeasibility(model, preferDir);
      if (infeasibility) {
	++numInfs;
	DcoObjectInt * intObject = dynamic_cast<DcoObjectInt*>(object);
	if (intObject) {
	  infObjects.push_back(intObject);
	  if (!selectNow) {
	    minCount =
	      ALPS_MIN(intObject->pseudocost().getDownCount(),
		       intObject->pseudocost().getUpCount());
	    if (minCount < 1) {
	      firstObjects.push_back(intObject);
	    }
	  }
	  intObject = NULL;
	}
	else {
	  // TODO: currently all are integer objects.
	}
      }
    }
    if (numInfs) {
#ifdef DISCO_DEBUG_MORE
      std::cout << "REL: numInfs = " << numInfs
		<< std::endl;
#endif
      break;
    }
    else if (pass==0) {
      // The first pass and is IP feasible.
#ifdef DISCO_DEBUG
      std::cout << "REL: given a feasible sol" << std::endl;
#endif
      roundAgain = false;
      CoinWarmStartBasis * ws =
	dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
      if (!ws) {
	break;
      }
      // Force solution values within bounds
      for (int i=0; i<numCols; ++i) {
	lpX = origSolution[i];
	if (lpX < colLB[i]) {
	  origSolution[i] = colLB[i];
	  roundAgain = true;
	  ws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
	}
	else if (lpX > colUB[i]) {
	  origSolution[i] = colUB[i];
	  roundAgain = true;
	  ws->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
	}
      }
      if (roundAgain) {
	// Need resolve and do the second round selection.
	solver->setWarmStart(ws);
	delete ws;
	// Resolve.
	solver->resolve();
	if (!solver->isProvenOptimal()) {
	  // Become infeasible, can do nothing.
	  bStatus = -2;
	  goto TERM_CREATE;
	}
	else {
	  // Save new lp solution.
	  std::copy(solver->getColSolution(),
		    solver->getColSolution()+numCols, origSolution);
	  objValue = solver->getObjSense() * solver->getObjValue();
	}
      }
      else {
	delete ws;
	break;
      }
    }
  } // EOF 2 pass
  //--------------------------------------------------
  // If we have a set of first time object,
  // branch up and down to initialize pseudo-cost.
  //--------------------------------------------------
  numFirsts = static_cast<int> (firstObjects.size());
  if (numFirsts>0) {
    CoinWarmStart * ws = solver->getWarmStart();
    solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
    int maxIter = ALPS_MAX(model->getAveIterations(), 50);
    solver->setIntParam(OsiMaxNumIterationHotStart, MAX_ITER);
    solver->markHotStart();
    lbInd = new int [numFirsts];
    ubInd = new int [numFirsts];
    newLB = new double [numFirsts];
    newUB = new double [numFirsts];
    for (int i=0; i<numFirsts && bStatus!=-2; ++i) {
      colInd = firstObjects[i]->columnIndex();
      lpX = origSolution[colInd];
      DcoStrongBranch(model, objValue, colInd, lpX,
		      origColLB, origColUB,
		      downKeep, downGood, downDeg,
		      upKeep, upGood, upDeg);
      if(!downKeep && !upKeep) {
	// Both branch can be fathomed
	bStatus = -2;
      }
      else if (!downKeep) {
	// Down branch can be fathomed.
	lbInd[numLowerTightens] = colInd;
	newLB[numLowerTightens++] = ceil(lpX);
	//break;
      }
      else if (!upKeep) {
	// Up branch can be fathomed.
	ubInd[numUpperTightens] = colInd;
	newUB[numUpperTightens++] = floor(lpX);
	// break;
      }
      // Update pseudocost.
      if(downGood) {
	firstObjects[i]->pseudocost().update(-1, downDeg, lpX);
	model->setSharedObjectMark(firstObjects[i]->getObjectIndex());
      }
      // todo(aykut) is it downGood or upGood?
      if(upGood) {
	firstObjects[i]->pseudocost().update(1, upDeg, lpX);
	model->setSharedObjectMark(firstObjects[i]->getObjectIndex());
      }
    }

    //--------------------------------------------------
    // Set new bounds in lp solver for resolving
    //--------------------------------------------------

    if (bStatus != -2) {
      if (numUpperTightens > 0) {
	bStatus = -1;
	for (int i=0; i<numUpperTightens; ++i) {
	  solver->setColUpper(ubInd[i], newUB[i]);
	}
      }
      if (numLowerTightens > 0) {
	bStatus = -1;
	for (int i=0; i<numLowerTightens; ++i) {
	  solver->setColLower(lbInd[i], newLB[i]);
	}
      }
    }
    //--------------------------------------------------
    // Unmark hotstart and recover LP solver.
    //--------------------------------------------------
    solver->unmarkHotStart();
    solver->setColSolution(origSolution);
    solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
    solver->setWarmStart(ws);
    delete ws;
  }

  //std::cout << "REL: bStatus = " << bStatus << std::endl;

  if (bStatus < 0) {
    // Infeasible or monotone.
    goto TERM_CREATE;
  }
  else {
    // All object's pseudocost have been initialized.
    // Sort them, and do strong branch for the unreliable one
    // NOTE: it set model->savedLpSolution.
    sumDeg = 0.0;
    for (int i=0; i<numInfs; ++i) {
      score = infObjects[i]->pseudocost().getScore();
      sumDeg += score;
      std::pair<const double, DcoObjectInt*> sa(score, infObjects[i]);
      sortedObjects.insert(sa);
#ifdef DISCO_DEBUG_MORE
      std::cout << "col[" << infObjects[i]->columnIndex() << "]="
		<< score << ", "<< std::endl;
#endif
    }
    int numNotChange = 0;
    std::multimap< double, DcoObjectInt*, DcoPseuoGreater>::iterator pos;
    CoinWarmStart * ws = solver->getWarmStart();
    solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
    //todo(aykut) parametrize the following 50
    int maxIter = ALPS_MAX(model->getAveIterations(), MAX_ITER);
    solver->setIntParam(OsiMaxNumIterationHotStart, maxIter);
    solver->markHotStart();
    DcoObjectInt *bestObject = NULL;
    double bestScore = -10.0;
    for (pos=sortedObjects.begin(); pos!=sortedObjects.end(); ++pos) {
      DcoObjectInt * intObject  = pos->second;
      colInd = intObject->columnIndex();
#ifdef DISCO_DEBUG_MORE
      std::cout << "col[" << colInd << "]: "
		<< "score=" << pos->first
		<< ", upCount=" << intObject->pseudocost().getUpCount()
		<<", downCount="<< intObject->pseudocost().getDownCount()
		<< std::endl;
#endif
      // Check if reliable.
      int objReliability=ALPS_MIN(intObject->pseudocost().getUpCount(),
				  intObject->pseudocost().getDownCount());
      //todo(aykut): Should reability_ be a parameter?
      if (objReliability < reliability_) {
	// Unrelible object. Do strong branching.
	lpX = origSolution[colInd];
	DcoStrongBranch(model, objValue, colInd, lpX,
			origColLB, origColUB,
			downKeep, downGood, downDeg,
			upKeep, upGood, upDeg);
	// Update pseudocost.
	if(downGood) {
	  intObject->pseudocost().update(-1, downDeg, lpX);
	}
	if(upGood) {
	  intObject->pseudocost().update(1, upDeg, lpX);
	}
      }
      // Compare with the best.
      if (intObject->pseudocost().getScore() > bestScore) {
	bestScore = intObject->pseudocost().getScore();
	bestObject = intObject;
	// Reset
	numNotChange = 0;
      }
      else {
	// If best doesn't change for "lookAhead" comparisons, then
	// the best is reliable.
	if (++numNotChange > lookAhead) {
	  if (bestObject->pseudocost().getUpCost() >
	      bestObject->pseudocost().getDownCost()) {
	    preferDir = 1;
	  }
	  else {
	    preferDir = -1;
	  }
	  break;
	}
      }
    }
    solver->unmarkHotStart();
    solver->setColSolution(origSolution);
    solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
    solver->setWarmStart(ws);
    delete ws;
    model->setSolEstimate(objValue + sumDeg);
    assert(bestObject != NULL);
    bestBranchObject_ = bestObject->createBranchObject(model, preferDir);
    numBranchObjects_ = 1;
    branchObjects_ = new BcpsBranchObject* [1];
    branchObjects_[0] = bestBranchObject_;
  }
 TERM_CREATE:
  //------------------------------------------------------
  // Cleanup.
  //------------------------------------------------------
  delete [] lbInd;
  delete [] ubInd;
  delete [] newLB;
  delete [] newUB;
  delete[] origColLB;
  delete[] origColUB;
  delete [] origSolution;
  return bStatus;
}

//#############################################################################
