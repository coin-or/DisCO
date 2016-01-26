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

#include <cassert>
#include <iostream>
#include <utility>
#include <cmath>
#include <vector>

#include "CoinUtility.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"

#include <CglConicGD1.hpp>
#if defined(__OA__)
#include <CglConicIPM.hpp>
#include <OsiMosekSolverInterface.hpp>
#endif
#include "AlpsKnowledge.h"
#include "AlpsEnumProcessT.h"
#include "AlpsKnowledgeBroker.h"
#include "AlpsTreeNode.h"

#include "BcpsBranchStrategy.h"

#include "Dco.hpp"
#include "DcoBranchObjectInt.hpp"
#include "DcoBranchObjectBilevel.hpp"
#include "DcoConstraint.hpp"
#include "DcoHelp.hpp"
#include "DcoTreeNode.hpp"
#include "DcoModel.hpp"
#include "DcoNodeDesc.hpp"
#include "DcoModel.hpp"
#include "DcoObjectInt.hpp"
#include "DcoParams.hpp"
#include "DcoSolution.hpp"
//#include "DcoVariable.hpp"

#define REMOVE_SLACK 0
#define DISCO_SLACK_MAX 4

//#############################################################################

AlpsTreeNode*
DcoTreeNode::createNewTreeNode(AlpsNodeDesc *&desc) const
{
  DcoModel* model = dynamic_cast<DcoModel*>(desc_->getModel());
  BcpsBranchStrategy *strategy = model->branchStrategy();

  DcoBranchingStrategy type =
    static_cast<DcoBranchingStrategy>(strategy->getType());

  switch(type){

  case DcoBranchingStrategyMaxInfeasibility:
  case DcoBranchingStrategyPseudoCost:
  case DcoBranchingStrategyReliability:
  case DcoBranchingStrategyStrong:
    {
      double estimate = solEstimate_;

      // Set solution estimate for this nodes.
      // double solEstimate = quality_ + sum_i{min{up_i, down_i}}
      int branchDir=dynamic_cast<DcoNodeDesc *>(desc)->getBranchedDir();
      int branchInd=dynamic_cast<DcoNodeDesc *>(desc)->getBranchedInd();
      double lpX = dynamic_cast<DcoNodeDesc *>(desc)->getBranchedVal();
      double f = lpX - floor(lpX);
      assert(f > 0.0);

      int objInd = model->getIntObjIndices()[branchInd];
      DcoObjectInt *obj = dynamic_cast<DcoObjectInt *>(model->objects(objInd));

      if (branchDir == -1) {
	estimate -= (1.0-f) * obj->pseudocost().getUpCost();
      }
      else {
	estimate -= f * obj->pseudocost().getDownCost();
      }

#ifdef DISCO_DEBUG_MORE
      printf("DISCO:createNewTreeNode: quality=%g, solEstimate=%g\n",
	     quality_, solEstimate_);
#endif
      break;
    }

  case DcoBranchingStrategyBilevel:

    break;

  }

  // Create a new tree node
  DcoTreeNode *node = new DcoTreeNode(desc);
  desc = NULL;

  return node;
}

//#############################################################################

// NOTE: if rampup,
// - parent must be explicit if not NULL,
// - this node is explicit.

int
DcoTreeNode::process(bool isRoot, bool rampUp)
{
  DcoReturnStatus returnStatus = DcoReturnStatusUnknown;
  DcoLpStatus lpStatus = DcoLpStatusUnknown;
  int j, k = -1;
  int numCols, numRows, numCoreCols, numCoreRows;
  int numStartRows, origNumStartRows;
  int maxNumCons, tempNumCons = 0;

  int numIntInfs = 0;
  int numObjInfs = 0;

  int voilatedNumCons = 0;
  int origNumOldCons = 0;
  int currNumOldCons = 0;
  int newNumCons = 0;
  int maxNewNumCons = 0;

  int numAppliedCons = 0;
  int cutStrategy;
  int conicCutStrategy;

  // Only autmatic stategy has depth limit.
  int maxConstraintDepth = 20;

  int numPassesLeft = 0;
  int bStatus = -1;

  double cutoff = ALPS_INC_MAX;
  double parentObjValue = getQuality();
  double preObjValue = -ALPS_OBJ_MAX;
  double improvement = 100.0;

  double *currLpSolution = NULL;

  bool keepOn = true;
  bool needBranch = false;
  bool lpFeasible = false;
  bool foundSolution = false;
  bool genConsHere = false;
  bool genConicCutsHere = false;
  bool shareCon = false;
  bool shareVar = false;

  CoinWarmStartBasis::Status rowStatus;
  DcoConstraint *aCon = NULL;
  BcpsObject **newConstraints = NULL;

  DcoSolution *ipSol = NULL;

  int numDelRows = 0;
  int *delRow = NULL;
  int *oldConsPos = NULL;

  std::vector<int> delIndices;

  DcoModel* model = dynamic_cast<DcoModel*>(desc_->getModel());

  AlpsPhase phase = knowledgeBroker_->getPhase();

  int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);
  int hubMsgLevel = model->AlpsPar()->entry(AlpsParams::hubMsgLevel);
  int workerMsgLevel = model->AlpsPar()->entry(AlpsParams::workerMsgLevel);

  DcoParams * DcoPar = model->DcoPar();

  int maxPass = DcoPar->entry(DcoParams::cutPass);
  int quickCutPass = DcoPar->entry(DcoParams::quickCutPass);

  double tailOffTol = DcoPar->entry(DcoParams::tailOff);

  if (msgLevel >= 100){
    std::cout << std::endl;
    std::cout << "**********************************" << std::endl;
    std::cout << "* NOW PROCESSING NODE " << index_ << std::endl;
    std::cout << "**********************************" << std::endl;
    printf("Initial Bound:  %.3f\n", quality_);
  }

  //if ( (quickCutPass > 0 ) && (phase == AlpsPhaseRampup) ) {
  if (quickCutPass > 0) {
    std::cout << "+++ Quick branching; pass is "<<quickCutPass<<std::endl;
    maxPass = quickCutPass + 1;
  }
  else if (maxPass < ALPS_INT_MAX) {
    // Add one to solve initial subprolem.
    ++maxPass;
  }

  shareCon = DcoPar->entry(DcoParams::shareConstraints);
  shareVar = DcoPar->entry(DcoParams::shareVariables);

  cutoff = model->getCutoff();

  numPassesLeft = model->getNumBranchResolve();

  // std::cout << "numPassesLeft = " << numPassesLeft << std::endl;

  //------------------------------------------------------
  // Check if this can be fathomed by objective cutoff.
  //------------------------------------------------------

#if 0
  std::cout << "parentObjValue = " << parentObjValue
	    << "cutoff = " << cutoff << std::endl;
#endif

  if (parentObjValue > cutoff && !isRoot) {
    setStatus(AlpsNodeStatusFathomed);
    //std::cout << "fathom!" <<std::endl;
    goto TERM_PROCESS;
  }
  else {
    model->setActiveNode(this);
    model->addNumNodes();
  }

  //------------------------------------------------------
  // Get model information and parameters.
  //------------------------------------------------------

  numCols = model->solver()->getNumCols();
  numRows = model->solver()->getNumRows();
  numCoreCols = model->getNumCoreVariables();
  numCoreRows = model->getNumCoreConstraints();
  numAppliedCons =  numRows - numCoreRows;

  maxNumCons = model->getMaxNumCons();

  currLpSolution = new double [numCols];

  //------------------------------------------------------
  // Decides if can generate constraints.
  //------------------------------------------------------

  // Mark if this node is root or not.
  model->isRoot_ = isRoot;

  genConsHere = false;
  genConicCutsHere = false;

  cutStrategy = model->getCutStrategy();
  conicCutStrategy = model->getConicCutStrategy();

  assert(cutStrategy != DcoCutStrategyNotSet);
  assert(conicCutStrategy != DcoConicCutStrategyNotSet);

  if (cutStrategy == DcoCutStrategyNone) {
    genConsHere = false;
  }
  else if (cutStrategy == DcoCutStrategyRoot) {
    // The original root only
    if (isRoot && (index_ == 0)) genConsHere = true;
  }
  else if (cutStrategy == DcoCutStrategyAuto) {
    if (depth_ < maxConstraintDepth) {
      if (!diving_ || isRoot) genConsHere = true;
    }
  }
  else if (cutStrategy == DcoCutStrategyPeriodic) {
    genConsHere = true;
  }
  else {
    genConsHere = true;
  }

  if (genConsHere && (phase == AlpsPhaseRampup)) {
    if (!(DcoPar->entry(DcoParams::cutRampUp))) {
      genConsHere = false;
    }
  }

  // will we generate conic cuts
  if (conicCutStrategy == DcoConicCutStrategyNone) {
    genConicCutsHere = false;
  }
  else if (conicCutStrategy == DcoConicCutStrategyRoot) {
    // The original root only
    if (isRoot && (index_ == 0)) genConicCutsHere = true;
  }
  else if (conicCutStrategy == DcoConicCutStrategyAuto) {
    // what does maxConstraintDepth mean for conic cut case?
    if (depth_ < maxConstraintDepth) {
      if (!diving_ || isRoot) genConicCutsHere = true;
    }
  }
  else if (conicCutStrategy == DcoConicCutStrategyPeriodic) {
    genConicCutsHere = true;
  }
  else {
    genConicCutsHere = true;
  }

  if (genConicCutsHere && (phase == AlpsPhaseRampup)) {
    if (!(DcoPar->entry(DcoParams::cutRampUp))) {
      genConicCutsHere = false;
    }
  }

#if defined(__OA__)
  // when OA is used we may generate cuts in every node.
  genConsHere = true;
#endif

  //--------------------------------------------------
  // Call HeurisBounding before solving for the first node.
  //--------------------------------------------------

  if (model->getNumNodes() == 1) {
    int heurStatus = callHeuristics(model, true); // before root
    if (heurStatus == 1) {
      cutoff = model->getCutoff();
    }
    else if (heurStatus == 2) {
      // Fathom this node
      goto TERM_PROCESS;
    }
  }

  //======================================================
  // Restore, load and solve the subproblem.
  // (1) LP infeasible
  //     a. set status to be fathom.
  // (2) LP feasible
  //     a. MILP feasible. Check whether need update incumbent.
  //     b. LP feasible but not MIP feasible. Check whether can be
  //        fathomed, if not, choose a branch variable.
  //======================================================

  //------------------------------------------------------
  // Extract info from this node and load subproblem into lp solver.
  //------------------------------------------------------
  installSubProblem(model);

  //------------------------------------------------------
  // Bounding, heuristic searching, constraint generating.
  //------------------------------------------------------

  numStartRows = model->solver()->getNumRows();
  origNumStartRows = numStartRows;

  origNumOldCons = numStartRows - numCoreRows;
  currNumOldCons = origNumOldCons;

#if 0
  std::cout << "PROCESS: genConsHere =" << genConsHere
	    << ", cut strategy =" << model->getCutStrategy()
	    << ", numCoreRows =" << numCoreRows
	    << ", numStartRows =" << numStartRows
	    << ", currNumOldCons =" << currNumOldCons
	    << ", index_ = " << index_ << std::endl;
#endif

  if (currNumOldCons > 0) {
    oldConsPos = new int [currNumOldCons];
    for (k = 0; k < currNumOldCons; ++k) {
      oldConsPos[k] = k;
    }
  }

  if (genConsHere) {
    if (maxNumCons > ALPS_INT_MAX - 100) {
      maxNewNumCons = 10000;
    }
    else {
      maxNewNumCons = maxNumCons;
    }
    newConstraints = new BcpsObject* [maxNewNumCons];
  }

  model->boundingPass_ = 0;
  while (keepOn && (model->boundingPass_ < maxPass)) {
    ++(model->boundingPass_);
    keepOn = false;

    //--------------------------------------------------
    // Bounding to get the quality of this node.
    //--------------------------------------------------

    if ( ((knowledgeBroker_->getProcType() == AlpsProcessTypeMaster) ||
	  (knowledgeBroker_->getProcType() == AlpsProcessTypeSerial)) &&
	 isRoot && (model->boundingPass_ == 1) ) {
      if (msgLevel >= 50) {
	model->blisMessageHandler()->message(DISCO_ROOT_PROCESS,
					     model->blisMessages())
	  << model->getNumRows()
	  << model->getNumCols()
	  << CoinMessageEol;
	model->solver()->messageHandler()->setLogLevel(1);
      }
      else {
	model->solver()->messageHandler()->setLogLevel(-1);
      }
      getKnowledgeBroker()->tempTimer().start();
    }

    // Get lower bound.
    lpStatus = static_cast<DcoLpStatus> (bound(model));
    if (model->boundingPass_ == 1) {
      int iter = model->solver()->getIterationCount();
      model->addNumIterations(iter);
      if (isRoot) {
	getKnowledgeBroker()->tempTimer().stop();
	if (((knowledgeBroker_->getProcType()==AlpsProcessTypeMaster)||
	     (knowledgeBroker_->getProcType()==AlpsProcessTypeSerial))
	    && (msgLevel >= 50)) {
	  model->solver()->messageHandler()->setLogLevel(0);
	  model->blisMessageHandler()->message(DISCO_ROOT_TIME,
					       model->blisMessages())
	    << getKnowledgeBroker()->tempTimer().getCpuTime()
	    << CoinMessageEol;
	}
      }
    }

    switch(lpStatus) {
    case DcoLpStatusOptimal:
      // Check if IP feasible
      ipSol = model->feasibleSolution(numIntInfs, numObjInfs);
      if (ipSol) {
	// IP feasible
	model->storeSolution(DcoSolutionTypeBounding, ipSol);
	// Update cutoff
	cutoff = model->getCutoff();
	setStatus(AlpsNodeStatusFathomed);
	goto TERM_PROCESS;
      }
      else {
	if (quality_ > cutoff) {
	  setStatus(AlpsNodeStatusFathomed);
	  goto TERM_PROCESS;
	}
	// call reduced cost fix when solver is COLA, do not call for other
	// solver, ie. mosek, cplex,
	reducedCostFix(model);
	//------------------------------------------
	// Check if tailoff or have to keep on.
	//------------------------------------------

	// at this point there is no solution both integer and conic feasible.
	// if there is a fractional variable, branch.
	if (numIntInfs) {
	  needBranch = true;
	  keepOn = false;
	}
	else {
	  // we have an integer solution that is not conic feasible
#if defined(__OA__)
#if defined(DISCO_DEBUG)
	  std::cout << "===============================================================================" << std::endl;
	  std::cout << "Node index " << index_ << std::endl;
	  std::cout << "Integer feasibility obtained." << std::endl;
	  std::cout << "Generate cuts for conic feasibility." << std::endl;
	  std::cout << "Number of infeasible conic constraints: " << numObjInfs << std::endl;
	  std::cout << "===============================================================================" << std::endl;
#endif
	  keepOn = true;
	  // this will call cut generators
	  // in ipm cut generator we will check if the solution is integer. We will call
	  // cut procedure if it is integer.
#else
	  keepOn = false;
	  needBranch = true;
#endif
	}
	if (model->boundingPass_ > 1) {
	  improvement = quality_ - preObjValue;
	  if (improvement > tailOffTol) {
	    // NOTE: still need remove slacks, although
	    //       tailoff.
	    // we do not want this. Conic cuts will improve the bound
	    // but we want to branch if there is at least one fractional variable.
	    //keepOn = true;
	  }

#if 0
	  std::cout << "PROCESS: boundingPass_["
		    << model->boundingPass_ << "], improvement="
		    << improvement << ", tailOffTol=" << tailOffTol
		    << ", preObj=" << preObjValue
		    << ", newObj=" << quality_
		    << std::endl;
#endif
	}
	else {
	  keepOn = true;
	}

	if (model->boundingPass_ == maxPass) {
	  if (!numIntInfs && !numObjInfs) {
	    /* No fractional objects, have to process one more pass.
	       For applications like VRP. */
	    ++maxPass;
	    keepOn = true;
	  }
	}

	//------------------------------------------
	// Update previous objective value.
	//------------------------------------------

	preObjValue = quality_;

	//------------------------------------------
	// Remove non-core slack constraints.
	//------------------------------------------

	numRows = model->getNumRows();

	if ( genConsHere &&
	     //(improvement > tailOffTol) &&
	     //(numRows > numStartRows) ) {
	     (numRows > numCoreRows) ) {

#if 1
	  if ( (numStartRows + newNumCons != numRows) ||
	       (numCoreRows+currNumOldCons +newNumCons != numRows) ){

	    // // todo(aykut) this is done blindly and will create future problems
	    // std::cout << "ERROR: numRows=" << numRows
	    //           << "; numCoreRows=" << numCoreRows
	    //           << "; numStartRows=" << numStartRows
	    //           << "; newNumCons=" << newNumCons
	    //           << "; currNumOldCons=" << currNumOldCons
	    //           << std::endl;

	    // assert(numRows - numStartRows == newNumCons);
	  }
#endif

	  int *oldDelMark = NULL;
	  if (currNumOldCons > 0) {
	    oldDelMark = new int [currNumOldCons];
	    CoinZeroN(oldDelMark, currNumOldCons);
	  }
	  int *newDelMark = NULL;
	  if (newNumCons > 0) {
	    newDelMark = new int [newNumCons];
	    CoinZeroN(newDelMark, newNumCons);
	  }

	  const CoinWarmStartBasis* ws=
	    dynamic_cast<CoinWarmStartBasis*>
	    (model->solver()->getWarmStart());

	  // Make sure delIndices is empty.
	  assert(delIndices.size()==0);

#if REMOVE_SLACK
	  for (k = numCoreRows; k < numRows; ++k) {
	    rowStatus = ws->getArtifStatus(k);

	    if (rowStatus == CoinWarmStartBasis::basic) {
	      int count;
	      if (k < numStartRows) {
		DcoConstraint *tmpCon =
		  model->oldConstraints()[(k-numCoreRows)];
		count = tmpCon->getNumInactive() + 1;
		tmpCon->setNumInactive(count);
		if (tmpCon->getNumInactive() > DISCO_SLACK_MAX){
		  oldDelMark[(k-numCoreRows)] = 1;
		  delIndices.push_back(k);
		}
	      }
	      else {
		BcpsObject* tmpCon =
		  newConstraints[(k-numStartRows)];
		count = tmpCon->getNumInactive() + 1;
		tmpCon->setNumInactive(count);
		if (tmpCon->getNumInactive() > DISCO_SLACK_MAX){
		  newDelMark[(k-numStartRows)] = 1;
		  delIndices.push_back(k);
		}
	      }
	    }
	  }
#endif
	  numDelRows = static_cast<int> (delIndices.size());

	  if (msgLevel > 100) {
	    std::cout << "PROCESS: new cuts=" << newNumCons
		      << ", delete slack cuts=" << numDelRows
		      << ", numRows=" << numRows
		      << ", numStartRows=" <<numStartRows
		      << std::endl;
	  }


	  if (numDelRows > 0) {
	    delRow = new int [numDelRows];
	    for (k = 0; k < numDelRows; ++k) {
	      delRow[k] = delIndices[k];
#ifdef DISCO_DEBUG
	      std::cout << "REMOVE: slack row " << delRow[k]
			<< std::endl;
#endif
	    }

	    //----------------------------------
	    // Delete from lp solver.
	    //----------------------------------

	    model->solver()->deleteRows(numDelRows, delRow);

	    delete [] delRow;
	    delRow = NULL;
	    delIndices.clear();

	    //----------------------------------
	    // Delete from the old cut position array.
	    //----------------------------------

	    tempNumCons = 0;
	    for (k = 0; k < currNumOldCons; ++k) {
	      if (oldDelMark[k] != 1) {
		// Survived
		oldConsPos[tempNumCons++] = oldConsPos[k];
	      }
	    }
	    currNumOldCons = tempNumCons;
	    numStartRows = numCoreRows + currNumOldCons;

	    //----------------------------------
	    // Delete from new cut vector.
	    //----------------------------------

	    //std::cout << std::endl;
	    tempNumCons = 0;
	    for (k = 0; k < newNumCons; ++k) {
	      if (newDelMark[k] == 1) {
		// Deleted
#ifdef DISCO_DEBUG_MORE
		std::cout << "delete cut " << k
			  << ", size="
			  << dynamic_cast<DcoConstraint*>(newConstraints[k])->getSize()
			  << std::endl;
#endif

		delete newConstraints[k];
		newConstraints[k] = NULL;
	      }
	      else {
		// Survived
		newConstraints[tempNumCons++] = newConstraints[k];
	      }
	    }
	    //assert(tempNumCons + numDelRows == newNumCons);
	    numAppliedCons -= newNumCons;
	    numAppliedCons += tempNumCons;
	    newNumCons = tempNumCons;

	    //----------------------------------
	    // Resolve to update solution info in lp solver.
	    //----------------------------------

	    int easy = 2;
	    model->solver()->setHintParam(OsiDoInBranchAndCut,
					  true, OsiHintDo, &easy);
	    model->solver()->resolve();
	    model->solver()->setHintParam(OsiDoInBranchAndCut,
					  true, OsiHintDo, NULL) ;



#ifdef DISCO_DEBUG
	    if (model->solver()->getIterationCount() != 0) {
	      // TODO: maybe some cuts become slack again
#ifdef DISCO_DEBUG
	      std::cout << "SLACK: resolve changed solution!"
			<< ", iter="
			<< model->solver()->getIterationCount()
			<< std::endl;
#endif
	    }
	    else {
#ifdef DISCO_DEBUG
	      std::cout<<"SLACK: resolve don't changed solution!"
		       << std::endl;
#endif
	    }

	    double tempOV = model->solver()->getObjValue();
	    double ovDiff = fabs(quality_ - tempOV);

	    if (ovDiff /(1.0 + tempOV) > 1.0e-3) {
	      std::cout << "ERROR: SLACK: quality_("<<quality_
			<< ") != tempOV(" << tempOV
			<< ")" << std::endl;
	      assert(0);
	    }
	    else {
	      std::cout << "AFTER SLACK: quality_("<<quality_
			<< ") == tempOV(" << tempOV
			<< ")" << std::endl;
	    }
#endif
	  }
	  delete ws;
	  delete [] newDelMark;
	  delete [] oldDelMark;
	}
      }

      break;
    case DcoLpStatusAbandoned:
      assert(0);
      returnStatus = DcoReturnStatusErrLp;
      goto TERM_PROCESS;
    case DcoLpStatusDualInfeasible:
      //todo(aykut) what happens if the relaxation is unbounded?
      // this is the current bug


      // FIXME: maybe also primal infeasible
#ifdef DISCO_DEBUG
      assert(0);
#endif
      returnStatus = DcoReturnStatusUnbounded;
      goto TERM_PROCESS;
    case DcoLpStatusPrimalInfeasible:
      setStatus(AlpsNodeStatusFathomed);
      quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
      returnStatus = DcoReturnStatusInfeasible;
      goto TERM_PROCESS;
    case DcoLpStatusDualObjLim:
      setStatus(AlpsNodeStatusFathomed);
      quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
      returnStatus = DcoReturnStatusOverObjLim;
      goto TERM_PROCESS;
    case DcoLpStatusPrimalObjLim:
    case DcoLpStatusIterLim:
      /* Can't say much, need branch */
      needBranch = true;
#ifdef DISCO_DEBUG
      assert(0);
#endif
      returnStatus = DcoReturnStatusBranch;
      goto TERM_BRANCH;
      break;
    default:
#ifdef DISCO_DEBUG
      std::cout << "PROCESS: unknown status "  <<  lpStatus << std::endl;
      assert(0);
#endif
      break;
    }

// run generation/heuristics when solver is cola or OA.
#if defined(__OA__) || defined(__COLA__)
    //--------------------------------------------------
    // Call heuristics.
    //--------------------------------------------------

#if 0
    std::cout << "Process " << getKnowledgeBroker()->getProcRank()
	      << " : model->boundingPass_=" << model->boundingPass_
	      << " ; maxPass = " << maxPass
	      << " ; keepOn = " << keepOn << std::endl;
#endif

    // todo(aykut) disable heuristics
    // if (keepOn) {
    //   int heurStatus = callHeuristics(model, false);
    //   if (heurStatus == 1) {
    //	cutoff = model->getCutoff();
    //   }
    //   else if (heurStatus == 2) {
    //	// Fathom this node
    //	goto TERM_PROCESS;
    //   }
    // }

    //--------------------------------------------------
    // Generate constraints.
    //--------------------------------------------------

#ifdef DISCO_DEBUG
    std::cout << "keepOn = " << keepOn
	      << ", geneConsHere = " << genConsHere
	      << ", numAppliedCons = " << numAppliedCons
	      << ", maxNumCons = " << maxNumCons
	      << ", boundingPass = " << model->boundingPass_
	      << ", maxPass =" << maxPass
	      << std::endl;
#endif
    if ( keepOn && genConsHere &&
	 (numAppliedCons < maxNumCons) &&
	 (model->boundingPass_ < maxPass) ) {
      //OsiCuts newOsiCuts;
      BcpsConstraintPool newConPool;

      memcpy(currLpSolution,
	     model->getLpSolution(),
	     numCols * sizeof(double));

      // Get violated constraints that are from other processes.
      tempNumCons = newConPool.getNumConstraints();
      getViolatedConstraints(model, currLpSolution,
			     *(model->constraintPoolReceive()));
      // todo(aykut) fix typo voilated
      voilatedNumCons = newConPool.getNumConstraints() - tempNumCons;

      // Generate constraints (only if no violated).
      if (voilatedNumCons == 0) {
	lpStatus = static_cast<DcoLpStatus>
	  (generateConstraints(model, newConPool));

	if (lpStatus != DcoLpStatusOptimal) {
	  setStatus(AlpsNodeStatusFathomed);
	  quality_ = -ALPS_OBJ_MAX; // Remove it as soon as possilbe
	  goto TERM_PROCESS;
	}
      }

      tempNumCons = newConPool.getNumConstraints();
#ifdef DISCO_DEBUG
      std::cout << "Generated " << tempNumCons << " many cuts." << std::endl;
#endif

      if (tempNumCons > 0) {
	// Select and install new constraints
	applyConstraints(model, currLpSolution, newConPool);
	// Some weak/parallel/dense constraints might be discarded.
	tempNumCons = newConPool.getNumConstraints();
	if (tempNumCons > 0) {
	  keepOn = true;
	}
	else {
	  keepOn = false;
	}

	// Move cuts from pool to array newConstraints.
	for (k = 0; k < tempNumCons; ++k) {
	  aCon = dynamic_cast<DcoConstraint *>
	    (newConPool.getConstraint(k));
	  newConstraints[newNumCons++] = aCon;
	  if (newNumCons >= maxNewNumCons) {
	    // No space, need resize
#ifdef DISCO_DEBUG
	    std::cout << "NEWCUT: resize, maxNewNumCons = "
		      << maxNewNumCons << std::endl;
#endif
	    maxNewNumCons *= 2;
	    BcpsObject **tempNews = new BcpsObject* [maxNewNumCons];
	    memcpy(tempNews,
		   newConstraints,
		   newNumCons * sizeof(BcpsObject *));
	    delete [] newConstraints;
	    newConstraints = tempNews;
	  }

	  // Make a copy to send pool if share
	  if (shareCon && (voilatedNumCons == 0)) {
	    if (aCon->getValidRegion() == BcpsValidGlobal) {
	      model->constraintPoolSend()->
		addConstraint(new DcoConstraint(*aCon));
	    }
#if 0
	    std::cout << "+++ Num of send new constraint = "
		      << model->constraintPoolSend()->getNumConstraints()
		      << std::endl;
#endif
	  }
	}

	newConPool.clear();
	numAppliedCons += tempNumCons;
      }
      else { // Didn't generate any new constraints.
	keepOn = false;
#if defined(__OA__)
	// if OA algorithm is used and we did not generated any cuts
	// check if the current sol is integer
	// if it is integer, this means ipm method did not result any cuts.
	// we can use CglConicOA to generate cuts.
#ifdef DISCO_DEBUG_MORE
	std::cout << "Solution is integer and conic infeasible. " << std::endl;
	std::cout << "Failed to generate cuts." << std::endl;
#endif
#endif
      }
    }
    else { // Don't allow to generate constraints.
      keepOn = false;
    }
#endif // end of __COLA__ defined block
  } // EOF bounding/cutting/heuristics loop


  // todo(aykut) This is the proper place to generate and add conic cuts

    //------------------------------------------------------
    // Select branching object
    //------------------------------------------------------

  //FIXME: Turn this into a function
 TERM_BRANCH:

#ifdef DISCO_DEBUG_MORE
  printf("needBranch = %d\n", needBranch);
#endif

  //FIXME: Isn't needBranch always true here?
  if (needBranch) {

    bStatus = -1;

    while (bStatus == -1) {
      foundSolution = false;
      if(getKnowledgeBroker()->getProcRank() == -1) {
	std::cout << "*** I AM RANK ONE: before choose:bStatus = "
		  << bStatus << std::endl;
      }
      // todo(aykut) what if problem becomes integer feasible after the
      // branch on the prev iteration
      bStatus = selectBranchObject(model,
				   foundSolution,
				   numPassesLeft);
      --numPassesLeft;

      if (bStatus == -1) {
	// bStatus is -1 if a varaible is fixed to some value.
	lpFeasible = model->resolve();
	//resolved = true ;
#ifdef DISCO_DEBUG_MORE
	printf("Resolve since some col fixed, Obj value %g, numRows %d, cutoff %g\n",
	       model->solver()->getObjValue(),
	       model->solver()->getNumRows(),
	       cutoff);
#endif
	if (lpFeasible) {
	  // Update new quality.
	  quality_ = model->solver()->getObjValue();
	  if (quality_ > cutoff) {
	    bStatus = -2;
	  }
	  // Check if feasible at the other branch due to random LP
	  ipSol = model->feasibleSolution(numIntInfs, numObjInfs);
	  if (ipSol) {
	    // IP feasible
	    model->storeSolution(DcoSolutionTypeBounding, ipSol);
	    // Update cutoff
	    cutoff = model->getCutoff();
	    setStatus(AlpsNodeStatusFathomed);
	    goto TERM_PROCESS;
	  }
	  else if (!numIntInfs && numObjInfs) {
	    // problem is integer feasible but not conic feasible
	    // this can happen only if OA is used.
#ifdef DISCO_DEBUG_MORE
	    std::cout << "Some col is fixed. Problem is integer feasible but not conic feasible." << std::endl;
	    std::cout << "Fix all integer variables, solve problem with IPM..." << std::endl;
#endif
	    // fix all integer variables, load problem to ipm solver,
	    // if feasible check solution quality, if better store it, else
	    // fathom.
#if defined(__OA__)
	    OsiConicSolverInterface * ipm_solver = new OsiMosekSolverInterface();
	    integerFix(model, ipm_solver);
	    if (ipm_solver->isProvenOptimal()) {
	      quality_ = ipm_solver->getObjValue();
	      if (quality_ > cutoff) {
		bStatus = -2;
	      }
	      else {
		// Check if feasible at the other branch due to random LP
		ipSol = model->feasibleSolution(numIntInfs, numObjInfs);
		// IP feasible
		model->storeSolution(DcoSolutionTypeBounding, ipSol);
		// Update cutoff
		cutoff = model->getCutoff();
		setStatus(AlpsNodeStatusFathomed);
		goto TERM_PROCESS;
	      }
	    }
	    else if (ipm_solver->isProvenPrimalInfeasible() || ipm_solver->isProvenDualInfeasible()) {
	      // problem is infeasible when integer variables are fixed.
	      setStatus(AlpsNodeStatusFathomed);
	      goto TERM_PROCESS;
	    }
	    else {
	      std::cout << "IPM failed!" << std::endl;
	      throw std::exception();
	    }
	    delete ipm_solver;
#endif
	  }
	  if (msgLevel >= 100){
	    printf("Final Bound:  %.3f\n", quality_);
	  }
	}
	else {
	  // Should not happen. No, it will happen when other
	  // branch is ip feasible, and cause this branch to fathom
	  // when resolving. Test enigma.
	  // assert(0);
	  bStatus = -2;
	  setStatus(AlpsNodeStatusFathomed);
	}
      }

      if(getKnowledgeBroker()->getProcRank() == -1) {
	std::cout << "*** I AM RANK ONE: bStatus = " << bStatus
		  << std::endl;
      }
    }

    assert(bStatus != -1);

    //----------------------------------------------------
    // If found a branching object:
    // 1. Record basis
    // 2. Record soft var bound difference.
    // 3. Record add/del constraints.
    // NOTE: Hard var bound differences have been recorded when branch().
    //       startXXXXX have branching bounds for this node
    //----------------------------------------------------

    if (bStatus >= 0) {

#ifdef DISCO_DEBUG_MORE
      DcoBranchObjectInt *branchObject =
	dynamic_cast<DcoBranchObjectInt *>(branchObject_);
      std::cout << "SetPregnant: branchedOn = "
		<< branchObject->getObjectIndex()
		<< std::endl;
#endif
      //--------------------------------------------------
      // Mark as pregnant.
      //--------------------------------------------------

      setStatus(AlpsNodeStatusPregnant);

      //--------------------------------------------------
      // Save basis.
      //--------------------------------------------------

      CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
	(model->solver()->getWarmStart());
      DcoNodeDesc *desc = dynamic_cast<DcoNodeDesc *>(desc_);
      desc->setBasis(ws);

      //----------------------------------------------
      // Save variable/constraint bound, non-core constraints
      // and non-core variable.
      // The sizes of hard variable bound vectors are numCols.
      // The sizes of soft variable bound vectors are number
      // of modified.
      //----------------------------------------------

      int *tempVarLBPos = model->tempVarLBPos();
      int *tempVarUBPos = model->tempVarUBPos();
      //int *tempConLBPos = model->tempConLBPos();
      //int *tempConUBPos = model->tempConUBPos();

      int numModSoftColLB = 0;
      int numModSoftColUB = 0;
      const double *currColLB = model->solver()->getColLower();
      const double *currColUB = model->solver()->getColUpper();
      //const double *currRowLB = model->solver()->getRowLower();
      //const double *currRowUB = model->solver()->getRowUpper();

      double *startColLB = model->startVarLB();
      double *startColUB = model->startVarUB();
      //double *startRowLB = model->startConLB();
      //double *startRowUB = model->startConUB();

      //DcoConstraint **tempCons = NULL;


#ifdef DISCO_DEBUG_MORE
      // Debug survived old constraints.
      for (k = 0; k < currNumOldCons; ++k) {
	int oldPos = oldConsPos[k];
	DcoConstraint *aCon = model->oldConstraints()[oldPos];
	assert(aCon);
	std::cout << "SAVE: DBG: oldPos=" << oldPos
		  << ", k=" << k << ", len=" << aCon->getSize()
		  << ", node=" << index_ << std::endl;
      }
#endif

      //----------------------------------------------
      // Decide if save explicit decription.
      //----------------------------------------------

      DcoParams * DcoPar = model->DcoPar();
      int difference = DcoPar->entry(DcoParams::difference);

      if (difference == -1) {
	if (depth_ % 30 == 0 || isRoot || (phase == AlpsPhaseRampup)) {
	  explicit_ = 1;
	  //std::cout << "SAVE: node "<< index_ <<" explicitly, "
	  //  << "depth=" << depth_ << std::endl;
	}
	else {
	  explicit_ = 0;
	  //std::cout << "SAVE: node "<< index_ <<" relatively, "
	  //  << "depth=" << depth_ << std::endl;
	}
      }
      else if (difference == 0) {
	explicit_ = 1;
      }
      else {
	explicit_ = 0;
      }

      //explicit_ = 1;

      if (explicit_ || (phase == AlpsPhaseRampup) ) {
	// NOTE: full hard bound has been stored.

	int index;
	int numModify = 0;
	int numSoftVarLowers = 0;
	int numSoftVarUppers = 0;
	double value;

	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];

	double *fVarSoftLB = new double [numCols];
	double *fVarSoftUB = new double [numCols];
	int *fVarSoftLBInd = new int [numCols];
	int *fVarSoftUBInd = new int [numCols];

	for (k = 0; k < numCols; ++k) {
	  fVarSoftLB[k] = ALPS_DBL_MAX;
	  fVarSoftUB[k] = -ALPS_DBL_MAX;
	  fVarHardLB[k] = ALPS_DBL_MAX;
	  fVarHardUB[k] = -ALPS_DBL_MAX;
	  fVarHardLBInd[k] = k;
	  fVarHardUBInd[k] = k;
	}

	//------------------------------------------
	// Build the path to an explicit node.
	//------------------------------------------

	AlpsTreeNode *parent = parent_;

	// First push this node since it has branching hard bounds.
	std::vector<AlpsTreeNode*> leafToRootPath;
	leafToRootPath.push_back(this);
	DcoNodeDesc* pathDesc = NULL;

	if (phase != AlpsPhaseRampup) {
	  while(parent) {
	    leafToRootPath.push_back(parent);
	    if (parent->getExplicit()) {
	      // Reach an explicit node, then stop.
	      break;
	    }
	    else {
	      parent = parent->getParent();
	    }
	  }
	}

#ifdef DISCO_DEBUG_MORE
	std::cout << "SAVE: EXP: path len = "<<leafToRootPath.size()
		  << std::endl;
#endif
	//------------------------------------------
	// Summarize bounds.
	//------------------------------------------

	for(j = static_cast<int> (leafToRootPath.size() - 1); j > -1; --j) {

	  pathDesc = dynamic_cast<DcoNodeDesc*>
	    ((leafToRootPath.at(j))->getDesc());

	  //--------------------------------------
	  // Full variable hard bounds.
	  //--------------------------------------

	  numModify = pathDesc->getVars()->lbHard.numModify;
	  for (k = 0; k < numModify; ++k) {
	    index = pathDesc->getVars()->lbHard.posModify[k];
	    value = pathDesc->getVars()->lbHard.entries[k];
	    fVarHardLB[index] = value;
	  }

	  numModify = pathDesc->getVars()->ubHard.numModify;
	  for (k = 0; k < numModify; ++k) {
	    index = pathDesc->getVars()->ubHard.posModify[k];
	    value = pathDesc->getVars()->ubHard.entries[k];
	    fVarHardUB[index] = value;
	  }

	  //--------------------------------------
	  // Full variable soft bounds.
	  //--------------------------------------

	  numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef DISCO_DEBUG_MORE
	  std::cout << "SAVE: EXP: j=" << j << ", numModify soft lb="
		    << numModify << std::endl;
#endif
	  for (k = 0; k < numModify; ++k) {
	    index = pathDesc->getVars()->lbSoft.posModify[k];
	    value = pathDesc->getVars()->lbSoft.entries[k];
	    fVarSoftLB[index] = value;
	  }

	  numModify = pathDesc->getVars()->ubSoft.numModify;
#ifdef DISCO_DEBUG_MORE
	  std::cout << "SAVE: EXP: j=" << j << ", numModify soft ub="
		    << numModify << std::endl;
#endif
	  for (k = 0; k < numModify; ++k) {
	    index = pathDesc->getVars()->ubSoft.posModify[k];
	    value = pathDesc->getVars()->ubSoft.entries[k];
	    fVarSoftUB[index] = value;
	  }

	} // EOF of for(path)

	//------------------------------------------
	// Collect modified soft bounds at this node.
	// NOTE: Do this after collecting previous soft bounds.
	//------------------------------------------

	numModSoftColLB = 0;
	numModSoftColUB = 0;
	for (k = 0; k < numCoreCols; ++k) {
	  if (currColLB[k] != startColLB[k]) {
	    fVarSoftLB[k] = currColLB[k];
	    ++numModSoftColLB;
#ifdef DISCO_DEBUG_MORE
	    printf("Col %d, soft lb change, start %g, curr %g\n",
		   k, startColLB[k], currColLB[k]);
#endif

	  }
	  if (currColUB[k] != startColUB[k]) {
	    fVarSoftUB[k] = currColUB[k];
	    ++numModSoftColUB;
	  }
	}

#ifdef DISCO_DEBUG_MORE
	std::cout << "SAVE: EXP: THIS: numModSoftColLB = "<<numModSoftColLB
		  << ", numModSoftColUB = " << numModSoftColUB << std::endl;
#endif

	//--------------------------------------
	// Debug if bounds are consistant.
	//--------------------------------------

#ifdef DISCO_DEBUG
	for (k = 0; k < numCols; ++k) {

	  //std::cout << "EXP: COL[" << k <<"]: "
	  //      <<"hardLB=" << fVarHardLB[k]
	  //      <<", hardUB=" << fVarHardUB[k] << std::endl;

	  // Hard lower bound should not greater than
	  // hard upper bound.
	  if (fVarHardLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
	    printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
		   k, fVarHardLB[k], fVarHardUB[k]);
	    assert(0);
	  }

	  if (fVarSoftLB[k] < ALPS_BND_MAX) {
	    // Soft lower changed, and should not greater
	    // than its hard upper bound.
	    if (fVarSoftLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
	      printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
		     k, fVarSoftLB[k], fVarHardUB[k]);
	      assert(0);
	    }
	  }

	  if (fVarSoftUB[k] > -ALPS_BND_MAX) {
	    // Soft upper changed, and should not less
	    // than its hard lower bound.
	    if (fVarSoftUB[k] < fVarHardLB[k] - ALPS_GEN_TOL) {
	      printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
		     k, fVarHardLB[k], fVarSoftUB[k]);
	      assert(0);
	    }
	  }
	}
#endif

	//------------------------------------------
	// Record hard variable bounds. FULL set.
	//------------------------------------------

	desc->assignVarHardBound(numCols,
				 fVarHardLBInd,
				 fVarHardLB,
				 numCols,
				 fVarHardUBInd,
				 fVarHardUB);

	//------------------------------------------
	// Recode soft variable bound. Modified.
	//------------------------------------------

	for (k = 0; k < numCols; ++k) {
	  if (fVarSoftLB[k] < ALPS_BND_MAX) {
	    fVarSoftLBInd[numSoftVarLowers] = k;
	    fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
	  }
	  if (fVarSoftUB[k] > -ALPS_BND_MAX) {
	    fVarSoftUBInd[numSoftVarUppers] = k;
	    fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
	  }
	}


#ifdef DISCO_DEBUG_MORE
	// Print soft bounds.
	std::cout << "SAVE: EXP: numSoftVarLowers=" << numSoftVarLowers
		  << ", numSoftVarUppers=" << numSoftVarUppers
		  << std::endl;
	for (k = 0; k < numSoftVarLowers; ++k) {
	  std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
		    << fVarSoftLB[k] << std::endl;
	}
	std::cout << "------------------" << std::endl;
	for (k = 0; k < numSoftVarUppers; ++k) {
	  std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
		    << fVarSoftUB[k] << std::endl;
	}
	std::cout << "------------------" << std::endl << std::endl;
#endif

	//if ( (numSoftVarUppers > 0) || (numSoftVarLowers > 0) ) {

	// Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
	desc->assignVarSoftBound(numSoftVarLowers,
				 fVarSoftLBInd,
				 fVarSoftLB,
				 numSoftVarUppers,
				 fVarSoftUBInd,
				 fVarSoftUB);

	//------------------------------------------
	// Full set of active non-core constraints.
	//------------------------------------------
	// old constraint: model->oldConstraints_, currNumOldCons.
	// new constraint: newConstraints, newNumCons.

	BcpsObject **toAddCons = new BcpsObject * [currNumOldCons +
						   newNumCons];
	if (currNumOldCons > 0) {
	  // Hard copy of the survived old constraints.
	  for (k = 0; k < currNumOldCons; ++k) {
	    int oldPos = oldConsPos[k];
	    DcoConstraint *aCon = model->oldConstraints()[oldPos];
	    assert(aCon);
#ifdef DISCO_DEBUG
	    std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons
		      << ", k=" << k << ", len=" << aCon->getSize()
		      << ", node=" << index_ << std::endl;
#endif

	    DcoConstraint *newCon = new DcoConstraint(*aCon);
	    toAddCons[k] = newCon;
	  }
	}

	if (newNumCons > 0) {
	  // Include new constraints.
	  memcpy(toAddCons + currNumOldCons,
		 newConstraints,
		 newNumCons * sizeof(BcpsObject *));
	}

	//------------------------------------------
	// Save in description. Add first delete exiting, then add.
	// Pointers in model->oldConstraints_ can be dangling.
	// It is not safe to use model->oldConstraints_ after adding.

	// If this node is the root of a subtree, and before processing
	// it has a list of cuts, then model->oldConstraints_
	// stores pointer to the cuts when installing.

	// Need update model->constraints_ here OR do not be smart!
	// 1/6/06: Choose to be dumn.
	//------------------------------------------

	//------------------------------------------
	// Generating constraints,
	// also means that slack ones might been removed.
	//------------------------------------------

	int numTotal = currNumOldCons + newNumCons;
	desc->setAddedConstraints(numTotal, toAddCons);

#ifdef DISCO_DEBUG
	std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons
		  << ", newNumCons=" << newNumCons
		  << std::endl;
#endif

	//------------------------------------------
	// Full set of active non-core variables.
	//------------------------------------------
	// Todo.

	//------------------------------------------
	// Clear path vector.
	//------------------------------------------

	leafToRootPath.clear();
	assert(leafToRootPath.size() == 0);
      }
      else { // Relative.
	//------------------------------------------
	// Record soft bound changes for core vars.
	//------------------------------------------

	// Variable bound change
	numModSoftColLB = 0;
	numModSoftColUB = 0;
	for (k = 0; k < numCoreCols; ++k) {
	  if (currColLB[k] != startColLB[k]) {
	    tempVarLBPos[numModSoftColLB] = k;
	    /* startColLB as a temporary storage vector */
	    startColLB[numModSoftColLB] = currColLB[k];

#ifdef DISCO_DEBUG_MORE
	    printf("Col %d, soft lb change, start %g, curr %g\n",
		   k, startColLB[k], currColLB[k]);
#endif

	    ++numModSoftColLB;
	  }
	  if (currColUB[k] != startColUB[k]) {
	    tempVarUBPos[numModSoftColUB] = k;
	    startColUB[numModSoftColUB] = currColUB[k];
	    ++numModSoftColUB;
	  }
	}

#ifdef DISCO_DEBUG_MORE
	std::cout << "SAVE: REL: numModSoftColLB = "
		  << numModSoftColLB
		  << ", numModSoftColUB = "
		  << numModSoftColUB
		  << std::endl;
#endif

	if (numModSoftColLB > 0 || numModSoftColUB > 0) {
#ifdef DISCO_DEBUG
	  //assert(0);
#endif
	  desc->setVarSoftBound(numModSoftColLB,
				tempVarLBPos,
				startColLB,
				numModSoftColUB,
				tempVarUBPos,
				startColUB);
	}

	//------------------------------------------
	// TODO: Constraint bounds change.
	//------------------------------------------

#if 0
	for (k = 0; k < numCoreRows; ++k) {
	  if (currRowLB[k] != startRowLB[k]) {
	    tempConLBPos[numModSoftRowLB] = k;
	    startRowLB[numModSoftRowLB] = currRowLB[k];
	    ++numModSoftRowLB;
	  }
	  if (currRowUB[k] != startRowUB[k]) {
	    tempConUBPos[numModSoftRowUB] = k;
	    startRowUB[numModSoftRowUB] = currRowUB[k];
	    ++numModSoftRowUB;
	  }
	}
	if (numModSoftRowLB > 0 || numModSoftRowUB > 0) {
	  desc->setConSoftBound(numModSoftRowLB,
				tempConLBPos,
				startRowLB,
				numModSoftRowUB,
				tempConUBPos,
				startRowUB);
	}
#endif

	if (genConsHere) {
	  // NOTE: explicit_ can NOT do this if, since genConsHere maybe
	  //       false here, but there are maybe cons from parents.

	  //--------------------------------------
	  // Record add constraints.
	  //--------------------------------------

	  // Apend will copy old, then add new.
	  // If this node has a list of cuts before pointers in
	  // model->oldConstraints() will be kept. Safe!
	  if (newNumCons > 0) {
	    desc->appendAddedConstraints(newNumCons,
					 newConstraints);
	  }

	  //--------------------------------------
	  // Record deleted constraint positions.
	  //--------------------------------------

	  int *oldLeft = new int [origNumOldCons];
	  int leftCon;
	  CoinZeroN(oldLeft, origNumOldCons);

	  for (k = 0; k < currNumOldCons; ++k) {
	    leftCon = oldConsPos[k];
	    assert(leftCon >= 0 && leftCon < origNumOldCons);
	    oldLeft[leftCon] = 1;
	  }
	  //
	  leftCon = 0;
	  for (k = 0; k < origNumOldCons; ++k) {
	    if (oldLeft[k] == 0) {
	      // Deleted. Now oldLeft stores delete position.
	      oldLeft[leftCon++] = k;
	    }
	    //FIXME: clean
	    //assert(k < 15196924);
	  }
	  desc->delConstraints(leftCon, oldLeft);

#ifdef DISCO_DEBUG
	  std::cout << "PROCESS: ADD: new cuts=" << newNumCons
		    << ", numRows=" << model->solver()->getNumRows()
		    << ", numStartRows="<< numStartRows
		    << ", origNumStartRows="<< origNumStartRows
		    << ", num removed=" << leftCon << std::endl;
#endif

	}// EOF of if(genConsHere)
      } // EOF of relative
    }
    else if (bStatus == -2) {

#if 0
      std::cout << "bStatus = -2, fathom this node!" << std::endl;
#endif
      //branchObject->getDown()[0], branchObject->getDown()[1]);

      setStatus(AlpsNodeStatusFathomed);
    }
    else {
      throw CoinError("No branch object found", "process",
		      "DcoTreeNode");
    }
  }

  //------------------------------------------------------
  // End of process()
  //------------------------------------------------------

 TERM_PROCESS:

  bool printCutStat = false;
  if (genConsHere) {
    if ( ((getKnowledgeBroker()->getProcType()==AlpsProcessTypeMaster)||
	  (getKnowledgeBroker()->getProcType()==AlpsProcessTypeSerial)) &&
	 (msgLevel >= 10) ) {
      printCutStat = true;
      //std::cout << "+++++master print cut stats"<< std::endl;
    }
    else if ( (getKnowledgeBroker()->getProcType()==AlpsProcessTypeHub) &&
	      (hubMsgLevel > 0) ) {
      printCutStat = true;
    }
    else if ((getKnowledgeBroker()->getProcType()==AlpsProcessTypeWorker)&&
	     (workerMsgLevel > 0)) {
      printCutStat = true;
    }
  }

  if (printCutStat == true) {
    printCutStat = false;
    if ( (msgLevel >= 200) || (index_ == 0 && msgLevel >= 50) ) {
      printCutStat = true;
    }
  }

  if (printCutStat) {
    int numT = model->numCutGenerators();
    for (k = 0; k < numT; ++k) {
      if ( model->cutGenerators(k)->calls() > -1) {
	model->blisMessageHandler()->message(DISCO_CUT_STAT_NODE,
					     model->blisMessages())
	  << index_
	  << model->cutGenerators(k)->name()
	  << model->cutGenerators(k)->calls()
	  << model->cutGenerators(k)->numConsGenerated()
	  << model->cutGenerators(k)->time()
	  << model->cutGenerators(k)->strategy()
	  << CoinMessageEol;
      }
    }

    numT = model->numHeuristics();
    for (k = 0; k < numT; ++k) {
      if ( model->heuristics(k)->calls() > -1) {
	model->blisMessageHandler()->message(DISCO_HEUR_STAT_NODE,
					     model->blisMessages())
	  << index_
	  << model->heuristics(k)->name()
	  << model->heuristics(k)->calls()
	  << model->heuristics(k)->numSolutions()
	  << model->heuristics(k)->time()
	  << model->heuristics(k)->strategy()
	  << CoinMessageEol;
      }
    }
  }

  delete [] currLpSolution;

  if (msgLevel >= 100){
    printf("Final Bound:  %.3f\n", model->solver()->getObjValue());
  }

  if (status_ == AlpsNodeStatusFathomed) {
    // Delete new cuts since no use anymore.
    for (k = 0; k < newNumCons; ++k) {
      delete newConstraints[k];
    }
  }
  delete [] newConstraints;
  delete [] oldConsPos;

  model->isRoot_ = false;

#ifdef DISCO_DEBUG_MORE
  // Debug survived old constraints.
  //int currNumOldCons = model->getNumOldConstraints();
  for (k = 0; k < currNumOldCons; ++k) {
    DcoConstraint *aCon = model->oldConstraints()[k];
    assert(aCon);
    std::cout << "SAVE: DBG: TERM: "
	      << "k=" << k << ", len=" << aCon->getSize()
	      << ", node=" << index_ << std::endl;
  }
#endif

  return returnStatus;
}

//#############################################################################

std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
DcoTreeNode::branch() {

  //------------------------------------------------------
  // Change one var hard bound and record the change in nodedesc:
  // THINK: how about constraint bounds? When to update?
  // TODO: how about other SOS object, etc.?
  //------------------------------------------------------

  AlpsPhase phase = knowledgeBroker_->getPhase();

  double objVal = quality_;

  DcoNodeDesc* childDesc = NULL;

  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
    childNodeDescs;

  DcoModel* model = dynamic_cast<DcoModel*>(desc_->getModel());

  int numCols = model->getNumCols();

#ifdef DISCO_DEBUG_MORE
  // Debug survived old constraints.
  int currNumOldCons = model->getNumOldConstraints();
  for (int k = 0; k < currNumOldCons; ++k) {
    DcoConstraint *aCon = model->oldConstraints()[k];
    assert(aCon);
    std::cout << "BRANCH: DBG: "
	      << "k=" << k << ", len=" << aCon->getSize()
	      << ", node=" << index_ << std::endl;
  }
#endif

  //------------------------------------------------------
  // Get branching object. TODO: Assume integer branching object.
  //------------------------------------------------------

  DcoNodeDesc* thisDesc = dynamic_cast<DcoNodeDesc*>(desc_);

  DcoBranchingObjectType type =
    static_cast<DcoBranchingObjectType>(branchObject_->getType());

  switch(type){

  case DcoBranchingObjectTypeInt:
    {
      DcoBranchObjectInt *branchObject =
	dynamic_cast<DcoBranchObjectInt *>(branchObject_);

      int objInd = branchObject->getObjectIndex();

      double bValue = branchObject->getValue();

      DcoObjectInt *obj =
	dynamic_cast<DcoObjectInt *>(model->objects(objInd));
      int branchVar = obj->columnIndex();

#ifdef DISCO_DEBUG
      if ( (branchVar < 0) || (branchVar >= numCols) ) {
	std::cout << "ERROR: BRANCH(): branchVar = " << branchVar
		  << "; numCols = " << numCols  << std::endl;
	throw CoinError("branch index is out of range",
			"branch", "DcoTreeNode");
      }
#endif

#ifdef DISCO_DEBUG
      printf("BRANCH(): on %d, phase %d\n", branchVar, phase);
      printf("DOWN: lb %g, up %g\n",
	     branchObject->getDown()[0], branchObject->getDown()[1]);
      printf("UP  : lb %g, up %g\n",
	     branchObject->getUp()[0], branchObject->getUp()[1]);
#endif

      //======================================================
      //------------------------------------------------------
      // Create down-branch node description.
      //------------------------------------------------------
      //======================================================

      childDesc = new DcoNodeDesc(model);

      if (phase == AlpsPhaseRampup) {

	//--------------------------------------------------
	// Store a full description since each node will be the
	// root of a subtree.
	// NOTE: this desc must be explicit during rampup.
	//--------------------------------------------------

	int index, k;
	int numModify = -1;
	double value;

	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];

	double *fVarSoftLB = NULL;
	double *fVarSoftUB = NULL;
	int *fVarSoftLBInd = NULL;
	int *fVarSoftUBInd = NULL;

	//--------------------------------------------------
	// Full hard variable bounds.
	//--------------------------------------------------

	numModify = thisDesc->getVars()->lbHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	  index = thisDesc->getVars()->lbHard.posModify[k];
	  assert(index == k);
	  value = thisDesc->getVars()->lbHard.entries[k];
	  fVarHardLB[k] = value;
	  fVarHardLBInd[k] = index;
	}

	numModify = thisDesc->getVars()->ubHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	  index = thisDesc->getVars()->ubHard.posModify[k];
	  assert(index == k);
	  value = thisDesc->getVars()->ubHard.entries[k];
	  fVarHardUB[k] = value;
	  fVarHardUBInd[k] = index;
	}

	// Branching bounds.
	fVarHardLB[branchVar] = branchObject->getDown()[0];
	fVarHardUB[branchVar] = branchObject->getDown()[1];


	childDesc->assignVarHardBound(numCols,
				      fVarHardLBInd,
				      fVarHardLB,
				      numCols,
				      fVarHardUBInd,
				      fVarHardUB);

	//--------------------------------------------------
	// Soft variable bounds.
	//--------------------------------------------------

	int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
	assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
	if (numSoftVarLowers > 0) {
	  fVarSoftLB = new double [numSoftVarLowers];
	  fVarSoftLBInd = new int [numSoftVarLowers];
	  for (k = 0; k < numSoftVarLowers; ++k) {
	    index = thisDesc->getVars()->lbSoft.posModify[k];
	    value = thisDesc->getVars()->lbSoft.entries[k];
	    fVarSoftLB[k] = value;
	    fVarSoftLBInd[k] = index;
	  }
	}

	int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
	assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
	if (numSoftVarUppers > 0) {
	  fVarSoftUB = new double [numSoftVarUppers];
	  fVarSoftUBInd = new int [numSoftVarUppers];
	  for (k = 0; k < numSoftVarUppers; ++k) {
	    index = thisDesc->getVars()->ubSoft.posModify[k];
	    value = thisDesc->getVars()->ubSoft.entries[k];
	    fVarSoftUB[k] = value;
	    fVarSoftUBInd[k] = index;
	  }
	}

#ifdef DISCO_DEBUG_MORE
	// Print soft bounds.
	std::cout << "\nBRANCH: numSoftVarLowers=" << numSoftVarLowers
		  << ", numSoftVarUppers=" << numSoftVarUppers
		  << std::endl;
	for (k = 0; k < numSoftVarLowers; ++k) {
	  std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
		    << fVarSoftLB[k] << std::endl;
	}
	std::cout << "------------------" << std::endl;
	for (k = 0; k < numSoftVarUppers; ++k) {
	  std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
		    << fVarSoftUB[k] << std::endl;
	}
	std::cout << "------------------" << std::endl << std::endl;
#endif

	// Assign it anyway so to transfer ownership of
	// memory(fVarSoftLBInd,etc.)
	childDesc->assignVarSoftBound(numSoftVarLowers,
				      fVarSoftLBInd,
				      fVarSoftLB,
				      numSoftVarUppers,
				      fVarSoftUBInd,
				      fVarSoftUB);

	//--------------------------------------------------
	// Full set of non-core constraints.
	// NOTE: non-core constraints have been saved in description
	//       when process() during ramp-up.
	//--------------------------------------------------

	BcpsObject **tempCons = NULL;
	int tempInt = 0;

	tempInt = thisDesc->getCons()->numAdd;
	if (tempInt > 0) {
	  tempCons = new BcpsObject* [tempInt];
	  for (k = 0; k < tempInt; ++k) {
	    DcoConstraint *aCon = dynamic_cast<DcoConstraint *>
	      (thisDesc->getCons()->objects[k]);

	    assert(aCon);
	    assert(aCon->getSize() > 0);
	    assert(aCon->getSize() < numCols);
	    DcoConstraint *newCon = new DcoConstraint(*aCon);
	    tempCons[k] = newCon;
	  }
	}


#if 0
	else {
	  // No cons or only root cons.
	  tempInt = model->getNumOldConstraints();

	  if (tempInt > 0) {
	    tempCons = new BcpsObject* [tempInt];
	  }
	  for (k = 0; k < tempInt; ++k) {
	    DcoConstraint *aCon = model->oldConstraints()[k];
	    assert(aCon);
	    assert(aCon->getSize() > 0);
	    assert(aCon->getSize() < numCols);
	    DcoConstraint *newCon = new DcoConstraint(*aCon);
	    tempCons[k] = newCon;
	  }
	}
#endif

#ifdef DISCO_DEBUG_MORE
	std::cout << "BRANCH: down: tempInt=" << tempInt <<std::endl;
#endif
	// Fresh desc, safely add.
	childDesc->setAddedConstraints(tempInt, tempCons);
      }else{

	//--------------------------------------------------
	// Relative: Only need to record hard var bound change.
	// NOTE: soft var bound changes are Record after
	// selectBranchObject.
	//--------------------------------------------------

	childDesc->setVarHardBound(1,
				   &branchVar,
				   &(branchObject->getDown()[0]),
				   1,
				   &branchVar,
				   &(branchObject->getDown()[1]));
      }

      childDesc->setBranchedDir(-1);
      childDesc->setBranchedInd(objInd);
      childDesc->setBranchedVal(bValue);

      // Copy warm start.
      CoinWarmStartBasis *ws = thisDesc->getBasis();
      CoinWarmStartBasis * newWs;
      // check whether warm start basis is available
      if (ws==0) {
	newWs = 0;
      }
      else {
	newWs = new CoinWarmStartBasis(*ws);
      }
      childDesc->setBasis(newWs);

      childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					      (childDesc),
					      AlpsNodeStatusCandidate,
					      objVal));

      //======================================================
      //------------------------------------------------------
      // Create up-branch node description.
      //------------------------------------------------------
      //======================================================

      childDesc = new DcoNodeDesc(model);

      if (phase == AlpsPhaseRampup) {

	//--------------------------------------------------
	// Store a full description since each node will be the
	// root of a subtree.
	// NOTE: parent must be explicit during rampup.
	//--------------------------------------------------

	int index, k;
	int numModify = -1;
	double value;

	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];

	double *fVarSoftLB = NULL;
	double *fVarSoftUB = NULL;
	int *fVarSoftLBInd = NULL;
	int *fVarSoftUBInd = NULL;

	//--------------------------------------------------
	// Full hard variable bounds.
	//--------------------------------------------------

	numModify = thisDesc->getVars()->lbHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	  index = thisDesc->getVars()->lbHard.posModify[k];
	  assert(index == k);
	  value = thisDesc->getVars()->lbHard.entries[k];
	  fVarHardLB[k] = value;
	  fVarHardLBInd[k] = index;
	}

	numModify = thisDesc->getVars()->ubHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	  index = thisDesc->getVars()->ubHard.posModify[k];
	  assert(index == k);
	  value = thisDesc->getVars()->ubHard.entries[k];
	  fVarHardUB[k] = value;
	  fVarHardUBInd[k] = index;
	}

	// Branching bounds.
	fVarHardLB[branchVar] = branchObject->getUp()[0];
	fVarHardUB[branchVar] = branchObject->getUp()[1];

	childDesc->assignVarHardBound(numCols,
				      fVarHardLBInd,
				      fVarHardLB,
				      numCols,
				      fVarHardUBInd,
				      fVarHardUB);

	//--------------------------------------------------
	// Soft variable bounds.
	//--------------------------------------------------

	int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
	assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
	if (numSoftVarLowers > 0) {
	  fVarSoftLB = new double [numSoftVarLowers];
	  fVarSoftLBInd = new int [numSoftVarLowers];
	  for (k = 0; k < numSoftVarLowers; ++k) {
	    index = thisDesc->getVars()->lbSoft.posModify[k];
	    value = thisDesc->getVars()->lbSoft.entries[k];
	    fVarSoftLB[k] = value;
	    fVarSoftLBInd[k] = index;
	  }
	}

	int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
	assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
	if (numSoftVarUppers > 0) {
	  fVarSoftUB = new double [numSoftVarUppers];
	  fVarSoftUBInd = new int [numSoftVarUppers];
	  for (k = 0; k < numSoftVarUppers; ++k) {
	    index = thisDesc->getVars()->ubSoft.posModify[k];
	    value = thisDesc->getVars()->ubSoft.entries[k];
	    fVarSoftUB[k] = value;
	    fVarSoftUBInd[k] = index;
	  }
	}

	// Assign it anyway so to transfer ownership of
	// memory(fVarSoftLBInd,etc.)
	childDesc->assignVarSoftBound(numSoftVarLowers,
				      fVarSoftLBInd,
				      fVarSoftLB,
				      numSoftVarUppers,
				      fVarSoftUBInd,
				      fVarSoftUB);

	//--------------------------------------------------
	// Full set of non-core constraints.
	// NOTE: non-core constraints have been saved in description
	//       when process() during ramp-up.
	//--------------------------------------------------

	BcpsObject **tempCons = NULL;
	int tempInt = 0;

	tempInt = thisDesc->getCons()->numAdd;
	if (tempInt > 0) {
	  tempCons = new BcpsObject* [tempInt];

	  for (k = 0; k < tempInt; ++k) {
	    DcoConstraint *aCon = dynamic_cast<DcoConstraint *>
	      (thisDesc->getCons()->objects[k]);

	    assert(aCon);
	    assert(aCon->getSize() > 0);
	    assert(aCon->getSize() <= numCols);
	    DcoConstraint *newCon = new DcoConstraint(*aCon);
	    tempCons[k] = newCon;
	  }
	}

#if 0
	else {
	  // No cons or only root cons.
	  tempInt = model->getNumOldConstraints();
	  if (tempInt > 0) {
	    tempCons = new BcpsObject* [tempInt];
	  }
	  for (k = 0; k < tempInt; ++k) {
	    DcoConstraint *aCon = model->oldConstraints()[k];
	    assert(aCon);
	    assert(aCon->getSize() > 0);
	    assert(aCon->getSize() <= numCols);
	    DcoConstraint *newCon = new DcoConstraint(*aCon);
	    tempCons[k] = newCon;
	  }
	}
#endif

#ifdef DISCO_DEBUG_MORE
	std::cout << "BRANCH: up: tempInt=" << tempInt <<std::endl;
#endif
	// Fresh desc, safely add.
	childDesc->setAddedConstraints(tempInt, tempCons);
      }else{

	//--------------------------------------------------
	// Relative: Only need to record hard var bound change.
	// NOTE: soft var bound changes are Record after selectBranchObject.
	//--------------------------------------------------

	childDesc->setVarHardBound(1,
				   &branchVar,
				   &(branchObject->getUp()[0]),
				   1,
				   &branchVar,
				   &(branchObject->getUp()[1]));
      }

      childDesc->setBranchedDir(1);
      childDesc->setBranchedInd(objInd);
      childDesc->setBranchedVal(bValue);

      // Copy warm start.
      CoinWarmStartBasis * newWs2;
      // check whether warm start basis is available
      if (ws==0) {
	newWs2 = 0;
      }
      else {
	newWs2 = new CoinWarmStartBasis(*ws);
      }
      childDesc->setBasis(newWs2);
      childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					      (childDesc),
					      AlpsNodeStatusCandidate,
					      objVal));

      // Change node status to branched.
      status_ = AlpsNodeStatusBranched;
    }

    break;

  case DcoBranchingObjectTypeBilevel :
    {
      DcoBranchObjectBilevel *branchObject =
	dynamic_cast<DcoBranchObjectBilevel *>(branchObject_);
      std::deque<int> *branchingSet = branchObject->getBranchingSet();

      CoinWarmStartBasis *ws = thisDesc->getBasis();
      std::deque<int>::iterator ptr1, ptr2;
      int size = static_cast<int>(branchingSet->size());
      int *indices = new int[size];
      double *values = new double[size];
      values[0] = 1;
      int i;
      for (i = 0, ptr1 = branchingSet->begin();
	   ptr1 != branchingSet->end(); i++, ptr1++){
	indices[i] = *ptr1;
	values[i] = 1;
	childDesc = new DcoNodeDesc(model);
	childDesc->setVarHardBound(i+1, indices, values,
				   i+1, indices, values);
	CoinWarmStartBasis *newWs = new CoinWarmStartBasis(*ws);
	childDesc->setBasis(newWs);
	childNodeDescs.push_back(CoinMakeTriple(
						static_cast<AlpsNodeDesc *>(childDesc),
						AlpsNodeStatusCandidate,
						quality_));
	values[i] = 0;
      }
      delete[] indices;
      delete[] values;
    }

    break;

  default :

    std::cout << "Unknown branching object type" << std::endl;
  }
  return childNodeDescs;
}

//#############################################################################

/* FIXME: need rewrite from scratch */
/* 0: find a branch var, -1 no branch var (should not happen) */

int DcoTreeNode::selectBranchObject(DcoModel *model,
				    bool& foundSol,
				    int numPassesLeft) {
  int bStatus = 0;
  BcpsBranchStrategy * strategy = 0;
  AlpsPhase phase = knowledgeBroker_->getPhase();
  if(branchObject_) {
    delete branchObject_;
    branchObject_ = NULL;
  }
  //------------------------------------------------------
  // Get branching strategy.
  //------------------------------------------------------
  if (phase == AlpsPhaseRampup) {
    strategy = model->rampUpBranchStrategy();
  }
  else {
    strategy = model->branchStrategy();
  }

  if (!strategy) {
    throw CoinError("No branch strategy.", "selectBranchObject()",
		    "DcoTreeNode");
  }
  //------------------------------------------------------
  // Create branching object candidates.
  //-----------------------------------------------------
  bStatus = strategy->createCandBranchObjects(numPassesLeft,
					      model->getCutoff());
  //------------------------------------------------------
  // Select the best branching objects.
  //-----------------------------------------------------
  if (bStatus >= 0) {
    branchObject_ = strategy->bestBranchObject();
    if (branchObject_==0) {
      throw CoinError("No branch object created.", "selectBranchObject()",
		      "DcoTreeNode");
    }
  }
  if (!model->branchStrategy()) {
    delete strategy;
  }
  return bStatus;
}

//#############################################################################

int
DcoTreeNode::bound(BcpsModel *model)
{

  DcoLpStatus status = DcoLpStatusUnknown;
  DcoModel *m = dynamic_cast<DcoModel *>(model);
  // Bounding
  m->solver()->resolve();
  if (m->solver()->isAbandoned()) {
    status = DcoLpStatusAbandoned;
  }
  else if (m->solver()->isProvenOptimal()) {
    // todo(aykut) if obj val is greater than 1e+30 we consider problem
    // infeasible. This is due to our cut generation when the problem is
    // infeasible. We add a high lower bound (1e+30) to the objective
    // function. This can be improved.
    if (m->solver()->getObjValue()>=1e+30) {
      status = DcoLpStatusPrimalInfeasible;
    }
    else {
      status = DcoLpStatusOptimal;
      DcoNodeDesc *desc = dynamic_cast<DcoNodeDesc*>(desc_);
      double objValue = m->solver()->getObjValue() *
	m->solver()->getObjSense();
      int dir = desc->getBranchedDir();
      if (dir != 0) {
	double objDeg = objValue - quality_;
	int objInd = desc->getBranchedInd();
	double lpX = desc->getBranchedVal();
	DcoObjectInt *intObject =
	  dynamic_cast<DcoObjectInt *>(m->objects(objInd));
	intObject->pseudocost().update(dir, objDeg, lpX);
	m->setSharedObjectMark(intObject->getObjectIndex());
      }
      // Update quality of this nodes.
      quality_ = objValue;
    }
  }
  else if (m->solver()->isProvenPrimalInfeasible()) {
    status = DcoLpStatusPrimalInfeasible;
  }
  else if (m->solver()->isProvenDualInfeasible()) {
    status = DcoLpStatusDualInfeasible;
  }
  else if (m->solver()->isPrimalObjectiveLimitReached()) {
    status = DcoLpStatusPrimalObjLim;
  }
  else if (m->solver()->isDualObjectiveLimitReached()) {
    status = DcoLpStatusDualObjLim;
  }
  else if (m->solver()->isIterationLimitReached()) {
    status = DcoLpStatusIterLim;
  }
  else {
    std::cout << "UNKNOWN LP STATUS" << std::endl;
    assert(0);
  }
  return status;
}

//#############################################################################
int DcoTreeNode::installSubProblem(BcpsModel *m) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  int i, k;
  int index;
  double value;
  DcoModel *model = dynamic_cast<DcoModel *>(m);
  assert(model);
  DcoNodeDesc *desc = dynamic_cast<DcoNodeDesc*>(desc_);
  int numModify = 0;
  int numCoreCols = model->getNumCoreVariables();
  int numCoreRows = model->getNumCoreConstraints();
  int numCols = model->solver()->getNumCols();
  int numRows = model->solver()->getNumRows();
  //double *varSoftLB = NULL;
  //double *varSoftUB = NULL;
  double *colHardLB = NULL;
  double *colHardUB = NULL;
  //double *conSoftLB = NULL;
  //double *conSoftUB = NULL;
  //double *conHardLB = NULL;
  //double *conHardUB = NULL;
  double *startColLB = model->startVarLB();
  double *startColUB = model->startVarUB();
  double *startRowLB = model->startConLB();
  double *startRowUB = model->startConUB();
  CoinFillN(startColLB, numCoreCols, -ALPS_DBL_MAX);
  CoinFillN(startColUB, numCoreCols, ALPS_DBL_MAX);
  CoinFillN(startRowLB, numCoreRows, -ALPS_DBL_MAX);
  CoinFillN(startRowUB, numCoreRows, ALPS_DBL_MAX);
  int numOldRows = 0;
  int tempInt = 0;
  DcoConstraint *aCon = NULL;
  int nodeID = -1;
  nodeID = getIndex();
  //std::cout << "nodeID=" << nodeID << std::endl;
  AlpsPhase phase = knowledgeBroker_->getPhase();
  //======================================================
  // Restore subproblem:
  //  1. Remove noncore var/con
  //  2. Travel back to root and correct differencing to
  //     full var/con bounds into model->startXXX
  //  3. Set col bounds
  //  4. Set row bounds (is this necessary?)
  //  5. Add contraints except cores
  //  6. Add variables except cores
  //  7. Set basis (should not need modify)
  //======================================================
  //------------------------------------------------------
  // Remove old constraints from lp solver.
  //------------------------------------------------------
  int numDelRows = numRows - numCoreRows;
  if (numDelRows > 0) {
    int *indices = new int [numDelRows];
    if (indices == NULL) {
      throw CoinError("Out of memory", "installSubProblem", "DcoTreeNode");
    }
    for (i = 0; i < numDelRows; ++i) {
      indices[i] = numCoreRows + i;
    }
    model->solver()->deleteRows(numDelRows, indices);
    delete [] indices;
    indices = NULL;
  }
  //--------------------------------------------------------
  // Travel back to a full node, then collect diff (add/rem col/row,
  // hard/soft col/row bounds) from the node full to this node.
  //----------------------------
  // Note: if we store full set of logic/agorithm col/row, then
  //       no diff are needed for col/row
  //--------------------------------------------------------
  //--------------------------------------------------------
  // Collect differencing bounds. Branching bounds of this node
  // are ALSO collected.
  //--------------------------------------------------------
  DcoNodeDesc* pathDesc = NULL;
  AlpsTreeNode *parent = parent_;
  /* First push this node since it has branching hard bounds.
     NOTE: during rampup, this desc has full description when branch(). */
  std::vector<AlpsTreeNode*> leafToRootPath;
  leafToRootPath.push_back(this);
  if (phase != AlpsPhaseRampup) {
    while(parent) {
      leafToRootPath.push_back(parent);
      if (parent->getExplicit()) {
	// Reach an explicit node, then stop.
	break;
      }
      else {
	parent = parent->getParent();
      }
    }
  }
  //------------------------------------------------------
  // Travel back from this node to the explicit node to
  // collect full description.
  //------------------------------------------------------
  for(i = static_cast<int> (leafToRootPath.size() - 1); i > -1; --i) {
    //--------------------------------------------------
    // NOTE: As away from explicit node, bounds become
    //       tighter and tighter.
    //--------------------------------------------------
    pathDesc = dynamic_cast<DcoNodeDesc*>((leafToRootPath.at(i))->
					  getDesc());
    colHardLB = pathDesc->getVars()->lbHard.entries;
    colHardUB = pathDesc->getVars()->ubHard.entries;
    //--------------------------------------------------
    // Adjust bounds according to hard var lb/ub.
    // If rampup or explicit, collect hard bounds so far.
    //--------------------------------------------------
    numModify = pathDesc->getVars()->lbHard.numModify;
    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->lbHard.posModify[k];
      value = pathDesc->getVars()->lbHard.entries[k];
      // Hard bounds do NOT change according to soft bounds, so
      // here need CoinMax.
      startColLB[index] = CoinMax(startColLB[index], value);
    }
    numModify = pathDesc->getVars()->ubHard.numModify;
    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->ubHard.posModify[k];
      value = pathDesc->getVars()->ubHard.entries[k];
      startColUB[index] = CoinMin(startColUB[index], value);
    }
    //--------------------------------------------------
    // Adjust bounds according to soft var lb/ub.
    // If rampup or explicit, collect soft bounds so far.
    //--------------------------------------------------
    numModify = pathDesc->getVars()->lbSoft.numModify;
    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->lbSoft.posModify[k];
      value = pathDesc->getVars()->lbSoft.entries[k];
      startColLB[index] = CoinMax(startColLB[index], value);
    }
    numModify = pathDesc->getVars()->ubSoft.numModify;

    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->ubSoft.posModify[k];
      value = pathDesc->getVars()->ubSoft.entries[k];
      startColUB[index] = CoinMin(startColUB[index], value);
    }
    //--------------------------------------------------
    // TODO: Modify hard/soft row lb/ub.
    //--------------------------------------------------
    //--------------------------------------------------
    // Collect active non-core constraints at parent.
    //--------------------------------------------------
    //----------------------------------------------
    // First collect all generated cuts, then remove
    // deleted.
    //----------------------------------------------
    tempInt = pathDesc->getCons()->numAdd;
    int maxOld = model->getOldConstraintsSize();
    for (k = 0; k < tempInt; ++k) {
      aCon = dynamic_cast<DcoConstraint *>
	(pathDesc->getCons()->objects[k]);
      assert(aCon);
      assert(aCon->getSize() > 0);
      assert(aCon->getSize() < 100000);
      (model->oldConstraints())[numOldRows++] = aCon;
      if (numOldRows >= maxOld) {
	// Need resize
	maxOld *= 2;
	DcoConstraint **tempCons = new DcoConstraint* [maxOld];
	memcpy(tempCons,
	       model->oldConstraints(),
	       numOldRows * sizeof(DcoConstraint *));
	model->delOldConstraints();
	model->setOldConstraints(tempCons);
	model->setOldConstraintsSize(maxOld);
      }
    }
    //----------------------------------------------
    // Remove those deleted.
    // NOTE: model->oldConstraints_ stores all previously
    // generated active constraints at parent.
    //----------------------------------------------
    tempInt = pathDesc->getCons()->numRemove;
    if (tempInt > 0) {
      int tempPos;
      int *tempMark = new int [numOldRows];
      CoinZeroN(tempMark, numOldRows);
      for (k = 0; k < tempInt; ++k) {
	tempPos = pathDesc->getCons()->posRemove[k];
	tempMark[tempPos] = 1;
      }
      tempInt = 0;
      for (k = 0; k < numOldRows; ++k) {
	if (tempMark[k] != 1) {
	  // Survived.
	  (model->oldConstraints())[tempInt++]=
	    (model->oldConstraints())[k];
	}
      }
      if (tempInt + pathDesc->getCons()->numRemove != numOldRows) {
	std::cout << "INSTALL: tempInt=" << tempInt
		  <<", numRemove="<<pathDesc->getCons()->numRemove
		  << ", numOldRows=" << numOldRows << std::endl;
	assert(0);
      }
      // Update number of old non-core constraints.
      numOldRows = tempInt;
      delete [] tempMark;
    }
  } // EOF leafToRootPath.
  //--------------------------------------------------------
  // Clear path vector.
  //--------------------------------------------------------
  leafToRootPath.clear();
  assert(leafToRootPath.size() == 0);
  //--------------------------------------------------------
  // Adjust column bounds in lp solver
  //--------------------------------------------------------
  for(i = 0; i < numCols; ++i) {
    model->solver()->setColBounds(i, startColLB[i], startColUB[i]);
  }
  //--------------------------------------------------------
  // TODO: Set row bounds
  //--------------------------------------------------------
  //--------------------------------------------------------
  // Add old constraints, which are collect from differencing.
  //--------------------------------------------------------
  // If removed cuts due to local cuts.
  model->setNumOldConstraints(numOldRows);
  if (numOldRows > 0) {
    const OsiRowCut ** oldOsiCuts = new const OsiRowCut * [numOldRows];
    for (k = 0; k < numOldRows; ++k) {
      OsiRowCut * acut = (model->oldConstraints()[k])->createOsiRowCut();
      oldOsiCuts[k] = acut;
    }
    model->solver()->applyRowCuts(numOldRows, oldOsiCuts);
    for (k = 0; k < numOldRows; ++k) {
      delete oldOsiCuts[k];
    }
    delete [] oldOsiCuts;
    oldOsiCuts = NULL;
  }
  //--------------------------------------------------------
  // Add parent variables, which are collect from differencing.
  //--------------------------------------------------------
  //--------------------------------------------------------
  // Set basis
  //--------------------------------------------------------
  CoinWarmStartBasis *pws = desc->getBasis();
  if (pws != NULL) {
    model->solver()->setWarmStart(pws);
  }
  return status;
}

//#############################################################################

void
DcoTreeNode::getViolatedConstraints(DcoModel *model,
				    const double *LpSolution,
				    BcpsConstraintPool & conPool)
{
  int k;
  int numCons = model->constraintPoolReceive()->getNumConstraints();
  DcoConstraint *blisCon = NULL;
  std::vector<DcoConstraint *> conVector;

  // Check violation and move voilated constraints to conPool
  for (k = 0; k < numCons; ++k) {
    blisCon = dynamic_cast<DcoConstraint *>(conPool.getConstraint(k));

    if (blisCon->violation(LpSolution) > ALPS_SMALL_4) {
      conPool.addConstraint(blisCon);
    }
    else {
      conVector.push_back(blisCon);
    }
  }

  if (numCons > 0) {
    std::cout << "Has constraints " << numCons
	      << "; violated " << numCons-conVector.size()
	      << std::endl;
  }

  if ((int)conVector.size() != numCons) {
    // There violated constraints. Remove them from conPool.
    conPool.clear();
    numCons = static_cast<int>(conVector.size());
    for (k = 0; k < numCons; ++k) {
      conPool.addConstraint(conVector[k]);
    }
    conVector.clear();
  }
}

//#############################################################################

int
DcoTreeNode::generateConstraints(DcoModel *model,BcpsConstraintPool &conPool)
{
  int i, j, numCGs;
  DcoLpStatus status = DcoLpStatusOptimal;
  int preNumCons = 0;
  int newNumCons = 0;
  DcoCutStrategy strategy = DcoCutStrategyRoot;

  // Only autmatic stategy has depth limit.
  int maxConstraintDepth = 20;

  bool mustResolve = false;

  double genCutTime;

  numCGs = model->numCutGenerators();
  int ipm_fails = 0;
  // whether ipm cut generator exists and before oa cut generator
  int ipm_exists = 0;
  for (i = 0 ; i < numCGs; ++i) {
    std::string name = model->cutGenerators(i)->name();
    if (name.compare("ipm_gen")==0) {
      // ipm is before oa cut generator, it exists.
      ipm_exists = 1;
      break;
    }
    else if (name.compare("oa_gen")==0) {
      // we get oa before ipm
      break;
    }
  }
  for (i = 0 ; i < numCGs; ++i) {

    //----------------------------------------------------
    // Check if call this generator.
    //----------------------------------------------------

    strategy =  model->cutGenerators(i)->strategy();

    bool useThisCutGenerator = false;
    if (strategy == DcoCutStrategyNone) {
      useThisCutGenerator = false;
    }
    else if (strategy == DcoCutStrategyRoot) {
      if (model->isRoot_ && (index_ == 0)) useThisCutGenerator = true;
    }
    else if (strategy == DcoCutStrategyAuto) {
      if (depth_ < maxConstraintDepth) {
	if (!diving_ || model->isRoot_) useThisCutGenerator = true;
      }
    }
    else if (strategy == DcoCutStrategyPeriodic) {
      // Num of nodes is set at the beginning of process().
      if ((model->getNumNodes()-1) %
	  model->cutGenerators(i)->cutGenerationFreq() == 0) {
	useThisCutGenerator = true;
      }
    }
    else {
      throw CoinError("Unknown cut generation strategy",
		      "generateConstraints", "DcoTreeNode");
    }
    // todo(aykut)
    // if cut generator is ipm, check if the solution is integer
    // use it if the solution is integer feasible only.
#if defined(__OA__)
    DcoSolution *ipSol = NULL;
    int numIntInfs;
    int numObjInfs;
    std::string name = model->cutGenerators(i)->name();
    if (name.compare("ipm_gen")==0) {
      ipSol = model->feasibleSolution(numIntInfs, numObjInfs);
      if (numIntInfs==0) {
	if (numObjInfs==0) {
	  std::cerr << "Current solution is feasible! Why are we trying to generate a cut?" << std::endl;
	  throw std::exception();
	}
	else {
	  // integer solution does not satisfy conic constraints
	  // use ipm cut generator
	  useThisCutGenerator = true;
	}
      }
      else {
	// solution is not integer
	useThisCutGenerator = false;
      }
    }
    else if (name.compare("oa_gen")==0) {
      // if ipm does not exists generate cuts using oa
      if (ipm_exists) {
	// what if ipm is not added as a generator?
	if (ipm_fails) {
	  // ipm did not yield any cuts, use CglConicOA to generate cuts.
	  useThisCutGenerator = true;
	}
	else {
	  useThisCutGenerator = false;
	}
      }
      else {
	useThisCutGenerator = true;
      }
    }
    // generate ipm cuts 10% of the time
    if (name.compare("ipm_gen")==0 && useThisCutGenerator) {
      // generate Random Number
      int rand_number = rand()%100;
      if (rand_number > -1) {
	ipm_fails = 2;
	useThisCutGenerator = false;
      }
    }
#endif

#if 0
    std::cout<<"CUTGEN: " << model->cutGenerators(i)->name()
	     <<": useThisCutGenerator ="<<useThisCutGenerator
	     <<", diving =" << diving_
	     << ", strategy =" << strategy
	     << ", num of nodes =" << model->getNumNodes()
	     <<std::endl;
#endif
    // int numIntInfs;
    // int numObjInfs;
    // DcoSolution const * ipSol = model->feasibleSolution(numIntInfs, numObjInfs);
    // std::string ipm_strat("ipm_gen");
    // std::cout << "cut generator " << model->cutGenerators(0)->name() << std::endl;
    // if (model->cutGenerators(0)->name()==ipm_strat) {
    //   std::cout << "***** set ipm cut generator to false. *****" << std::endl;
    //   useThisCutGenerator = false;
    //   if (!numIntInfs && numObjInfs) {
    //	std::cout << "***** number of integer feasibility is " << numIntInfs;
    //	std::cout << " cut generator is set to true. *****" << std::endl;
    //	useThisCutGenerator = true;
    //   }
    //   if (model->boundingPass_>=2) {
    //	useThisCutGenerator = false;
    //   }
    // }

    //----------------------------------------------------
    // Generator constraints.
    //----------------------------------------------------
    if (useThisCutGenerator) {
      preNumCons = conPool.getNumConstraints();
      genCutTime = CoinCpuTime();
      // Call constraint generator
      mustResolve = model->cutGenerators(i)->generateConstraints(conPool);
#if defined(__OA__)
      if (conPool.getNumConstraints()==0 and name.compare("ipm_gen")==0) {
	ipm_fails = 1;
      }
#endif
      genCutTime = CoinCpuTime() - genCutTime;
      // Statistics
      model->cutGenerators(i)->addTime(genCutTime);
      model->cutGenerators(i)->addCalls(1);
      newNumCons = conPool.getNumConstraints() - preNumCons;
      if (newNumCons == 0) {
	model->cutGenerators(i)->addNoConsCalls(1);
      }
      else {
	// Reset to 0
	int noCuts = model->cutGenerators(i)->noConsCalls();
	model->cutGenerators(i)->addNoConsCalls(-noCuts);
	model->cutGenerators(i)->addNumConsGenerated(newNumCons);
      }
      if (mustResolve) {
	// TODO: Only probing will return ture.
	status = static_cast<DcoLpStatus> (bound(model));
	if (status == DcoLpStatusOptimal) {
#ifdef DISCO_DEBUG
	  std::cout << "CUTGEN: after probing, this node survived."
		    << std::endl;
#endif
	}
	else {
#ifdef DISCO_DEBUG
	  std::cout<<"CUTGEN: after probing, this node can fathomed."
		   << std::endl;
#endif
	  break;
	}
      }
      //------------------------------------------------
      // Modify control.
      // NOTE: only modify if user choose automatic.
      //------------------------------------------------
      if (model->getCutStrategy() == DcoCutStrategyNone) {
	for (j = 0; j < numCGs; ++j) {
	  strategy =  model->cutGenerators(j)->strategy();
	  if (strategy != DcoCutStrategyNone) {
	    break;
	  }
	}
	if (j == numCGs) {
	  model->setCutStrategy(DcoCutStrategyNone);
	}
      }
    }
  }
  return status;
}

//#############################################################################

DcoReturnStatus
DcoTreeNode::applyConstraints(DcoModel *model,
			      const double *solution,
			      BcpsConstraintPool & conPool)
{
  DcoReturnStatus status = DcoReturnStatusOk;
  int i, k;

  int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);

  int numRowCuts = conPool.getNumConstraints();

  int numToAdd = numRowCuts;
  int numAdded = 0;

  if (numRowCuts > 0) {
    DcoParams * DcoPar = model->DcoPar();
    double scaleConFactor = DcoPar->entry(DcoParams::scaleConFactor);

    if (numToAdd > 0) {

      DcoConstraint *blisCon = NULL;

      if (msgLevel > 100) {
	printf("\nAPPLYCUT: Select cuts to be added in LP from %d candidates\n",
	       numRowCuts);
      }

      int numRowsNow = model->solver()->getNumRows();
      int numCols = model->solver()->getNumCols();
      CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
	(model->solver()->getWarmStart());

      // Tranform constraints to Osi cut so that easily add them to LP.
      const OsiRowCut ** addedCuts = new const OsiRowCut * [numToAdd];

      for (i = 0 ; i < numToAdd ; ++i) {
	bool keep = true;
	blisCon = dynamic_cast<DcoConstraint *>(conPool.getConstraint(i));

	//------------------------------------------
	// Remove:
	//  - empty cuts
	//  - dense cuts
	//  - bad scaled cuts
	//  - weak cuts
	//  - parallel cuts
	//------------------------------------------

	int length = blisCon->getSize();
	const double *elements = blisCon->getValues();
	const int *indices = blisCon->getIndices();

	bool check = true;
	while (check) { // while is used to turn off checking.
	  //--------------------------------------
	  // Empty.
	  //--------------------------------------

	  if (length <= 0) {
	    keep = false;

#if 0
	    std::cout << "APPLYCUT: A empty cut." << std::endl;
#endif
	    break;
	  }

	  //--------------------------------------
	  // Dense.
	  //--------------------------------------

	  if(length > model->getDenseConCutoff()){
	    keep = false;
	    if (msgLevel > 100) {
	      std::cout << "APPLYCUT: Discard a dense cut. length = "
			<< length << ", cutoff = "
			<< model->getDenseConCutoff() << std::endl;
	    }
	    break;
	  }

	  //--------------------------------------
	  // Compuate scale factor.
	  //--------------------------------------

	  int index;
	  double activity = 0.0;

	  double maxElem = 0.0;
	  double minElem = ALPS_DBL_MAX;
	  double scaleFactor;

	  for (k = 0; k < length; ++k) {
	    // todo(aykut) skip if the value of the element is 0.
	    // should this be an exact 0.0?
	    if (elements[k]==0.0) {
	      continue;
	    }
	    if (fabs(elements[k]) > maxElem) {
	      maxElem = fabs(elements[k]);
	    }
	    if (fabs(elements[k]) < minElem) {
	      minElem = fabs(elements[k]);
	    }
	    index = indices[k];
	    activity += elements[k] * solution[index];
	  }
	  if(minElem != 0.0) {
	    scaleFactor = maxElem/minElem;
	  }
	  else {
	    assert(0);
	    scaleFactor = ALPS_DBL_MAX;
	  }

#ifdef DISCO_DEBUG
	  std::cout << "APPLYCUT: scaleFactor=" << scaleFactor
		    << ", maxElem=" << maxElem
		    << ", minElem=" << minElem << std::endl;
#endif
	  if (scaleFactor > scaleConFactor) {
	    if (msgLevel > 100) {
	      std::cout<< "APPLYCUT: Discard a bad scaled cut"
		       << std::endl;
	    }
	    keep = false;
	    break;
	  }

	  //--------------------------------------
	  // Weak.
	  //--------------------------------------

	  double rowLower = CoinMax(blisCon->getLbHard(),
				    blisCon->getLbSoft());
	  double rowUpper = CoinMin(blisCon->getUbHard(),
				    blisCon->getUbSoft());
	  double violation = -9.87; // Any negative number is OK

	  if (rowLower > -ALPS_INFINITY) {
	    violation = rowLower - activity;
	  }
	  if (rowUpper < ALPS_INFINITY) {
	    violation = CoinMax(violation, activity-rowUpper);
	  }

	  if (violation < 1.0e-6) {
	    // Found a weak cuts.
	    if (msgLevel > 100) {
	      std::cout<< "WARNING: APPLYCUT: applied a weak cut, violation="
		       << violation << std::endl;
	    }
	    //keep = false;
	    break;
	  }

	  //--------------------------------------
	  // Parallel cuts.
	  //--------------------------------------

	  bool paral = parallel(model,
				conPool,
				i,
				blisCon);
	  if (paral) {
	    if (msgLevel > 100) {
	      std::cout<< "APPLYCUT: Discard a parallel cut"
		       << std::endl;
	    }
	    keep = false;
	    break;
	  }

	  //--------------------------------------
	  // Check once and stop.
	  //--------------------------------------

	  check = false;
	}//while

	if (keep) {
	  addedCuts[numAdded++] = blisCon->createOsiRowCut();
	}
	else {
	  conPool.deleteConstraint(i);
	  --i;
	  --numToAdd;
	}
      }

      assert(numToAdd == numAdded);

      if (msgLevel > 100) {
	printf("APPLYCUT: After selecting, added %d cuts to LP and discared %d cuts\n",
	       numAdded, numRowCuts - numAdded);
      }

      //----------------------------------------------
      // Add cuts to lp and adjust basis.
      //----------------------------------------------

      if (numAdded > 0) {
	model->solver()->applyRowCuts(numAdded, addedCuts);
	ws->resize(numRowsNow + numToAdd, numCols);
	for (i = 0 ; i < numToAdd; ++i) {
	  ws->setArtifStatus(numRowsNow + i,
			     CoinWarmStartBasis::basic);
	}
	if (model->solver()->setWarmStart(ws) == false) {
	  throw CoinError("Fail setWarmStart() after cut installation.",
			  "applyConstraints","DcoTreeNode");
	}

	for (k = 0; k < numAdded; ++k) {
	  delete addedCuts[k];
	}
      }

      delete [] addedCuts;
      delete ws;
    }
  }

  return status;
}

//#############################################################################

DcoReturnStatus
DcoTreeNode::reducedCostFix(DcoModel *model)
{
  int i, var;
  DcoReturnStatus status = DcoReturnStatusOk;

  int numFixedUp = 0;
  int numFixedDown = 0;
  int numTighten = 0;

  double movement;
  double newBound;
  double boundDistance;
  double dj;

  int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);

  const double *lb = model->solver()->getColLower();
  const double *ub = model->solver()->getColUpper();
  const double *solution = model->solver()->getColSolution();
  const double *reducedCost = model->solver()->getReducedCost();

  double cutoff = model->getCutoff();

  if (cutoff >= ALPS_OBJ_MAX) return status;

  double lpObjValue = model->solver()->getObjValue() *
    model->solver()->getObjSense();
  double epInt = 1.0e-5;

  int numIntegers = model->getNumIntObjects();
  const int *intIndices = model->getIntColIndices();

  for (i = 0; i < numIntegers; ++i) {
    var = intIndices[i];

    dj = reducedCost[var];

    if (fabs(dj) < epInt) continue;

    boundDistance = ub[var] - lb[var];
    if (boundDistance < epInt) continue;

    movement = floor((cutoff - lpObjValue) / fabs(dj));

    if (solution[var] > ub[var] - epInt) {
      /* At upper bound */
      if (movement < boundDistance) {
	/* new lower bound. If movement is 0, then fix. */
	newBound = ub[var] - movement;
	newBound = CoinMin(newBound, ub[var]);

	if (msgLevel > 300) {
	  printf("RED-FIX: dj %g, lb %.10g, ub %.10g, newBound %.10g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
	}

	if (movement <= ALPS_ZERO) {
	  ++numFixedUp;
	}
	else if (newBound < ub[var]){
	  ++numTighten;
	}
	model->solver()->setColLower(var, newBound);
      }
    }
    else if (solution[var] < lb[var] + epInt) {
      /* At lower bound */
      if (movement < boundDistance) {
	newBound = lb[var] + movement;
	newBound = CoinMax(newBound, lb[var]);

	if (msgLevel > 300) {
	  printf("RED-FIX: dj %g, lb %g, ub %g, newBound %g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
	}

	if (movement <= ALPS_ZERO) {
	  ++numFixedDown;
	}
	else if(newBound > lb[var] ){
	  ++numTighten;
	}
	/* new upper bound. If movement is 0, then fix. */
	model->solver()->setColUpper(var, newBound);
      }
    }
  }

  //int change = numFixedUp + numFixedDown + numTighten;
  //model->reducedCostFixed_ += change;

  if (msgLevel > 200) {
    if (numFixedUp > 0 || numFixedDown > 0 || numTighten > 0) {
      printf("reducedCostFix: numFixedUp = %d, numFixedDown = %d, numTighten %d\n", numFixedUp, numFixedDown, numTighten);
    }
  }

  return status;
}

//#############################################################################

AlpsEncoded*
DcoTreeNode::encode() const
{
#ifdef DISCO_DEBUG
  std::cout << "DcoTreeNode::encode()--start to encode node "
	    << index_ << std::endl;
#endif

  AlpsReturnStatus status = AlpsReturnStatusOk;

  // NOTE: "AlpsKnowledgeTypeNode" is used as type name.
  AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeNode);

  // Encode decription.
  status = desc_->encode(encoded);

  // Encode Alps portion.
  status = encodeAlps(encoded);

  // Encode Bcps portion.
  status = encodeBcps(encoded);

  // Nothing to encode for Dco portion.

  return encoded;
}

//#############################################################################

AlpsKnowledge*
DcoTreeNode::decode(AlpsEncoded& encoded) const
{
  AlpsReturnStatus status = AlpsReturnStatusOk;
  DcoTreeNode* treeNode = NULL;

  DcoModel *model = dynamic_cast<DcoModel*>(desc_->getModel());

  //------------------------------------------------------
  // Unpack decription.
  //------------------------------------------------------

  AlpsNodeDesc* nodeDesc = new DcoNodeDesc(model);
  status = nodeDesc->decode(encoded);

  //------------------------------------------------------
  // Unpack node.
  //------------------------------------------------------

  // Unpack Alps portion.
  treeNode = new DcoTreeNode(nodeDesc);
  nodeDesc = NULL;

  treeNode->decodeAlps(encoded);

  // Unpack Bcps portion.
  int type = 0;
  encoded.readRep(type);
  if (type == DcoBranchingObjectTypeInt) {
    // branchObject_ is simple integer.
    DcoBranchObjectInt *bo = new DcoBranchObjectInt();
    status = bo->decode(encoded);

    // Set bo in treeNode.
    treeNode->setBranchObject(bo);
    bo = NULL;
  }

  // Nothing to unpack for Dco portion.

  return treeNode;
}

//#############################################################################
void DcoTreeNode::convertToExplicit() {
  if (explicit_) {
    return;
  }
  // Convert to explicit
  explicit_ = 1;
  DcoModel* model = dynamic_cast<DcoModel*>(desc_->getModel());
  DcoNodeDesc *desc = dynamic_cast<DcoNodeDesc *>(desc_);
  DcoConstraint *aCon = NULL;
  int numCols = model->solver()->getNumCols();
  int i, k, index;
  int tempInt;
  int numModify = 0;
  int numSoftVarLowers = 0;
  int numSoftVarUppers = 0;
  double value;
  double *fVarHardLB = new double [numCols];
  double *fVarHardUB = new double [numCols];
  int *fVarHardLBInd = new int [numCols];
  int *fVarHardUBInd = new int [numCols];
  double *fVarSoftLB = new double [numCols];
  double *fVarSoftUB = new double [numCols];
  int *fVarSoftLBInd = new int [numCols];
  int *fVarSoftUBInd = new int [numCols];
  for (k = 0; k < numCols; ++k) {
    fVarSoftLB[k] = ALPS_DBL_MAX;
    fVarSoftUB[k] = -ALPS_DBL_MAX;
    fVarHardLB[k] = ALPS_DBL_MAX;
    fVarHardUB[k] = -ALPS_DBL_MAX;
    fVarHardLBInd[k] = k;
    fVarHardUBInd[k] = k;
  }
  int numOldCons = 0;
  int maxOld = model->getOldConstraintsSize();
  BcpsObject ** oldConstraints = new BcpsObject* [maxOld];
  //--------------------------------------------------
  // Travel back to a full node, then collect diff (add/rem col/row,
  // hard/soft col/row bounds) from the node full to this node.
  //--------------------------------------------------------
  DcoNodeDesc* pathDesc = NULL;
  AlpsTreeNode *parent = parent_;
  std::vector<AlpsTreeNode*> leafToRootPath;
  leafToRootPath.push_back(this);
  while(parent) {
    leafToRootPath.push_back(parent);
    if (parent->getExplicit()) {
      // Reach an explicit node, then stop.
      break;
    }
    else {
      parent = parent->getParent();
    }
  }
  //------------------------------------------------------
  // Travel back from this node to the explicit node to
  // collect full description.
  //------------------------------------------------------
  for(i = static_cast<int> (leafToRootPath.size() - 1); i > -1; --i) {
    pathDesc = dynamic_cast<DcoNodeDesc*>((leafToRootPath.at(i))->
					  getDesc());
    //--------------------------------------
    // Full variable hard bounds.
    //--------------------------------------
    numModify = pathDesc->getVars()->lbHard.numModify;
    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->lbHard.posModify[k];
      value = pathDesc->getVars()->lbHard.entries[k];
      fVarHardLB[index] = value;
    }
    numModify = pathDesc->getVars()->ubHard.numModify;
    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->ubHard.posModify[k];
      value = pathDesc->getVars()->ubHard.entries[k];
      fVarHardUB[index] = value;
    }
    //--------------------------------------
    // Full variable soft bounds.
    //--------------------------------------
    numModify = pathDesc->getVars()->lbSoft.numModify;
    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->lbSoft.posModify[k];
      value = pathDesc->getVars()->lbSoft.entries[k];
      fVarSoftLB[index] = value;
    }
    numModify = pathDesc->getVars()->ubSoft.numModify;
    for (k = 0; k < numModify; ++k) {
      index = pathDesc->getVars()->ubSoft.posModify[k];
      value = pathDesc->getVars()->ubSoft.entries[k];
      fVarSoftUB[index] = value;
    }
    //----------------------------------------------
    // Collect all generated constraints, then remove deleted.
    //----------------------------------------------
    tempInt = pathDesc->getCons()->numAdd;
    for (k = 0; k < tempInt; ++k) {
      aCon = dynamic_cast<DcoConstraint *>
	(pathDesc->getCons()->objects[k]);
      assert(aCon);
      assert(aCon->getSize() > 0);
      assert(aCon->getSize() < 100000);
      oldConstraints[numOldCons++] = aCon;
      if (numOldCons >= maxOld) {
	// Need resize
	maxOld *= 2;
	BcpsObject **tempCons = new BcpsObject* [maxOld];
	memcpy(tempCons,
	       oldConstraints,
	       numOldCons * sizeof(BcpsObject *));
	delete [] oldConstraints;
	oldConstraints = tempCons;
	tempCons = NULL;
      }
    }
    //----------------------------------------------
    // Remove those deleted.
    // NOTE: oldConstraints stores all previously
    // generated active constraints at parent.
    //----------------------------------------------
    tempInt = pathDesc->getCons()->numRemove;
    if (tempInt > 0) {
      int tempPos;
      int *tempMark = new int [numOldCons];
      CoinZeroN(tempMark, numOldCons);
      for (k = 0; k < tempInt; ++k) {
	tempPos = pathDesc->getCons()->posRemove[k];
	tempMark[tempPos] = 1;
      }
      tempInt = 0;
      for (k = 0; k < numOldCons; ++k) {
	if (tempMark[k] != 1) {
	  // Survived.
	  oldConstraints[tempInt++] = oldConstraints[k];
	}
      }
      if (tempInt + pathDesc->getCons()->numRemove != numOldCons) {
	std::cout << "INSTALL: tempInt=" << tempInt
		  <<", numRemove="<<pathDesc->getCons()->numRemove
		  << ", numOldCons=" << numOldCons << std::endl;

	assert(0);
      }
      // Update number of old non-core constraints.
      numOldCons = tempInt;
      delete [] tempMark;
    }
  } // EOF for (path)
    //------------------------------------------
    // Record hard variable bounds. FULL set.
    //------------------------------------------
  desc->assignVarHardBound(numCols,
			   fVarHardLBInd,
			   fVarHardLB,
			   numCols,
			   fVarHardUBInd,
			   fVarHardUB);
  //------------------------------------------
  // Recode soft variable bound. Modified.
  //------------------------------------------
  for (k = 0; k < numCols; ++k) {
    if (fVarSoftLB[k] < ALPS_BND_MAX) {
      fVarSoftLBInd[numSoftVarLowers] = k;
      fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
    }
    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
      fVarSoftUBInd[numSoftVarUppers] = k;
      fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
    }
  }
  // Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
  desc->assignVarSoftBound(numSoftVarLowers,
			   fVarSoftLBInd,
			   fVarSoftLB,
			   numSoftVarUppers,
			   fVarSoftUBInd,
			   fVarSoftUB);
  //------------------------------------------
  // Recode added constraints.
  //------------------------------------------
  // First make a hard copy.
  for (k = 0; k < numOldCons; ++k) {
    aCon = dynamic_cast<DcoConstraint *>(oldConstraints[k]);
    assert(aCon);
    DcoConstraint *newCon = new DcoConstraint(*aCon);
    oldConstraints[k] = newCon;
  }
  // Add will first delete, then add. It is safe to use in parallel.
  desc->setAddedConstraints(numOldCons, oldConstraints);
  //------------------------------------------
  // Recode deleted constraints.
  //------------------------------------------
  desc->delConstraints(0, NULL);
  //--------------------------------------------------
  // Clear path vector.
  //--------------------------------------------------
  leafToRootPath.clear();
  assert(leafToRootPath.size() == 0);
}

//#############################################################################

// Not defined yet.
void
DcoTreeNode::convertToRelative()
{
  if(explicit_) {


  }
}

//#############################################################################

bool
DcoTreeNode::parallel(DcoModel *model,
		      BcpsConstraintPool &conPool,
		      int lastNew,
		      DcoConstraint *aCon)
{
  bool parallel = false;
  int k;
  double threshold = 0.999;

  //------------------------------------------------------
  // Compare with old cuts
  //------------------------------------------------------
#if 0  // couse error for VRP
  int numOldCons = model->getNumOldConstraints();
  for (k = 0; k < numOldCons; ++k) {
    DcoConstraint *aCon = model->oldConstraints()[k];
    assert(aCon);
    parallel = DcoParallelCutCon(rowCut,
				 aCon,
				 threshold);
    if (parallel) return parallel;
  }
#endif

  //------------------------------------------------------
  // Compare with new cuts
  //------------------------------------------------------

  for (k = 0; k < lastNew; ++k) {
    DcoConstraint *thisCon =
      dynamic_cast<DcoConstraint *>(conPool.getConstraint(k));
    parallel = DcoParallelConCon(aCon,
				 thisCon,
				 threshold);
    if (parallel) return parallel;
  }

  return parallel;
}

//#############################################################################

double
DcoTreeNode::estimateSolution(DcoModel *model,
			      const double *lpSolution,
			      double lpObjValue) const
{
  // lpObjective + sum_i{downCost_i*f_i + upCost_i*(1-f_i)}
  int k, col;
  int numInts= model->getNumIntObjects();

  double x, f, downC, upC, estimate = lpObjValue;

  DcoObjectInt *obj = NULL;

  for (k = 0; k < numInts; ++k) {
    obj = dynamic_cast<DcoObjectInt *>(model->objects(k));
    col = obj->columnIndex();
    x = lpSolution[col];
    f = CoinMax(0.0, x - floor(x));
    if (f > model->integerTol_) {
      downC = obj->pseudocost().getDownCost();
      upC = obj->pseudocost().getUpCost();
      estimate += (downC * f + upC * (1-f));
    }
  }
  return estimate;
}

//#############################################################################

int
DcoTreeNode::callHeuristics(DcoModel *model, bool onlyBeforeRoot)
{
  int status = 0;

  if (model->heurStrategy_ == DcoHeurStrategyNone) {
    return status;
  }

  int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);

  int foundSolution = false;
  int numCols = model->solver()->getNumCols();

  double heurObjValue = getKnowledgeBroker()->getIncumbentValue();
  double *heurSolution = new double [numCols];

  DcoSolution *bSol = NULL;

  for (int k = 0; k < model->numHeuristics(); ++k) {
    int heurStrategy = model->heuristics(k)->strategy();
    //std::cout << " call heur " << k << "; strategy = "
    //        << heurStrategy << std::endl;

    if (heurStrategy != DcoHeurStrategyNone) {
      if (onlyBeforeRoot) {
	// heuristics that can only be used before root.
	if (heurStrategy != DcoHeurStrategyBeforeRoot) {
	  continue;
	}
      }
      else {
	// regular heuristics
	if (heurStrategy == DcoHeurStrategyBeforeRoot) {
	  continue;
	}
      }

      getKnowledgeBroker()->tempTimer().start();
      foundSolution = false;
      foundSolution =
	model->heuristics(k)->searchSolution(heurObjValue,
					     heurSolution);
      getKnowledgeBroker()->tempTimer().stop();

      model->heuristics(k)->
	addTime(getKnowledgeBroker()->tempTimer().getCpuTime());

      if (foundSolution) {
	// Check if solution from heuristic is feasible.
	bSol = model->feasibleSolutionHeur(heurSolution);
      }

      if (bSol) {
	model->heuristics(k)->addNumSolutions(1);
	int noSols = model->heuristics(k)->noSolCalls();
	model->heuristics(k)->addNoSolCalls(-noSols);
	// Store the newly found blis solution.
	model->storeSolution(DcoSolutionTypeHeuristic, bSol);
	if (onlyBeforeRoot) {
	  status = 1;
	}
	else if (quality_ >  model->getCutoff()) {
	  setStatus(AlpsNodeStatusFathomed);
	  status = 2;
	  goto TERM_HEUR;
	}
	else {
	  status = 1;
	}
	if (heurStrategy == DcoHeurStrategyBeforeRoot &&
	    msgLevel > 200) {
	  model->blisMessageHandler()->message(DISCO_HEUR_BEFORE_ROOT,
					       model->blisMessages())
	    << (model->heuristics(k)->name())
	    << bSol->getQuality()
	    << CoinMessageEol;
	}
      }
      else {
	model->heuristics(k)->addNoSolCalls(1);
      }
    }
  }

 TERM_HEUR:

  if (heurSolution) delete [] heurSolution;

  return status;
}

//#############################################################################

bool DcoTreeNode::fractional_vars_exist() const {
  DcoModel const * model = dynamic_cast<DcoModel*>(desc_->getModel());
  double const * sol = model->getLpSolution();
  int numIntegers = model->getNumIntObjects();
  const int *intIndices = model->getIntColIndices();
  for (int i=0; i<numIntegers; ++i) {
    // check whether variable i is fractional?
  }
  return true;
}

#if defined(__OA__)
void DcoTreeNode::integerFix(DcoModel * model, OsiConicSolverInterface * ipm_solver) const {
  // build core problem and fix integer variables.
  OsiSolverInterface * si = model->solver();
  int num_cols = si->getNumCols();
  CoinPackedMatrix const * matrix = si->getMatrixByCol();
  double const * rowlb = si->getRowLower();
  double const * rowub = si->getRowUpper();
  double const * collb = si->getColLower();
  double const * colub = si->getColUpper();
  double const * obj = si->getObjCoefficients();
  ipm_solver->loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
  // add cones to the problem
  int num_cones = model->getNumCoreCones();
  OsiLorentzConeType const * cone_type = model->getConeTypes();
  int const * cone_size = model->getConeSizes();
  int const * const * members = model->getConeMembers();
  for (int i=0; i<num_cones; ++i) {
    ipm_solver->addConicConstraint(cone_type[i], cone_size[i], members[i]);
  }
  // fix integer variables
  double const * sol = si->getColSolution();
  int const * ind = model->getIntColIndices();
  for (int i=0; i<model->getNumIntObjects(); ++i) {
    int index = ind[i];
    double value = sol[index];
    ipm_solver->setColBounds(index, value, value);
  }
  ipm_solver->initialSolve();
  if (ipm_solver->isProvenOptimal()) {
    si->setColSolution(ipm_solver->getColSolution());
  }
}
#endif
