#include "DcoTreeNode.hpp"

#include <OsiRowCut.hpp>

#include "DcoNodeDesc.hpp"
#include "DcoMessage.hpp"

#include <vector>

DcoTreeNode::DcoTreeNode() {
}

DcoTreeNode::DcoTreeNode(AlpsNodeDesc *& desc) {
  desc_ = desc;
  desc = NULL;
}

DcoTreeNode::~DcoTreeNode() {
}

// create tree node from given description
AlpsTreeNode * DcoTreeNode::createNewTreeNode(AlpsNodeDesc *& desc) const {
  DcoNodeDesc * dco_node = dynamic_cast<DcoNodeDesc*>(desc);
  int branched_index = dco_node->getBranchedInd();
  double branched_value = dco_node->getBranchedVal();
  double frac = branched_value - floor(branched_value);
  DcoModel * model = dynamic_cast<DcoModel*>(desc->getModel());
  if (frac > 0.0) {
    model->dcoMessageHandler_->message(DISCO_NODE_BRANCHONINT,
				       *(model->dcoMessages_))
      << branched_index << CoinMessageEol;
  }
  // Create a new tree node
  DcoTreeNode * node = new DcoTreeNode(desc);
  desc = NULL;
  return node;
}

void DcoTreeNode::convertToExplicit() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

void DcoTreeNode::convertToRelative() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

int DcoTreeNode::generateConstraints(BcpsModel * model,
				     BcpsConstraintPool * conPool) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::generateVariables(BcpsModel * model,
				   BcpsVariablePool * varPool) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::chooseBranchingObject(BcpsModel * model) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::handleBoundingStatus(int status, bool & keepOn, bool & fathomed) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

/**
   This function will process the node.
   <ul>
   <il> If it is candidate or evaluated, decide to either generate cut or branch.
   <il> If it is pregnant, branch.
   <il> If it is branched, fathomed or discarded raise error.
   </ul>

 */
int DcoTreeNode::process(bool isRoot, bool rampUp) {
  DcoNodeDesc * desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  DcoModel * model = dynamic_cast<DcoModel*>(desc->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  AlpsNodeStatus status = getStatus();
  installSubProblem(model);
  if (status==AlpsNodeStatusCandidate or
      status==AlpsNodeStatusEvaluated) {
    // generate cuts or branch
    bound(model);
  }
  else if (status==AlpsNodeStatusPregnant) {
    // branch
  }
  else if (status==AlpsNodeStatusBranched or
	   status==AlpsNodeStatusFathomed or
	   status==AlpsNodeStatusDiscarded) {
    // this should not happen
    message_handler->message(DISCO_NODE_UNEXPECTEDSTATUS, *messages)
      << static_cast<int>(status) << CoinMessageEol;

  }
  return AlpsReturnStatusOk;
}

/** Bounding procedure to estimate quality of this node. */
// todo(aykut) Why do we need model as input?
// desc_ is a member of AlpsTreeNode, and model_ is a member of AlpsNodeDesc.
int DcoTreeNode::bound(BcpsModel * bcps_model) {
  DcoLpStatus lp_status;
  DcoModel * model = dynamic_cast<DcoModel*>(bcps_model);
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  AlpsNodeStatus node_status = getStatus();
  if (node_status==AlpsNodeStatusPregnant or
      node_status==AlpsNodeStatusBranched or
      node_status==AlpsNodeStatusFathomed or
      node_status==AlpsNodeStatusDiscarded) {
    message_handler->message(DISCO_NODE_UNEXPECTEDSTATUS, *messages)
      << static_cast<int>(node_status) << CoinMessageEol;
  }
  // solve problem loaded to the solver
  model->solver()->resolve();
  if (model->solver()->isAbandoned()) {
    lp_status = DcoLpStatusAbandoned;
  }
  else if (model->solver()->isProvenOptimal()) {
    // todo(aykut) if obj val is greater than 1e+30 we consider problem
    // infeasible. This is due to our cut generation when the problem is
    // infeasible. We add a high lower bound (1e+30) to the objective
    // function. This can be improved.
    if (model->solver()->getObjValue()>=1e+30) {
      lp_status = DcoLpStatusPrimalInfeasible;
    }
    else {
      lp_status = DcoLpStatusOptimal;
      double objValue = model->solver()->getObjValue() *
	model->solver()->getObjSense();
      // Update quality of this nodes.
      quality_ = objValue;
    }
  }
  else if (model->solver()->isProvenPrimalInfeasible()) {
    lp_status = DcoLpStatusPrimalInfeasible;
  }
  else if (model->solver()->isProvenDualInfeasible()) {
    lp_status = DcoLpStatusDualInfeasible;
  }
  else if (model->solver()->isPrimalObjectiveLimitReached()) {
    lp_status = DcoLpStatusPrimalObjLim;
  }
  else if (model->solver()->isDualObjectiveLimitReached()) {
    lp_status = DcoLpStatusDualObjLim;
  }
  else if (model->solver()->isIterationLimitReached()) {
    lp_status = DcoLpStatusIterLim;
  }
  else {
    message_handler->message(DISCO_SOLVER_UNKNOWN_STATUS, *messages)
      << CoinMessageEol;
  }
  return lp_status;
}

/**
    install subproblem to the solver. This involes changing corresponding
    variable bounds.
 */
int DcoTreeNode::installSubProblem(BcpsModel * bcps_model) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  DcoModel * model = dynamic_cast<DcoModel*>(bcps_model);
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  DcoNodeDesc * desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  // number of columns and rows
  int numCoreCols = model->getNumCoreVariables();
  int numCoreRows = model->getNumCoreLinearConstraints();
  int numCols = model->solver()->getNumCols();
  int numRows = model->solver()->getNumRows();
  // set col and row bounds
  double * colLB = model->colLB();
  double * colUB = model->colUB();
  double * rowLB = model->colLB();
  double * rowUB = model->colUB();
  CoinFillN(colLB, numCoreCols, -ALPS_DBL_MAX);
  CoinFillN(colUB, numCoreCols, ALPS_DBL_MAX);
  CoinFillN(rowLB, numCoreRows, -ALPS_DBL_MAX);
  CoinFillN(rowUB, numCoreRows, ALPS_DBL_MAX);
  //======================================================
  // Restore subproblem:
  //  1. Remove noncore columns and rows
  //  2. Travel back to root and correct differencing to
  //     full column/row bounds into model->startXXX
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
    int * indices = new int [numDelRows];
    if (indices==NULL) {
      message_handler->message(DISCO_OUT_OF_MEMORY, *messages)
	<< __FILE__ << __LINE__ << CoinMessageEol;
      throw CoinError("Out of memory", "installSubProblem", "DcoTreeNode");
    }
    for (int i=0; i<numDelRows; ++i) {
      indices[i] = numCoreRows + i;
    }
    model->solver()->deleteRows(numDelRows, indices);
    delete[] indices;
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
  /* First push this node since it has branching hard bounds.
     NOTE: during rampup, this desc has full description when branch(). */
  std::vector<AlpsTreeNode*> leafToRootPath;
  leafToRootPath.push_back(this);
  if (knowledgeBroker_->getPhase() != AlpsPhaseRampup) {
    AlpsTreeNode * parent = parent_;
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
  int numOldRows = 0;
  for(int i = static_cast<int> (leafToRootPath.size() - 1); i > -1; --i) {
    //--------------------------------------------------
    // NOTE: As away from explicit node, bounds become
    //       tighter and tighter.
    //--------------------------------------------------
    // node description of node i
    DcoNodeDesc * currDesc =
      dynamic_cast<DcoNodeDesc*>((leafToRootPath.at(i))->getDesc());
    //--------------------------------------------------
    // Adjust bounds according to hard var lb/ub.
    // If rampup or explicit, collect hard bounds so far.
    //--------------------------------------------------
    int index;
    double value;
    int numModify;
    numModify = currDesc->getVars()->lbHard.numModify;
    for (int k=0; k<numModify; ++k) {
      index = currDesc->getVars()->lbHard.posModify[k];
      value = currDesc->getVars()->lbHard.entries[k];
      // Hard bounds do NOT change according to soft bounds, so
      // here need CoinMax.
      colLB[index] = CoinMax(colLB[index], value);
    }
    numModify = currDesc->getVars()->ubHard.numModify;
    for (int k=0; k<numModify; ++k) {
      int index = currDesc->getVars()->ubHard.posModify[k];
      double value = currDesc->getVars()->ubHard.entries[k];
      colUB[index] = CoinMin(colUB[index], value);
    }
    //--------------------------------------------------
    // Adjust bounds according to soft var lb/ub.
    // If rampup or explicit, collect soft bounds so far.
    //--------------------------------------------------
    numModify = currDesc->getVars()->lbSoft.numModify;
    for (int k=0; k<numModify; ++k) {
      index = currDesc->getVars()->lbSoft.posModify[k];
      value = currDesc->getVars()->lbSoft.entries[k];
      colLB[index] = CoinMax(colLB[index], value);
    }
    numModify = currDesc->getVars()->ubSoft.numModify;
    for (int k=0; k<numModify; ++k) {
      index = currDesc->getVars()->ubSoft.posModify[k];
      value = currDesc->getVars()->ubSoft.entries[k];
      colUB[index] = CoinMin(colUB[index], value);
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
    int tempInt = currDesc->getCons()->numAdd;
    int maxOld = model->getOldConstraintsSize();
    for (int k=0; k<tempInt; ++k) {
      DcoConstraint * aCon = dynamic_cast<DcoConstraint *>
	(currDesc->getCons()->objects[k]);
      (model->oldConstraints())[numOldRows++] = aCon;
      if (numOldRows >= maxOld) {
	// Need resize
	maxOld *= 2;
	DcoConstraint **tempCons = new DcoConstraint* [maxOld];
	std::copy(model->oldConstraints(),
		  model->oldConstraints()+numOldRows, tempCons);
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
    tempInt = currDesc->getCons()->numRemove;
    if (tempInt > 0) {
      int tempPos;
      int * tempMark = new int [numOldRows];
      CoinZeroN(tempMark, numOldRows);
      for (int k=0; k<tempInt; ++k) {
	tempPos = currDesc->getCons()->posRemove[k];
	tempMark[tempPos] = 1;
      }
      tempInt = 0;
      for (int k=0; k<numOldRows; ++k) {
	if (tempMark[k] != 1) {
	  // Survived.
	  (model->oldConstraints())[tempInt++]=
	    (model->oldConstraints())[k];
	}
      }
      if (tempInt + currDesc->getCons()->numRemove != numOldRows) {
	std::cout << "INSTALL: tempInt=" << tempInt
		  <<", numRemove="<<currDesc->getCons()->numRemove
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
  model->solver()->setColLower(colLB);
  model->solver()->setColUpper(colUB);
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
    for (int k=0; k<numOldRows; ++k) {
      OsiRowCut * acut = (model->oldConstraints()[k])->createOsiRowCut(model);
      oldOsiCuts[k] = acut;
    }
    model->solver()->applyRowCuts(numOldRows, oldOsiCuts);
    for (int k=0; k<numOldRows; ++k) {
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

/** This method must be invoked on a \c pregnant node (which has all the
    information needed to create the children) and should create the
    children's decriptions. The stati of the children
    can be any of the ones \c process() can return. */
std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
DcoTreeNode::branch() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > res;
  return res;
}
