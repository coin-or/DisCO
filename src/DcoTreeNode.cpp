#include "DcoTreeNode.hpp"

#include <OsiRowCut.hpp>

#include "DcoNodeDesc.hpp"
#include "DcoMessage.hpp"
#include "DcoLinearConstraint.hpp"

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
  DcoNodeDesc * node_desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  DcoModel * model = dynamic_cast<DcoModel*>(node_desc->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
  throw std::exception();
}

void DcoTreeNode::convertToRelative() {
  DcoNodeDesc * node_desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  DcoModel * model = dynamic_cast<DcoModel*>(node_desc->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
  throw std::exception();
}

int DcoTreeNode::generateConstraints(BcpsModel * model,
				     BcpsConstraintPool * conPool) {
  DcoModel * disco_model = dynamic_cast<DcoModel*>(model);
  CoinMessageHandler * message_handler = disco_model->dcoMessageHandler_;
  CoinMessages * messages = disco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
  throw std::exception();
  return 0;
}

int DcoTreeNode::generateVariables(BcpsModel * model,
				   BcpsVariablePool * varPool) {
  DcoModel * disco_model = dynamic_cast<DcoModel*>(model);
  CoinMessageHandler * message_handler = disco_model->dcoMessageHandler_;
  CoinMessages * messages = disco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
  throw std::exception();
  return 0;
}

int DcoTreeNode::chooseBranchingObject(BcpsModel * model) {
  DcoModel * disco_model = dynamic_cast<DcoModel*>(model);
  CoinMessageHandler * message_handler = disco_model->dcoMessageHandler_;
  CoinMessages * messages = disco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
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
  AlpsNodeStatus status = getStatus();
  DcoModel * model = getModel();
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  if (status==AlpsNodeStatusPregnant) {
    // branch
    std::cout << "Gotta branch!" << std::endl;
  }
  else if (status==AlpsNodeStatusCandidate or
	   status==AlpsNodeStatusEvaluated) {
    boundingLoop(isRoot, rampUp);
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

int DcoTreeNode::boundingLoop(bool isRoot, bool rampUp) {
  AlpsNodeStatus status = getStatus();
  DcoNodeDesc * desc = getDesc();
  DcoModel * model = getModel();
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  bool keepBounding = true;
  bool fathomed = false;
  bool genConstraints = false;
  bool genVariables = false;
  BcpsConstraintPool * constraintPool = NULL;
  BcpsVariablePool * variablePool = NULL;
  installSubProblem(model);
  while (keepBounding) {
    keepBounding = false;
    // solve subproblem corresponds to this node
    DcoSubproblemStatus subproblem_status =
      static_cast<DcoSubproblemStatus> (bound(model));
    decide(subproblem_status, genConstraints, genVariables);
    if (getStatus()==AlpsNodeStatusFathomed and !keepBounding) {
      break;
      // node is fathomed, nothing to do.
    }
    else if (keepBounding and genConstraints) {
      generateConstraints(model, constraintPool);
    }
    else if (keepBounding and genVariables) {
      generateVariables(model, variablePool);
    }
  }
  return AlpsReturnStatusOk;
}

/** Bounding procedure to estimate quality of this node. */
// todo(aykut) Why do we need model as input?
// desc_ is a member of AlpsTreeNode, and model_ is a member of AlpsNodeDesc.
int DcoTreeNode::bound(BcpsModel * bcps_model) {
  DcoSubproblemStatus subproblem_status;
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
    subproblem_status = DcoSubproblemStatusAbandoned;
  }
  else if (model->solver()->isProvenOptimal()) {
    // todo(aykut) if obj val is greater than 1e+30 we consider problem
    // infeasible. This is due to our cut generation when the problem is
    // infeasible. We add a high lower bound (1e+30) to the objective
    // function. This can be improved.
    if (model->solver()->getObjValue()>=1e+30) {
      subproblem_status = DcoSubproblemStatusPrimalInfeasible;
    }
    else {
      subproblem_status = DcoSubproblemStatusOptimal;
      double objValue = model->solver()->getObjValue() *
	model->solver()->getObjSense();
      // Update quality of this nodes.
      quality_ = objValue;
    }
  }
  else if (model->solver()->isProvenPrimalInfeasible()) {
    subproblem_status = DcoSubproblemStatusPrimalInfeasible;
  }
  else if (model->solver()->isProvenDualInfeasible()) {
    subproblem_status = DcoSubproblemStatusDualInfeasible;
  }
  else if (model->solver()->isPrimalObjectiveLimitReached()) {
    subproblem_status = DcoSubproblemStatusPrimalObjLim;
  }
  else if (model->solver()->isDualObjectiveLimitReached()) {
    subproblem_status = DcoSubproblemStatusDualObjLim;
  }
  else if (model->solver()->isIterationLimitReached()) {
    subproblem_status = DcoSubproblemStatusIterLim;
  }
  else {
    message_handler->message(DISCO_SOLVER_UNKNOWN_STATUS, *messages)
      << CoinMessageEol;
  }
  return subproblem_status;
}

int DcoTreeNode::installSubProblem(BcpsModel * bcps_model) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  DcoModel * model = dynamic_cast<DcoModel*>(bcps_model);
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  DcoNodeDesc * desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  // number of columns and rows
  int numCoreCols = model->getNumCoreVariables();
  int numCoreLinearRows = model->getNumCoreLinearConstraints();
  int numSolverCols = model->solver()->getNumCols();
  // solver rows are all linear, for both Osi and OsiConic.
  int numSolverRows = model->solver()->getNumRows();
  // set col and row bounds
  double * colLB = model->colLB();
  double * colUB = model->colUB();
  double * rowLB = model->colLB();
  double * rowUB = model->colUB();
  CoinFillN(colLB, numCoreCols, -ALPS_DBL_MAX);
  CoinFillN(colUB, numCoreCols, ALPS_DBL_MAX);
  CoinFillN(rowLB, numCoreLinearRows, -ALPS_DBL_MAX);
  CoinFillN(rowUB, numCoreLinearRows, ALPS_DBL_MAX);
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
  int numDelRows = numSolverRows - numCoreLinearRows;
  if (numDelRows > 0) {
    int * indices = new int[numDelRows];
    if (indices==NULL) {
      message_handler->message(DISCO_OUT_OF_MEMORY, *messages)
	<< __FILE__ << __LINE__ << CoinMessageEol;
      throw CoinError("Out of memory", "installSubProblem", "DcoTreeNode");
    }
    for (int i=0; i<numDelRows; ++i) {
      indices[i] = numCoreLinearRows + i;
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
  std::vector<DcoConstraint*> old_cons;
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
    for (int k=0; k<tempInt; ++k) {
      DcoConstraint * aCon = dynamic_cast<DcoConstraint *>
	(currDesc->getCons()->objects[k]);
      old_cons.push_back(aCon);
    }
    //----------------------------------------------
    // Remove those deleted.
    // NOTE: old_cons stores all previously
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
	  old_cons[tempInt++] = old_cons[k];
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
  if (numOldRows > 0) {
    const OsiRowCut ** oldOsiCuts = new const OsiRowCut * [numOldRows];
    for (int k=0; k<numOldRows; ++k) {
      OsiRowCut * acut = old_cons[k]->createOsiRowCut(model);
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
  old_cons.clear();
  return status;
}

/** This method must be invoked on a \c pregnant node (which has all the
    information needed to create the children) and should create the
    children's decriptions. The stati of the children
    can be any of the ones \c process() can return. */
std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
DcoTreeNode::branch() {
  // get description of this node
  DcoNodeDesc * node_desc = getDesc();
  // get model the node belongs
  DcoModel * model = dynamic_cast<DcoModel*>(node_desc->getModel());
  // get message handler and messages for loging
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  // check node status
  if (getStatus()!=AlpsNodeStatusPregnant) {
      message_handler->message(DISCO_NODE_UNEXPECTEDSTATUS, *messages)
	<< static_cast<int>(getStatus()) << CoinMessageEol;
  }  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > res;
  return res;
}

DcoNodeDesc * DcoTreeNode::getDesc() const {
  return dynamic_cast<DcoNodeDesc*>(AlpsTreeNode::getDesc());
}

DcoModel * DcoTreeNode::getModel() const {
  return getDesc()->getModel();
}


void DcoTreeNode::processSetPregnant() {
  // get warm start basis from solver
  // todo(aykut) This does not help much if the underlying solver is an IPM
  // based solver.
  DcoModel * model = getModel();
  CoinWarmStartBasis * ws = dynamic_cast<CoinWarmStartBasis*>
    (model->solver()->getWarmStart());
  // store basis in the node desciption.
  getDesc()->setBasis(ws);
  // set status pregnant
  setStatus(AlpsNodeStatusPregnant);
}

void DcoTreeNode::decide(DcoSubproblemStatus subproblem_status,
			 bool & generateConstraints,
			 bool & generateVariables) {
  if (subproblem_status==DcoSubproblemStatusPrimalInfeasible or
      subproblem_status==DcoSubproblemStatusDualInfeasible) {
    // fathom
    setStatus(AlpsNodeStatusFathomed);
    return;
  }
  if (subproblem_status!=DcoSubproblemStatusOptimal) {
    std::cout << "Subproblem is not optimal. Donnow what to do!"  << std::endl;
    throw std::exception();
  }
  // subproblem is solved to optimality. Check feasibility of the solution.
  // iterate over objects and check their feasibility

  // what is in the subproblem and what is not in the subproblem?
}
