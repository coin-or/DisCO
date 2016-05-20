// CoinUtils headers
#include "CoinUtility.hpp"

// Osi headers
#include <OsiRowCut.hpp>

// Disco headers
#include "DcoTreeNode.hpp"
#include "DcoNodeDesc.hpp"
#include "DcoMessage.hpp"
#include "DcoLinearConstraint.hpp"
#include "DcoConicConstraint.hpp"
#include "Dco.hpp"
#include "DcoConGenerator.hpp"
#include "DcoSolution.hpp"
#include "DcoBranchObject.hpp"

// STL headers
#include <vector>

// define structs that will be used in creating nodes.
struct SparseVector {
  int * ind;
  double * val;
};

struct Bound {
  SparseVector lower;
  SparseVector upper;
};


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
  DcoModel * model = dynamic_cast<DcoModel*>(desc->getModel());

  double dist_to_floor = branched_value - floor(branched_value);
  double dist_to_ceil = ceil(branched_value) - branched_value;
  double frac = (dist_to_floor>dist_to_ceil) ? dist_to_ceil : dist_to_floor;
  double tol = model->dcoPar()->entry(DcoParams::integerTol);

  if (frac < tol) {
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

/// Generate constraints for the problem.
int DcoTreeNode::generateConstraints(BcpsModel * model,
                                     BcpsConstraintPool * conPool) {
  // todo(aykut) why do we need this status?
  DcoModel * disco_model = dynamic_cast<DcoModel*>(model);
  CoinMessageHandler * message_handler = disco_model->dcoMessageHandler_;
  CoinMessages * messages = disco_model->dcoMessages_;

  // number of constraint generators in model
  int num_cg = disco_model->numCutGenerators();
  for (int i=0; i<num_cg; ++i) {
    bool do_use = false;
    DcoConGenerator * cg = disco_model->conGenerators(i);
    // decide whether we should use this cut generator with respect to the
    // specified cut strategy
    decide_using_cg(do_use, cg);
    if (!do_use) {
      // jump to the next geenrator if we will not use this one.
      continue;
    }
    int pre_num_cons = conPool->getNumConstraints();
    double start_time = CoinCpuTime();
    // Call constraint generator
    bool must_resolve = cg->generateConstraints(*conPool);
    double cut_time = CoinCpuTime() - start_time;
    // Statistics
    cg->addTime(cut_time);
    cg->addCalls(1);
    int num_cons_generated = conPool->getNumConstraints() - pre_num_cons;
    if (num_cons_generated == 0) {
      cg->addNoConsCalls(1);
    }
    else {
      cg->addNumConsGenerated(num_cons_generated);
    }
    // if (must_resolve) {
    //   // todo(aykut) conic cuts may need to resolve?
    //   // TODO: Only probing will return ture.
    //   if (static_cast<DcoLpStatus> (bound(model)) != DcoLpStatusOptimal) {
    //     // probing resolve failed.
    //     std::cerr << "Probing resolve failed."  << std::endl;
    //     break;
    //   }
    // }
    //------------------------------------------------
    // Modify control.
    // NOTE: only modify if user choose automatic.
    //------------------------------------------------
  }
  // return value will make sense when DcoTreeNode::process is implemented in Bcps level.
  return 0;
}


void DcoTreeNode::decide_using_cg(bool & do_use,
                                  DcoConGenerator const * cg) const {
  DcoModel * model = getModel();
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  DcoCutStrategy strategy = cg->strategy();
  // todo(aykut) make it a parameter, Only autmatic stategy has depth limit.
  int maxConstraintDepth = 20;

  do_use = false;

  switch (strategy) {
  case DcoCutStrategyNone:
    do_use = false;
    break;
  case DcoCutStrategyRoot:
    if (depth_==0) {
      do_use = true;
    }
    break;
  case DcoCutStrategyAuto:
    if (depth_ < maxConstraintDepth) {
      if (!diving_ || depth_==0) {
        do_use = true;
      }

    }
    break;
  case DcoCutStrategyPeriodic:
    if (index_ % cg->frequency() == 0) {
      do_use = true;
    }
    break;
  default:
    message_handler->message(DISCO_UNKNOWN_CUTSTRATEGY, *messages)
      << strategy << CoinMessageEol;
  }
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

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "Processing node ";
  debug_msg << this;
  debug_msg << " index ";
  debug_msg << getIndex();
  debug_msg << " parent ";
  debug_msg << getParent();
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                           'G', DISCO_DLOG_PROCESS)
    << CoinMessageEol;
  // end of debug stuff

  if (status==AlpsNodeStatusCandidate or
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
  bool do_branch = false;
  bool genConstraints = false;
  bool genVariables = false;
  BcpsConstraintPool * constraintPool = new BcpsConstraintPool();
  BcpsVariablePool * variablePool = new BcpsVariablePool();
  installSubProblem(model);
  while (keepBounding) {
    keepBounding = false;
    // solve subproblem corresponds to this node
    DcoSubproblemStatus subproblem_status =
      static_cast<DcoSubproblemStatus> (bound(model));
    // decide what to do
    branchConstrainOrPrice(subproblem_status, keepBounding, do_branch,
                           genConstraints,
                           genVariables);
    if (getStatus()==AlpsNodeStatusFathomed) {
      // node is fathomed, nothing to do.
      break;
    }
    else if (keepBounding and genConstraints) {
      generateConstraints(model, constraintPool);
      // add constraints to the model
      applyConstraints(*constraintPool);
      // clear constraint pool
      constraintPool->freeGuts();
    }
    else if (keepBounding and genVariables) {
      generateVariables(model, variablePool);
      // add variables to the model
    }
    else if (keepBounding==false and do_branch==false) {
      // put node back into the list.
      // this means do not change the node status and end processing the node.
    }
    else if (keepBounding==false and do_branch) {
      // branch
      BcpsBranchStrategy * branchStrategy = model->branchStrategy();
      // todo(aykut) following should be a parameter
      // Maximum number of resolve during branching.
      int numBranchResolve = 10;
      // todo(aykut) why ub should be an input?
      branchStrategy->createCandBranchObjects(this);
      // prepare this node for branching, bookkeeping for differencing.
      // call pregnant setting routine
      processSetPregnant();
    }
    else {
      message_handler->message(9998, "Dco", "This should not happen. "
                               " branchConstrainOrPrice() is buggy.", 'E', 0)
        << CoinMessageEol;
    }

  }
  delete constraintPool;
  delete variablePool;
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
      // Update quality of this node
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
  // load the problem in this node to the solver
  //
  //
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
  AlpsReturnStatus status = AlpsReturnStatusOk;
  DcoModel * model = dynamic_cast<DcoModel*>(bcps_model);
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  DcoNodeDesc * desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  // get number of columns and rows
  int numCoreCols = model->getNumCoreVariables();
  int numCoreLinearRows = model->getNumCoreLinearConstraints();
  int numCoreConicRows = model->getNumCoreConicConstraints();
  // get number of columns and rows stored in the solver
  int numSolverCols = model->solver()->getNumCols();
  // solver rows are all linear, for both Osi and OsiConic.
  int numSolverRows = model->solver()->getNumRows();

  // 1. Remove noncore columns and rows
  // 1.1 Remove non-core rows from solver, i.e. cuts
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
  // 1.1 Remove non-core columns from solver
  // End of 1.

  //  2. Travel back to root and correct differencing to full column/row bounds
  // into model->startXXX
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

  // column and row lower bounds
  // todo(aykut) why do we need to write them to model->colLB_ fields?
  // is it enough to have a local array?
  double * colLB = model->colLB();
  double * colUB = model->colUB();
  double * rowLB = model->rowLB();
  double * rowUB = model->rowUB();
  CoinFillN(colLB, numCoreCols, -ALPS_DBL_MAX);
  CoinFillN(colUB, numCoreCols, ALPS_DBL_MAX);
  CoinFillN(rowLB, numCoreLinearRows, -ALPS_DBL_MAX);
  CoinFillN(rowUB, numCoreLinearRows, ALPS_DBL_MAX);
  // generate path to root from this
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
  // End of 2

  //  3. Set col bounds
  //--------------------------------------------------------
  // Adjust column bounds in lp solver
  //--------------------------------------------------------
  model->solver()->setColLower(colLB);
  model->solver()->setColUpper(colUB);
  // End of 3


  //  4. Set row bounds (is this necessary?)

  //  5. Add contraints except cores
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
  old_cons.clear();
  //  End of 5


  //  6. Add variables except cores

  //  7. Set basis (should not need modify)
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

  // write mps file after problem is installed.
  // std::cout << "Writing problem corresponding to node " << this << std::endl;
  // std::stringstream ss;
  // ss << this;
  // model->solver()->writeMps(ss.str().c_str(), "mps", 0.0);
  // end of mps file writing

  return status;
  //  End of 7
}

/** This method must be invoked on a \c pregnant node (which has all the
    information needed to create the children) and should create the
    children's decriptions. The stati of the children
    can be any of the ones \c process() can return. */
std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
DcoTreeNode::branch() {
  // get model the node belongs
  DcoModel * model = getModel();
  // get message handler and messages for loging
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "Branching node ";
  debug_msg << this;
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                           'G', DISCO_DLOG_BRANCH)
    << CoinMessageEol;
  // end of debug stuff

  // check node status, this should be a pregnant node.
  if (getStatus()!=AlpsNodeStatusPregnant) {
      message_handler->message(DISCO_NODE_UNEXPECTEDSTATUS, *messages)
        << static_cast<int>(getStatus()) << CoinMessageEol;
  }

  //todo(aykut) update scores
  //BcpsBranchStrategy * branchStrategy = model->branchStrategy();


  // get Alps phase
  AlpsPhase phase = knowledgeBroker_->getPhase();

  // get branch object
  DcoBranchObject const * branch_object =
    dynamic_cast<DcoBranchObject const *>(branchObject());

  // get index and value of branch variable.
  //int branch_var = model->relaxedCols()[branch_object->getObjectIndex()];
  int branch_var = branch_object->index();
  double branch_value = branch_object->value();

  // compute child nodes' warm start basis
  CoinWarmStartBasis * child_ws;
  child_ws = (getDesc()->getBasis()==0) ? 0 :
    new CoinWarmStartBasis(*getDesc()->getBasis());

  // create new node descriptions
  DcoNodeDesc * down_node = new DcoNodeDesc(model);
  DcoNodeDesc * up_node = new DcoNodeDesc(model);
  if (phase == AlpsPhaseRampup) {
    // Store a full description in the child nodes

    // Down Node
    // == copy description of this node to the down node
    copyFullNode(down_node);
    // == update the branching variable hard bounds for the down node
    // todo(aykut) do we need lower bound for the down node?
    // down_node->vars_->lbHard.entries[branch_var] =
    //   model->getColLower()[branch_var];
    down_node->vars()->ubHard.entries[branch_var] =
      branch_object->ubDownBranch();

    // Up Node
    // == copy description of this node to up node
    copyFullNode(up_node);
    // == update the branching variable hard bounds for the up node
    // todo(aykut) do we need upper bound for the up node?
    up_node->vars()->lbHard.entries[branch_var] =
      branch_object->lbUpBranch();
    // up_node->vars_->ubHard.entries[branch_var] =
    //   model->getColUpper()[branch_var];
  }
  else {
    // Store node description relative to the parent.
    // We need to add a hard bound for the branching variable.

    double ub_down_branch = branch_object->ubDownBranch();
    double lb_up_branch = branch_object->lbUpBranch();
    double lb = model->colLB()[branch_var];
    double ub = model->colUB()[branch_var];
    down_node->setVarHardBound(1,
                             &branch_var,
                             &lb,
                             1,
                             &branch_var,
                             &ub_down_branch);
    up_node->setVarHardBound(1,
                             &branch_var,
                             &lb_up_branch,
                             1,
                             &branch_var,
                             &ub);
  }

  // Down Node
  // == set other relevant fields of down node
  down_node->setBranchedDir(DcoNodeBranchDirectionDown);
  down_node->setBranchedInd(branch_object->index());
  down_node->setBranchedVal(branch_value);
  // == set warm start basis for the down node.
  down_node->setBasis(child_ws);

  // Up Node
  // == set other relevant fields of up node
  up_node->setBranchedDir(DcoNodeBranchDirectionUp);
  up_node->setBranchedInd(branch_object->index());
  up_node->setBranchedVal(branch_value);
  // == set warm start basis for the up node.
  up_node->setBasis(child_ws);

  // Alps does this. We do not need to change the status here
  //status_ = AlpsNodeStatusBranched;

  // create return value and push the down and up nodes.
  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > res;
  res.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc*>(down_node),
                               AlpsNodeStatusCandidate,
                               model->solver()->getObjValue()));
  res.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc*>(up_node),
                               AlpsNodeStatusCandidate,
                               model->solver()->getObjValue()));
  return res;
}

// Copies node description of this node to given child node.
// New node is explicitly stored in the memory (no differencing).
void DcoTreeNode::copyFullNode(DcoNodeDesc * child_node) const {
  // get model the node belongs
  DcoModel * model = getModel();
  // get message handler and messages for loging
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  // get description of this node
  DcoNodeDesc * node_desc = getDesc();

  // solver number of columns
  int num_cols = model->solver()->getNumCols();

  // this will keep bound information, lower and upper
  // it is used for both hard and soft bounds.
  Bound bound;
  bound.lower.ind = new int[num_cols];
  bound.lower.val = new double[num_cols];
  bound.upper.ind = new int[num_cols];
  bound.upper.val = new double[num_cols];

  // Hard bounds
  // == Hard lower bounds
  std::copy(node_desc->getVars()->lbHard.posModify,
            node_desc->getVars()->lbHard.posModify + num_cols,
            bound.lower.ind);
  std::copy(node_desc->getVars()->lbHard.entries,
            node_desc->getVars()->lbHard.entries + num_cols,
            bound.lower.val);
  // == Hard upper bounds
  std::copy(node_desc->getVars()->ubHard.posModify,
            node_desc->getVars()->ubHard.posModify + num_cols,
            bound.upper.ind);
  std::copy(node_desc->getVars()->ubHard.entries,
            node_desc->getVars()->ubHard.entries + num_cols,
            bound.upper.val);

  // assign hard bounds for the child node.
  // this takes ownership of the arrays.
  child_node->assignVarHardBound(num_cols,
                                 bound.lower.ind,
                                 bound.lower.val,
                                 num_cols,
                                 bound.upper.ind,
                                 bound.upper.val);

  // Soft bounds.
  // == Soft lower bounds
  // number of entries modified for soft lower bounds
  int sl_num_modify = node_desc->getVars()->lbSoft.numModify;
  // We can reuse bound here since previos arrays are owned by Bcps.
  bound.lower.ind = new int[sl_num_modify];
  bound.lower.val = new double[sl_num_modify];
  std::copy(node_desc->getVars()->lbSoft.posModify,
            node_desc->getVars()->lbSoft.posModify + sl_num_modify,
            bound.lower.ind);
  std::copy(node_desc->getVars()->lbSoft.entries,
            node_desc->getVars()->lbSoft.entries + sl_num_modify,
            bound.lower.val);

  // == Soft upper bounds
  // number of entries modified for soft upper bounds
  int su_num_modify = node_desc->getVars()->ubSoft.numModify;
  bound.upper.ind = new int[su_num_modify];
  bound.upper.val = new double[su_num_modify];
  std::copy(node_desc->getVars()->ubSoft.posModify,
            node_desc->getVars()->ubSoft.posModify + su_num_modify,
            bound.upper.ind);
  std::copy(node_desc->getVars()->ubSoft.entries,
            node_desc->getVars()->ubSoft.entries + su_num_modify,
            bound.upper.val);

  // Assign soft bounds for the child node.
  child_node->assignVarSoftBound(sl_num_modify,
                                 bound.lower.ind,
                                 bound.lower.val,
                                 su_num_modify,
                                 bound.upper.ind,
                                 bound.upper.val);
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

void DcoTreeNode::branchConstrainOrPrice(DcoSubproblemStatus subproblem_status,
                                         bool & keepBounding,
                                         bool & branch,
                                         bool & generateConstraints,
                                         bool & generateVariables) {
  DcoModel * model = getModel();
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;

  // check whether the subproblem solved properly, solver status should be
  // infeasible or optimal.
  if (subproblem_status==DcoSubproblemStatusPrimalInfeasible or
      subproblem_status==DcoSubproblemStatusDualInfeasible) {
    // subproblem is infeasible, fathom
    keepBounding = false;
    setStatus(AlpsNodeStatusFathomed);
    return;
  }
  if (subproblem_status!=DcoSubproblemStatusOptimal) {
    message_handler->message(DISCO_SOLVER_FAILED, *messages)
      << CoinMessageEol;
    throw std::exception();
  }

  // subproblem is solved to optimality. Check feasibility of the solution.
  bool feasible = true;
  int numColsInf;
  int numRowsInf;
  branch = true;
  DcoSolution * sol = model->feasibleSolution(numColsInf, numRowsInf);
  if (numColsInf) {
    feasible = false;
    // default strategy is to branch, for now
    keepBounding = false;
    generateVariables = false;
    generateConstraints = false;
    message_handler->message(0, "Dco", "Some columns are infeasible. "
                             "Decided to branch.", 'G', DISCO_DLOG_BRANCH)
      << CoinMessageEol;
    return;

  }
  else {
    // all columns are feasible
    // check conic feasibility
    message_handler->message(0, "Dco", "All columns are feasible. "
                             "Checking conic feasibility...",
                             'G', DISCO_DLOG_BRANCH)
      << CoinMessageEol;
    // numRowInf will have the number of infeasible conic constraints.
    // since we relax rows corresponding to conic constraints only.
    if (numRowsInf) {
      message_handler->message(0, "Dco", "Some conic constraints are "
                               "infeasible. Decided to generate cuts.",
                               'G', DISCO_DLOG_BRANCH)
        << CoinMessageEol;

      feasible = false;
      // all cols are feasible, keep bounding and generateConstraints
      keepBounding = true;
      generateVariables = false;
      generateConstraints = true;
      return;
    }
  }
  if (sol) {
    model->storeSolution(sol);
    // set node status to fathomed.
    message_handler->message(0, "Dco", "Node is feasible, fathoming... ",
                             'G', DISCO_DLOG_BRANCH)
      << CoinMessageEol;
    feasible = true;
    // default strategy is to branch, for now
    keepBounding = false;
    generateVariables = false;
    generateConstraints = false;
    setStatus(AlpsNodeStatusFathomed);
  }
  else {
    message_handler->message(9998, "Dco", "This should not happen.",
                             'E')
      << CoinMessageEol;
  }

  // if still feasible then check solution quality, store solution
  // set status to fathom
}

//todo(aykut) replace this with DcoModel::feasibleSolution????
void DcoTreeNode::checkRelaxedCols(int & numInf) {
  DcoModel * model = getModel();
  double const * sol = model->solver()->getColSolution();
  // get integer tolerance
  double tol = model->dcoPar()->entry(DcoParams::integerTol);
  numInf = 0;
  int numRelaxedCols = model->numRelaxedCols();
  int const * relaxedCols = model->relaxedCols();
  //std::vector<BcpsVariable*> & cols = getVariables();
  // iterate over relaxed columns
  for (int i=0; i<numRelaxedCols; ++i) {
    double value = sol[relaxedCols[i]];
    double down = floor(value);
    double up = ceil(value);
    if ((value-down)<tol or (up-value)<tol) {
      continue;
    }
    else {
      numInf++;
    }
  }
}

void DcoTreeNode::applyConstraints(BcpsConstraintPool const & conPool) {
  DcoModel * model = getModel();
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  double scale_par = model->dcoPar()->entry(DcoParams::scaleConFactor);

  double const * sol = model->solver()->getColSolution();

  // Tranform constraints to Osi cut so that easily add them to LP.
  int num_cuts = conPool.getNumConstraints();
  OsiRowCut const ** cuts_to_add = new OsiRowCut const * [num_cuts];
  int num_add = 0;
  std::vector<int> cuts_to_del;

  // iterate over cuts and
  //------------------------------------------
  // Remove:
  //  - empty cuts
  //  - dense cuts
  //  - bad scaled cuts
  //  - weak cuts
  //  - parallel cuts
  //------------------------------------------
  for (int i=0; i<num_cuts; ++i) {
    DcoLinearConstraint * curr_con =
      dynamic_cast<DcoLinearConstraint*>(conPool.getConstraint(i));
    int length = curr_con->getSize();
    double const * elements = curr_con->getValues();
    int const * indices = curr_con->getIndices();

    // check whether cut is empty.
    if (length <= 0) {
      cuts_to_del.push_back(i);
      continue;
    }
    // check whether cut is dense.
    if(length > 100) {
      // discard the cut
      cuts_to_del.push_back(i);
      continue;
    }

    // check cut scaling
    double activity = 0.0;
    double maxElem = 0.0;
    double minElem = ALPS_DBL_MAX;
    double scaleFactor;

    // iterate over coefficients
    for (int k=0; k<length; ++k) {
      if (elements[k]==0.0) {
        // if a coef is exactly 0, ignore it.
        // what if all coef is zero?
        cuts_to_del.push_back(i);
        continue;
      }
      if (fabs(elements[k]) > maxElem) {
        maxElem = fabs(elements[k]);
      }
      if (fabs(elements[k]) < minElem) {
        minElem = fabs(elements[k]);
      }
      activity += elements[k] * sol[indices[k]];
    }
    if (maxElem==0.0) {
      // all coefficients are 0, skip cut
      cuts_to_del.push_back(i);
      continue;
    }
    scaleFactor = maxElem/minElem;

    if (scaleFactor > scale_par) {
      // skip the cut since it is badly scaled.
      cuts_to_del.push_back(i);
      continue;
    }

    // Check whether the cut is weak.
    double rowLower = CoinMax(curr_con->getLbHard(),
                              curr_con->getLbSoft());
    double rowUpper = CoinMin(curr_con->getUbHard(),
                              curr_con->getUbSoft());
    // violation will be positive if the cut cuts the solution.
    double violation = -1.0;
    if (rowLower > -ALPS_INFINITY) {
      violation = rowLower - activity;
    }
    if (rowUpper < ALPS_INFINITY) {
      violation = CoinMax(violation, activity-rowUpper);
    }

    // todo(aykut) this should be a parameter.
    if (violation < 1.0e-8) {
      // cut is weak, skip it.
      cuts_to_del.push_back(i);
      continue;
    }

    // todo(aykut)
    // Check whether cut is parallel to an existing constraint or cut.
    bool parallel = false;
    if (parallel) {
      cuts_to_del.push_back(i);
      continue;
    }
    cuts_to_add[num_add++] = curr_con->createOsiRowCut(model);
  }

  // Add cuts to lp and adjust basis.
  CoinWarmStartBasis * ws = dynamic_cast<CoinWarmStartBasis*>
    (model->solver()->getWarmStart());

  int num_solver_rows = model->solver()->getNumRows();
  int num_solver_cols = model->solver()->getNumCols();

  if (num_add > 0) {
    model->solver()->applyRowCuts(num_add, cuts_to_add);
    ws->resize(num_solver_rows + num_add, num_solver_cols);
    for (int i=0; i<num_add; ++i) {
      ws->setArtifStatus(num_solver_rows + i,
                         CoinWarmStartBasis::basic);
    }
    if (model->solver()->setWarmStart(ws) == false) {
      throw CoinError("Fail setWarmStart() after cut installation.",
                      "applyConstraints","DcoTreeNode");
    }
    for (int k=0; k<num_add; ++k) {
      delete cuts_to_add[k];
    }
  }

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << num_add;
  debug_msg << " out of ";
  debug_msg << conPool.getNumConstraints();
  debug_msg << "  cuts added to the solver.";
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                           'G', DISCO_DLOG_CUT)
    << CoinMessageEol;
  // end of debug stuff


  delete[] cuts_to_add;
  delete ws;
}
