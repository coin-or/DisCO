/*===========================================================================*
 * This file is part of the Discrete Conic Optimization (DisCO) Solver.      *
 *                                                                           *
 * DisCO is distributed under the Eclipse Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *          Aykut Bulut, Lehigh University                                   *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and      *
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


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
#include "DcoVariable.hpp"
#include "Dco.hpp"
#include "DcoConGenerator.hpp"
#include "DcoSolution.hpp"
#include "DcoBranchObject.hpp"
#include "DcoHeuristic.hpp"
#include "DcoHeurRounding.hpp"

// STL headers
#include <vector>
#include <numeric>

extern std::map<DISCO_Grumpy_Msg_Type, char const *> grumpyMessage;
extern std::map<DcoNodeBranchDir, char> grumpyDirection;

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
  desc_ = NULL;
  bcpStats_.numBoundIter_ = 0;
  bcpStats_.numTotalCuts_ = 0;
  bcpStats_.numLastCuts_ = 0;
  bcpStats_.numOaCuts_ = 0;
  bcpStats_.numMilpCuts_ = 0;
  bcpStats_.totalImp_ = 0.0;
  bcpStats_.lastImp_ = 0.0;
  bcpStats_.lastObjVal_ = 0.0;
  bcpStats_.startObjVal_ = 0.0;
  bcpStats_.numMilpIter_ = 0;
  bcpStats_.inactive_.clear();
  bcpStats_.generatorIndex_.clear();
}

DcoTreeNode::DcoTreeNode(AlpsNodeDesc *& desc) : BcpsTreeNode() {
  desc_ = desc;
  desc = NULL;
  bcpStats_.numBoundIter_ = 0;
  bcpStats_.numTotalCuts_ = 0;
  bcpStats_.numLastCuts_ = 0;
  bcpStats_.numOaCuts_ = 0;
  bcpStats_.numMilpCuts_ = 0;
  bcpStats_.totalImp_ = 0.0;
  bcpStats_.lastImp_ = 0.0;
  bcpStats_.lastObjVal_ = 0.0;
  bcpStats_.startObjVal_ = 0.0;
  bcpStats_.numMilpIter_ = 0;
  bcpStats_.inactive_.clear();
  bcpStats_.generatorIndex_.clear();
}

DcoTreeNode::~DcoTreeNode() {
}

// create tree node from given description
AlpsTreeNode * DcoTreeNode::createNewTreeNode(AlpsNodeDesc *& desc) const {
  DcoNodeDesc * dco_node = dynamic_cast<DcoNodeDesc*>(desc);
  int branched_index = dco_node->getBranchedInd();
  double branched_value = dco_node->getBranchedVal();

  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());

  double dist_to_floor = branched_value - floor(branched_value);
  double dist_to_ceil = ceil(branched_value) - branched_value;
  double frac = (dist_to_floor>dist_to_ceil) ? dist_to_ceil : dist_to_floor;
  double tol = model->dcoPar()->entry(DcoParams::integerTol);

  if (frac < tol) {
    model->dcoMessageHandler_->message(DISCO_NODE_BRANCHONINT,
                                       *(model->dcoMessages_))
      << broker()->getProcRank()
      << branched_index << CoinMessageEol;
  }

  // Create a new tree node
  DcoTreeNode * node = new DcoTreeNode(desc);
  node->setBroker(broker_);
  desc = NULL;
  return node;
}

void DcoTreeNode::convertToExplicit() {
  DcoNodeDesc * node_desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  //CoinMessages * messages = model->dcoMessages_;

  if (explicit_) {
    return;
  }
  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "["
            << broker_->getProcRank()
            << "] converting "
            << this
            << " to explicit.";
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                           'G', DISCO_DLOG_MPI)
    << CoinMessageEol;
  // end of debug stuff

  explicit_ = 1;

  int num_cols = model->solver()->getNumCols();
  // this will keep bound information, lower and upper
  Bound hard_bound;
  // ownership is transferred to Bcps, no need to free.
  hard_bound.lower.ind = new int[num_cols];
  hard_bound.lower.val = new double[num_cols];
  hard_bound.upper.ind = new int[num_cols];
  hard_bound.upper.val = new double[num_cols];
  std::fill_n(hard_bound.lower.val, num_cols, ALPS_DBL_MAX);
  std::fill_n(hard_bound.upper.val, num_cols, -ALPS_DBL_MAX);

  Bound soft_bound;
  // ownership is transferred to Bcps, no need to free.
  soft_bound.lower.ind = new int[num_cols];
  soft_bound.lower.val = new double[num_cols];
  soft_bound.upper.ind = new int[num_cols];
  soft_bound.upper.val = new double[num_cols];
  std::fill_n(soft_bound.lower.val, num_cols, ALPS_DBL_MAX);
  std::fill_n(soft_bound.upper.val, num_cols, -ALPS_DBL_MAX);

  for (int k = 0; k < num_cols; ++k) {
    hard_bound.lower.ind[k] = k;
    hard_bound.upper.ind[k] = k;
    soft_bound.lower.ind[k] = k;
    soft_bound.upper.ind[k] = k;
  }

  //--------------------------------------------------
  // Travel back to a full node, then collect diff (add/rem col/row,
  // hard/soft col/row bounds) from the node full to this node.
  //--------------------------------------------------------

  AlpsTreeNode * parent = parent_;

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
  for(int i=static_cast<int> (leafToRootPath.size()-1); i>-1; --i) {
    DcoNodeDesc * curr = dynamic_cast<DcoNodeDesc*>((leafToRootPath.at(i))->
                                           getDesc());
    //--------------------------------------
    // Full variable hard bounds.
    //--------------------------------------
    int numModify = curr->getVars()->lbHard.numModify;
    for (int k=0; k<numModify; ++k) {
      int index = curr->getVars()->lbHard.posModify[k];
      double value = curr->getVars()->lbHard.entries[k];
      hard_bound.lower.val[index] = value;
    }

    numModify = curr->getVars()->ubHard.numModify;
    for (int k=0; k<numModify; ++k) {
      int index = curr->getVars()->ubHard.posModify[k];
      double value = curr->getVars()->ubHard.entries[k];
      hard_bound.upper.val[index] = value;
    }

    //--------------------------------------
    // Full variable soft bounds.
    //--------------------------------------
    numModify = curr->getVars()->lbSoft.numModify;
    for (int k=0; k<numModify; ++k) {
      int index = curr->getVars()->lbSoft.posModify[k];
      double value = curr->getVars()->lbSoft.entries[k];
      soft_bound.lower.val[index] = value;
    }

    numModify = curr->getVars()->ubSoft.numModify;
    for (int k=0; k<numModify; ++k) {
      int index = curr->getVars()->ubSoft.posModify[k];
      double value = curr->getVars()->ubSoft.entries[k];
      soft_bound.upper.val[index] = value;
    }



  } // EOF for (path)

  //------------------------------------------
  // Record hard variable bounds. FULL set.
  //------------------------------------------
  node_desc->assignVarHardBound(num_cols,
                                hard_bound.lower.ind,
                                hard_bound.lower.val,
                                num_cols,
                                hard_bound.upper.ind,
                                hard_bound.upper.val);

  //------------------------------------------
  // Recode soft variable bound. Modified.
  //------------------------------------------
  int numSoftVarLowers=0;
  int numSoftVarUppers=0;
  for (int k=0; k<num_cols; ++k) {
    if (soft_bound.lower.val[k] < ALPS_BND_MAX) {
      soft_bound.lower.ind[numSoftVarLowers] = k;
      soft_bound.lower.val[numSoftVarLowers++] = soft_bound.lower.val[k];
    }
    if (soft_bound.upper.val[k] > -ALPS_BND_MAX) {
      soft_bound.upper.ind[numSoftVarUppers] = k;
      soft_bound.upper.val[numSoftVarUppers++] = soft_bound.upper.val[k];
    }
  }
  // Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
  node_desc->assignVarSoftBound(numSoftVarLowers,
                                soft_bound.lower.ind,
                                soft_bound.lower.val,
                                numSoftVarUppers,
                                soft_bound.upper.ind,
                                soft_bound.upper.val);

  //--------------------------------------------------
  // Clear path vector.
  //--------------------------------------------------
  leafToRootPath.clear();
  assert(leafToRootPath.size() == 0);

}

void DcoTreeNode::convertToRelative() {
  //DcoNodeDesc * node_desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << broker()->getProcRank()
    << __FILE__ << __LINE__ << CoinMessageEol;
  throw std::exception();
}


// todo(aykut) I feel like branchConstraintPrice function and
// generateConstraint function should have more communication, since decisions
// of generating cuts and which type of cuts to generate should be made
// together.

/// type is bitmask to generate cuts
/// 001 generate OA
/// 010 generate MILP
/// 011 generate both
/// Generate constraints for the problem.
int DcoTreeNode::generateConstraints(BcpsConstraintPool * conPool, int type) {
  DcoModel * disco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = disco_model->dcoMessageHandler_;
  //CoinMessages * messages = disco_model->dcoMessages_;

  std::map<DcoConstraintType, DcoConGenerator*>::iterator it;
  std::map<DcoConstraintType, DcoConGenerator*> con_generators = disco_model->conGenerators();
  for (it=con_generators.begin(); it != con_generators.end(); ++it) {
    bool do_use = false;
    DcoConGenerator * cg = it->second;
    // decide whether we should use this cut generator with respect to the
    // specified cut strategy
    decide_using_cg(do_use, cg, type);
    if (!do_use) {
      // jump to the next geenrator if we will not use this one.
      continue;
    }
    int pre_num_cons = conPool->getNumConstraints();
    double start_time = CoinCpuTime();
    // Call constraint generator
    //bool must_resolve = cg->generateConstraints(*conPool);
    cg->generateConstraints(*conPool);
    double cut_time = CoinCpuTime() - start_time;
    // Statistics
    cg->stats().addTime(cut_time);
    cg->stats().addNumCalls(1);
    int num_cons_generated = conPool->getNumConstraints() - pre_num_cons;
    if (num_cons_generated == 0) {
      cg->stats().addNumNoConsCalls(1);
    }
    else {
      cg->stats().addNumConsGenerated(num_cons_generated);
    }

    // debug msg
    std::stringstream debug_msg;
    debug_msg << "[" << broker()->getProcRank() << "] Called "
              << cg->name() << ", generated "
              << num_cons_generated << " cuts in "
              << cut_time << " seconds.";
    message_handler->message(0, "Dco", debug_msg.str().c_str(),
                             'G', DISCO_DLOG_CUT)
      << CoinMessageEol;
    // end of debug
  }
  // return value will make sense when DcoTreeNode::process is implemented
  // in Bcps level.
  return 0;
}


//todo(aykut) if all columns are feasible, fix columns and call IPM. If the
// resulting solution is better update UB. Fathom otherwise.
void DcoTreeNode::decide_using_cg(bool & do_use,
                                  DcoConGenerator * cg,
                                  int type) const {
  DcoModel * model = dynamic_cast<DcoModel*> (broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  DcoCutStrategy strategy = cg->strategy();
  do_use = false;
  // if type & 1 generate OA
  // if type & 2 generate MILP
  if (cg->name().compare("OA") == 0) {
    do_use = (type & 1) ? true : false;
  }
  else {
    // MILP generator
    if (type & 2) {
      // generate
      if (strategy==DcoCutStrategyNone) {
        // do nothing,
        do_use = false;
      }
      else if (strategy==DcoCutStrategyRoot) {
        do_use = (depth_ == 0)  ? true : false;
      }
      else if (strategy==DcoCutStrategyAuto) {
        // moreover disable generator if it does not do well
        do_use = (index_ % cg->frequency() == 0) ? true : false;
      }
      else if (strategy==DcoCutStrategyPeriodic) {
        do_use = (index_ % cg->frequency() == 0) ? true : false;
      }
      else {
        message_handler->message(DISCO_UNKNOWN_CUTSTRATEGY, *messages)
          << broker()->getProcRank()
          << strategy << CoinMessageEol;
      }
    }
    else {
      do_use = false;
    }
  }
}


int DcoTreeNode::generateVariables(BcpsVariablePool * varPool) {
  DcoModel * disco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = disco_model->dcoMessageHandler_;
  CoinMessages * messages = disco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << broker()->getProcRank()
    << __FILE__ << __LINE__ << CoinMessageEol;
  throw std::exception();
  return 0;
}

int DcoTreeNode::chooseBranchingObject() {
  DcoModel * disco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = disco_model->dcoMessageHandler_;
  CoinMessages * messages = disco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << broker()->getProcRank()
    << __FILE__ << __LINE__ << CoinMessageEol;
  throw std::exception();
  return 0;
}

int DcoTreeNode::process(bool isRoot, bool rampUp) {
  AlpsNodeStatus status = getStatus();
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;

  // debug stuff
#ifdef DISCO_DEBUG
  std::stringstream debug_msg;
  debug_msg << "["
            << broker()->getProcRank()
            << "] "
            << "Processing node "
            << this
            << " index "
            << getIndex()
            << " parent "
            << getParent();
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                           'G', DISCO_DLOG_PROCESS)
    << CoinMessageEol;
#endif
  // end of debug stuff

  // std::cout << "Broker reports quality " << broker()->getIncumbentValue()
  //           << std::endl;

  double cutoff = model->dcoPar()->entry(DcoParams::cutoff);
  double sense = model->dcoPar()->entry(DcoParams::objSense);
  cutoff = sense*cutoff;
  cutoff = CoinMin(cutoff, broker()->getIncumbentValue());
  // check if this can be fathomed
  double rel_gap_limit = model->dcoPar()->entry(DcoParams::optimalRelGap);
  double abs_gap_limit = model->dcoPar()->entry(DcoParams::optimalAbsGap);
  double abs_gap = cutoff-getQuality();
  double rel_gap = abs_gap/fabs(cutoff);
  // std::cout << "cutoff " << cutoff << std::endl;
  // std::cout << "abs " << abs_gap << " limit " << abs_gap_limit << std::endl;
  // std::cout << "rel " << rel_gap << " limit " << rel_gap_limit << std::endl;
  if (rel_gap_limit>rel_gap or abs_gap_limit>abs_gap) {
#ifdef DISCO_DEBUG
    // debug message
    message_handler->message(DISCO_NODE_FATHOM_PARENTQ, *model->dcoMessages_)
      << broker()->getProcRank()
      << getIndex()
      << abs_gap
      << rel_gap
      << CoinMessageEol;

    // end of debug message
    // grumpy message
    model->dcoMessageHandler_->message(DISCO_GRUMPY_MESSAGE_MED,
                                       *model->dcoMessages_)
      << broker()->getProcRank()
      << broker()->timer().getTime()
      << grumpyMessage[DISCO_GRUMPY_FATHOMED]
      << getIndex()
      << getParentIndex()
      << grumpyDirection[dynamic_cast<DcoNodeDesc*>(desc_)->getBranchedDir()]
      << model->objSense()*getQuality()
      << CoinMessageEol;
#endif
    // end of grumpy message
    setStatus(AlpsNodeStatusFathomed);
    return AlpsReturnStatusOk;
  }

  if (status==AlpsNodeStatusCandidate or
      status==AlpsNodeStatusEvaluated) {
    boundingLoop(isRoot, rampUp);
  }
  else if (status==AlpsNodeStatusBranched or
           status==AlpsNodeStatusFathomed or
           status==AlpsNodeStatusDiscarded) {
    // this should not happen
    message_handler->message(DISCO_NODE_UNEXPECTEDSTATUS, *messages)
      << broker()->getProcRank()
      << static_cast<int>(status)
      << CoinMessageEol;
  }
  return AlpsReturnStatusOk;
}

int DcoTreeNode::boundingLoop(bool isRoot, bool rampUp) {
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  bool keepBounding = true;
  bool do_branch = false;
  int genConstraints = 0;
  bool genVariables = false;
  BcpsConstraintPool * constraintPool = new BcpsConstraintPool();
  BcpsVariablePool * variablePool = new BcpsVariablePool();
  installSubProblem();

  while (keepBounding) {
    keepBounding = false;
    // solve subproblem corresponds to this node
    BcpsSubproblemStatus subproblem_status = bound();
    // update number of iterations statistics
    model->addNumRelaxIterations();
    // update bcp statistics
    if (bcpStats_.numBoundIter_==0) {
      bcpStats_.startObjVal_ = model->solver()->getObjValue();
      bcpStats_.lastObjVal_ = model->solver()->getObjValue();
    }
    else {
      bcpStats_.totalImp_ = model->objSense()*
        (model->solver()->getObjValue()-bcpStats_.startObjVal_);
      bcpStats_.lastImp_ = model->objSense()*
        (model->solver()->getObjValue()-bcpStats_.lastObjVal_);
      bcpStats_.lastObjVal_ = model->solver()->getObjValue();
    }
#ifdef DISCO_DEBUG
    CoinMessages * messages = model->dcoMessages_;
    // debug print solver status
    message_handler->message(DISCO_SOLVER_STATUS, *messages)
      << broker()->getProcRank()
      << subproblem_status
      << model->solver()->getObjValue()
      << quality_
      << solEstimate_
      << CoinMessageEol;
    // end of debug stuff
    // grumpy message
    if ((subproblem_status==BcpsSubproblemStatusOptimal) &&
        (getStatus()==AlpsNodeStatusCandidate or
         getStatus()==AlpsNodeStatusEvaluated)) {
      double sum_inf = 0.0;
      int num_inf = 0;
      int const * relaxedCols = model->relaxedCols();
      for (int i=0; i<model->numRelaxedCols(); ++i) {
        int dir;
        double infeas = dynamic_cast<DcoVariable*>
          (model->getVariables()[relaxedCols[i]])->infeasibility(model, dir);
        if (infeas!=0.0) {
          num_inf++;
        }
        sum_inf += infeas;
      }
      message_handler->message(DISCO_GRUMPY_MESSAGE_LONG, *messages)
        << broker()->getProcRank()
        << broker()->timer().getTime()
        << grumpyMessage[DISCO_GRUMPY_CANDIDATE]
        << getIndex()
        << getParentIndex()
        << grumpyDirection[dynamic_cast<DcoNodeDesc*>(desc_)->getBranchedDir()]
        << model->objSense()*getQuality()
        << sum_inf
        << num_inf
        << CoinMessageEol;
    }
    // end of grumpy message
#endif

    // check if this can be fathomed
    double cutoff = model->dcoPar()->entry(DcoParams::cutoff);
    double sense = model->dcoPar()->entry(DcoParams::objSense);
    cutoff = sense*cutoff;
    double rel_gap_limit = model->dcoPar()->entry(DcoParams::optimalRelGap);
    double abs_gap_limit = model->dcoPar()->entry(DcoParams::optimalAbsGap);
    cutoff = CoinMin(cutoff, broker()->getIncumbentValue());
    double abs_gap = cutoff-getQuality();
    double rel_gap = abs_gap/fabs(cutoff);
    //std::cout << "abs " << abs_gap << " limit " << abs_gap_limit << std::endl;
    //std::cout << "rel " << rel_gap << " limit " << rel_gap_limit << std::endl;
    if (rel_gap_limit>rel_gap or abs_gap_limit>abs_gap) {
#ifdef DISCO_DEBUG
      // debug message
      message_handler->message(DISCO_NODE_FATHOM, *messages)
        << broker()->getProcRank()
        << getIndex()
        << abs_gap
        << rel_gap
        << CoinMessageEol;

      // grumpy message
      model->dcoMessageHandler_->message(DISCO_GRUMPY_MESSAGE_MED,
                                         *model->dcoMessages_)
        << broker()->getProcRank()
        << broker()->timer().getTime()
        << grumpyMessage[DISCO_GRUMPY_FATHOMED]
        << getIndex()
        << getParentIndex()
        << grumpyDirection[dynamic_cast<DcoNodeDesc*>(desc_)->getBranchedDir()]
        << model->objSense()*getQuality()
        << CoinMessageEol;
      // end of grumpy message
#endif
      setStatus(AlpsNodeStatusFathomed);
      break;
    }
    // call heuristics to search for a solution
    if (bcpStats_.numBoundIter_==0) {
      callHeuristics();
    }
#ifdef __OA__
    // update cut stats and clean if necessary
    checkCuts();
#endif

    // decide what to do
    branchConstrainOrPrice(subproblem_status, keepBounding, do_branch,
                           genConstraints,
                           genVariables);
#ifdef DISCO_DEBUG
    // debug message
    message_handler->message(DISCO_NODE_BCP_DECISION, *messages)
      << broker()->getProcRank()
      << keepBounding
      << do_branch
      << genConstraints
      << CoinMessageEol;
    // end of debug stuff
#endif
    if (getStatus()==AlpsNodeStatusFathomed) {
      // node is fathomed, nothing to do.
      break;
    }
    else if (keepBounding and genConstraints) {
      generateConstraints(constraintPool, genConstraints);
      // add constraints to the model
      applyConstraints(constraintPool);
      // clear constraint pool
      constraintPool->freeGuts();
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding and genVariables) {
      generateVariables(variablePool);
      // add variables to the model
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding==false and do_branch==false) {
      // put node back into the list.
      // this means do not change the node status and end processing the node.
      // set status to evaluated
      setStatus(AlpsNodeStatusEvaluated);
    }
    else if (keepBounding==false and do_branch) {
      // branch
      BcpsBranchStrategy * branchStrategy = model->branchStrategy();
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
    // increase iteration number
    bcpStats_.numBoundIter_++;
  }
  delete constraintPool;
  delete variablePool;
  return AlpsReturnStatusOk;
}

//1. update cut statistics
//2. remove inactive cuts
void DcoTreeNode::checkCuts() {
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  //CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  //CoinMessages * messages = model->dcoMessages_;
  int origNumRows = model->getNumCoreLinearConstraints();
  int initOAcuts = model->initOAcuts();
  int solverNumRows = model->solver()->getNumRows();
  int numCuts = solverNumRows - origNumRows;
  if (numCuts==0) {
    return;
  }
  CoinWarmStartBasis * ws =
    dynamic_cast<CoinWarmStartBasis*> (model->solver()->getWarmStart());
  if (ws==NULL) {
    // nothing to do if there is no warm start information
    std::cerr << "Disco warning: No warm start object exists in solver. "
              << "Unable to clean cuts." << std::endl;
    return;
  }
  double cutOaSlack = model->dcoPar()->entry(DcoParams::cutOaSlack2);
  BcpStats * st = &bcpStats_;
  // iterate over cuts and update inactivity statistics
  std::list<int>::iterator curr = st->inactive_.begin();
  for (int i=0; i<numCuts; ++i) {

      // double norm = 0.0;
      // {
      //   // compute norm of cut
      //   CoinPackedMatrix const * mat = solver_->getMatrixByRow();
      //   int start = mat->getVectorStarts()[i];
      //   int size = mat->getVectorLengths()[i];
      //   double const * value = mat->getElements() + start;
      //   for (int j=0; j<size; ++j) {
      //     norm += fabs(value[j]);
      //   }
      // }


    double slack = model->solver()->getRowUpper()[origNumRows+i] - model->solver()->getRowActivity()[origNumRows+i];
    // check whether ith cut is basic
    if (ws->getArtifStatus(origNumRows+i) == CoinWarmStartBasis::basic &&
        slack > cutOaSlack) {
      // print activity
      // std::cout << "activity: "
      //           << slack
      //           << std::endl;
      // cut is inactive, update stats
      if (curr==st->inactive_.end()) {
        // no entry in the list for cut i, insert
        st->inactive_.insert(curr, 1);
      }
      else {
        // update entry of cut i
        *curr += 1;
        curr++;
      }
    }
    else {
      // cut is active
      if (curr==st->inactive_.end()) {
        // insert if there is no entry for cut i yet
        st->inactive_.insert(curr, 0);
      }
      else {
        // reset stat for cut i to 0
        *curr = 0;
        curr++;
      }
    }
  }
  // iterate over cuts and remove the ones that are inactive for
  // 3 or more iterations
  //int slackLimit = model->dcoPar()->entry(DcoParams::cutOaSlackLimit);
  {
    std::list<int>::iterator curr = st->inactive_.begin();
    //std::list<int>::iterator curr2 = st->generatorIndex_.begin();
    int numDel = 0;
    int * delInd = new int[numCuts];
    int num_del_init_oa = 0;
    for (int i=0; i<numCuts; ++i) {
      if (*curr>3) {
        delInd[numDel++] = i+origNumRows;
        // check whether cut is an initial OA cut
        if (i<initOAcuts) {
          num_del_init_oa++;
        }
        // delete curr from inactivity_
        curr = st->inactive_.erase(curr);
        //curr2 = st->inactive_.erase(curr2);
      }
      else {
        curr++;
        //curr2++;
      }
    }
    // delete cuts from solver
    if (numDel) {
      // std::cout << "node " << index_
      //           << " iter: " << st->numBoundIter_
      //           << " removed: " << numDel
      //           << " remain: " << numCuts-numDel
      //           << std::endl;
      model->solver()->deleteRows(numDel, delInd);
      // resolve to correct status
      model->solver()->resolve();
    }
    delete[] delInd;
    model->decreaseInitOAcuts(num_del_init_oa);
  }
  delete ws;
}

void DcoTreeNode::callHeuristics() {
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  long unsigned int num_heur = model->numHeuristics();
  DcoSolution * sol = NULL;
  for (long unsigned int i=0; i<num_heur; ++i) {
    DcoHeuristic * curr = model->heuristics(i);
    double start_time = CoinCpuTime();
    sol = curr->searchSolution();
    double heur_time = CoinCpuTime() - start_time;
    // Statistics
    curr->stats().addTime(heur_time);
    curr->stats().addCalls(1);
    if (sol) {
      double cutoff = model->dcoPar()->entry(DcoParams::cutoff);
      double sense = model->dcoPar()->entry(DcoParams::objSense);
      cutoff = sense*cutoff;
      cutoff = CoinMin(cutoff, sol->getQuality());
      if (sol->getQuality() <= cutoff) {
        // set depth
        setDepth(depth_);
        // set index
        setIndex(index_);
        // Store in Alps pool
        broker()->addKnowledge(AlpsKnowledgeTypeSolution,
                               sol,
                               model->objSense() * sol->getQuality());
        double incum_value = broker()->getIncumbentValue();
        model->solver()->setDblParam(OsiDualObjectiveLimit,
                                     model->objSense()*incum_value);
        // debug log
        message_handler->message(DISCO_HEUR_SOL_FOUND, *messages)
          << broker()->getProcRank()
          << curr->name()
          << sol->getQuality();
      }
    }
    else {
      // debug message
      message_handler->message(DISCO_HEUR_NOSOL_FOUND, *messages)
        << broker()->getProcRank()
        << curr->name();
    }
  }

  if (sol) {
    //double comp_obj = 0.0;
    //int num_cols = model->solver()->getNumCols();
    //double const * obj_coef = model->solver()->getObjCoefficients();
    //comp_obj = std::inner_product(sol->getValues(), sol->getValues()+num_cols, obj_coef, 0.0);
    // std::cout << "==================== "
    //           << "objective "
    //           << "===================="
    //           << std::endl
    //           << comp_obj
    //           << std::endl;
  }
}

/// Bounding procedure to estimate quality of this node.
BcpsSubproblemStatus DcoTreeNode::bound() {
  BcpsSubproblemStatus subproblem_status = BcpsSubproblemStatusUnknown;
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
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
    subproblem_status = BcpsSubproblemStatusAbandoned;
  }
  else if (model->solver()->isProvenOptimal()) {
    // todo(aykut) if obj val is greater than 1e+30 we consider problem
    // infeasible. This is due to our cut generation when the problem is
    // infeasible. We add a high lower bound (1e+30) to the objective
    // function. This can be improved.
    if (model->solver()->getObjValue()>=1e+30) {
      subproblem_status = BcpsSubproblemStatusPrimalInfeasible;
      quality_ = ALPS_OBJ_MAX;
      solEstimate_ = ALPS_OBJ_MAX;
    }
    else {
      subproblem_status = BcpsSubproblemStatusOptimal;
      double new_quality = model->solver()->getObjValue() *
        model->solver()->getObjSense();
      // Update quality of this node if this bound is better than the
      // current quality_ of the node. quality_ is set to parent's quality_
      // when a node is created. The objValuexobjSense here might be larger
      // than the quality_ since not all conic cuts present in parent node is
      // present here.
      if (parent_==NULL) {
        // this is root node, update quality
        quality_ = new_quality;
      }
      else if (new_quality>quality_) {
        quality_ = new_quality;
        solEstimate_ = new_quality;
      }
    }
  }
  else if (model->solver()->isProvenPrimalInfeasible()) {
    subproblem_status = BcpsSubproblemStatusPrimalInfeasible;
  }
  else if (model->solver()->isProvenDualInfeasible()) {
    subproblem_status = BcpsSubproblemStatusDualInfeasible;
  }
  else if (model->solver()->isPrimalObjectiveLimitReached()) {
    subproblem_status = BcpsSubproblemStatusPrimalObjLim;
  }
  else if (model->solver()->isDualObjectiveLimitReached()) {
    subproblem_status = BcpsSubproblemStatusDualObjLim;
  }
  else if (model->solver()->isIterationLimitReached()) {
    subproblem_status = BcpsSubproblemStatusIterLim;
  }
  else {
    message_handler->message(DISCO_SOLVER_UNKNOWN_STATUS, *messages)
      << broker()->getProcRank()
      << CoinMessageEol;
  }
  return subproblem_status;
}

int DcoTreeNode::installSubProblem() {
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
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  DcoNodeDesc * desc = dynamic_cast<DcoNodeDesc*>(getDesc());
  // get number of columns and rows
  int numCoreCols = model->getNumCoreVariables();
  int numCoreLinearRows = model->getNumCoreLinearConstraints();
  //int numCoreConicRows = model->getNumCoreConicConstraints();
  // get number of columns and rows stored in the solver
  //int numSolverCols = model->solver()->getNumCols();
  // solver rows are all linear, for both Osi and OsiConic.
  int numSolverRows = model->solver()->getNumRows();
  int initOAcuts = model->initOAcuts();

  // 1. Remove noncore columns and rows
  // 1.1 Remove non-core rows from solver, i.e. cuts

  // we want to remove all cuts that are not added initialy after approximating
  // cones in #DcoModel::approximateCones(). Approximation cuts that are added
  // initially comes first in the order in solver rows.
  int numDelRows = numSolverRows - numCoreLinearRows - initOAcuts;
#ifndef __COLA__
  if (numDelRows > 0) {
    int * indices = new int[numDelRows];
    if (indices==NULL) {
      message_handler->message(DISCO_OUT_OF_MEMORY, *messages)
        << broker()->getProcRank()
        << __FILE__ << __LINE__ << CoinMessageEol;
      throw CoinError("Out of memory", "installSubProblem", "DcoTreeNode");
    }
    for (int i=0; i<numDelRows; ++i) {
      indices[i] = numCoreLinearRows + initOAcuts + i;
    }
    model->solver()->deleteRows(numDelRows, indices);
    delete[] indices;
    indices = NULL;
  }
#endif
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
  // this creates a bug unless you restore the bounds to original values.
  // this should be fixed.
  double * colLB = model->colLB();
  double * colUB = model->colUB();
  CoinFillN(colLB, numCoreCols, -ALPS_DBL_MAX);
  CoinFillN(colUB, numCoreCols, ALPS_DBL_MAX);
  // generate path to root from this
  std::vector<AlpsTreeNode*> leafToRootPath;
  leafToRootPath.push_back(this);
  if (broker_->getPhase() != AlpsPhaseRampup) {
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
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  // get message handler and messages for loging
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;

  // check node status, this should be a pregnant node.
  if (getStatus()!=AlpsNodeStatusPregnant) {
      message_handler->message(DISCO_NODE_UNEXPECTEDSTATUS, *messages)
        << broker()->getProcRank()
        << static_cast<int>(getStatus()) << CoinMessageEol;
  }

  // create return value and push the down and up nodes.
  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > res;


  // check if this can be fathomed
  double rel_gap_limit = model->dcoPar()->entry(DcoParams::optimalRelGap);
  double abs_gap_limit = model->dcoPar()->entry(DcoParams::optimalAbsGap);
  double abs_gap = broker()->getIncumbentValue()-getQuality();
  double rel_gap = abs_gap/fabs(broker()->getIncumbentValue());
  //std::cout << "abs " << abs_gap << " limit " << abs_gap_limit << std::endl;
  //std::cout << "rel " << rel_gap << " limit " << rel_gap_limit << std::endl;
  if (rel_gap_limit>rel_gap or abs_gap_limit>abs_gap) {
#ifdef DISCO_DEBUG
    message_handler->message(DISCO_NODE_FATHOM, *messages)
      << broker()->getProcRank()
      << getIndex()
      << abs_gap
      << rel_gap
      << CoinMessageEol;
#endif
    setStatus(AlpsNodeStatusFathomed);
    return res;
  }
  //todo(aykut) update scores
  //BcpsBranchStrategy * branchStrategy = model->branchStrategy();

  // get Alps phase
  AlpsPhase phase = broker()->getPhase();

  // get branch object
  DcoBranchObject const * branch_object =
    dynamic_cast<DcoBranchObject const *>(branchObject());

  assert(branch_object);
  // get index and value of branch variable.
  //int branch_var = model->relaxedCols()[branch_object->getObjectIndex()];
  int branch_var = branch_object->index();
  double branch_value = branch_object->value();

#ifdef DISCO_DEBUG
  message_handler->message(DISCO_NODE_BRANCH, *messages)
    << broker()->getProcRank()
    << index_
    << branch_var
    << branch_value
    << branch_object->score()
    << CoinMessageEol;
#endif

  // compute child nodes' warm start basis
  // ownership is transferred to node description object, no need to free.
  CoinWarmStartBasis * child_ws;
  child_ws = (getDesc()->getBasis()==0) ? 0 :
    new CoinWarmStartBasis(*getDesc()->getBasis());

  // create new node descriptions
  DcoNodeDesc * down_node = new DcoNodeDesc(model);
  down_node->setBroker(broker_);
  DcoNodeDesc * up_node = new DcoNodeDesc(model);
  up_node->setBroker(broker_);
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
    // todo(aykut) where does colLB and colUB get updated?
    // I think they should stay as they created.
    //double lb = model->colLB()[branch_var];
    //double ub = model->colUB()[branch_var];
    double lb = model->getVariables()[branch_var]->getLbHard();
    double ub = model->getVariables()[branch_var]->getUbHard();
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
#if defined(__OA__) || defined(__COLA__)
  down_node->setBasis(child_ws);
#endif
  // Up Node
  // == set other relevant fields of up node
  up_node->setBranchedDir(DcoNodeBranchDirectionUp);
  up_node->setBranchedInd(branch_object->index());
  up_node->setBranchedVal(branch_value);
  // == set warm start basis for the up node.
#if defined(__OA__) || defined(__COLA__)
  up_node->setBasis(child_ws);
#endif
  // Alps does this. We do not need to change the status here
  //status_ = AlpsNodeStatusBranched;

  // push the down and up nodes.
  res.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc*>(down_node),
                               AlpsNodeStatusCandidate,
                               quality_));
  res.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc*>(up_node),
                               AlpsNodeStatusCandidate,
                               quality_));
#ifdef DISCO_DEBUG
  // grumpy message
  int num_inf = 0;
  double sum_inf = 0.0;
  for (int i=0; i<model->branchStrategy()->numBranchObjects(); ++i) {
    double value = model->branchStrategy()->branchObjects()[i]->value();
    // get integer tolerance parameter
    double tolerance = model->dcoPar()->entry(DcoParams::integerTol);
    double dist_to_upper = ceil(value) - value;
    double dist_to_lower = value - floor(value);
    // return the minimum of distance to upper or lower
    double infeas;
    if (dist_to_upper>dist_to_lower) {
      infeas = dist_to_lower;
    }
    else {
      infeas = dist_to_upper;
    }
    if (infeas<tolerance) {
      infeas = 0.0;
    }
    else {
      num_inf++;
    }
    sum_inf += infeas;
  }
  message_handler->message(DISCO_GRUMPY_MESSAGE_LONG, *messages)
    << broker()->getProcRank()
    << broker()->timer().getTime()
    << grumpyMessage[DISCO_GRUMPY_BRANCHED]
    << getIndex()
    << getParentIndex()
    << grumpyDirection[dynamic_cast<DcoNodeDesc*>(desc_)->getBranchedDir()]
    << model->objSense()*getQuality()
    << sum_inf
    << num_inf
    << CoinMessageEol;
  // end of grumpy message
#endif

  setStatus(AlpsNodeStatusBranched);

  // todo(aykut) are these should be in alps level?
  //up_node->setSolEstimate(quality_);
  //down_node->setSolEstimate(quality_);
  return res;
}

// Copies node description of this node to given child node.
// New node is explicitly stored in the memory (no differencing).
void DcoTreeNode::copyFullNode(DcoNodeDesc * child_node) const {
  // get model the node belongs
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  // get message handler and messages for loging
  //CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  //CoinMessages * messages = model->dcoMessages_;
  // get description of this node
  DcoNodeDesc * node_desc = getDesc();

  // solver number of columns
  int num_cols = model->solver()->getNumCols();

  // this will keep bound information, lower and upper
  // it is used for both hard and soft bounds.
  // ownership is transferred, no need to free
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

void DcoTreeNode::processSetPregnant() {
  // get warm start basis from solver
  // todo(aykut) This does not help much if the underlying solver is an IPM
  // based solver.
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  CoinWarmStartBasis * ws = NULL;
  ws = dynamic_cast<CoinWarmStartBasis*>(model->solver()->getWarmStart());
  // store basis in the node desciption.
  getDesc()->setBasis(ws);
  // set status pregnant
  setStatus(AlpsNodeStatusPregnant);

#ifdef DISCO_DEBUG
  // grumpy message
  double sum_inf = 0.0;
  int num_inf = 0;
  for (int i=0; i<model->branchStrategy()->numBranchObjects(); ++i) {
    double value = model->branchStrategy()->branchObjects()[i]->value();
    // get integer tolerance parameter
    double tolerance = model->dcoPar()->entry(DcoParams::integerTol);
    double dist_to_upper = ceil(value) - value;
    double dist_to_lower = value - floor(value);
    // return the minimum of distance to upper or lower
    double infeas;
    if (dist_to_upper>dist_to_lower) {
      infeas = dist_to_lower;
    }
    else {
      infeas = dist_to_upper;
    }
    if (infeas<tolerance) {
      infeas = 0.0;
    }
    else {
      num_inf++;
    }
    sum_inf += infeas;
  }
  model->dcoMessageHandler_->message(DISCO_GRUMPY_MESSAGE_LONG, *model->dcoMessages_)
    << broker()->getProcRank()
    << broker()->timer().getTime()
    << grumpyMessage[DISCO_GRUMPY_PREGNANT]
    << getIndex()
    << getParentIndex()
    << grumpyDirection[dynamic_cast<DcoNodeDesc*>(desc_)->getBranchedDir()]
    << model->objSense()*getQuality()
    << sum_inf
    << num_inf
    << CoinMessageEol;
  // end of grumpy message
#endif
  if (ws) delete ws;
}

// generate constraints is a bitmask for the type of the constraint to
// be generated.
// 01 -> generate OA
// 10 -> generate MILP
// 11 -> generate both OA and MILP
//
//
void DcoTreeNode::branchConstrainOrPrice(BcpsSubproblemStatus subproblem_status,
                                         bool & keepBounding,
                                         bool & branch,
                                         int & generateConstraints,
                                         bool & generateVariables) {
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  // check whether the subproblem solved properly, solver status should be
  // infeasible or optimal.
  if (subproblem_status==BcpsSubproblemStatusPrimalInfeasible) {
    // subproblem is infeasible, fathom
#ifdef DISCO_DEBUG
    // debug message
    message_handler->message(DISCO_SOLVER_INFEASIBLE, *messages)
      << broker()->getProcRank()
      << getIndex()
      << CoinMessageEol;
    keepBounding = false;
    branch = false;
    generateConstraints = 0;
    generateVariables = false;
    // grumpy message
    message_handler->message(DISCO_GRUMPY_MESSAGE_SHORT,
                             *model->dcoMessages_)
      << broker()->getProcRank()
      << broker()->timer().getTime()
      << grumpyMessage[DISCO_GRUMPY_INFEASIBLE]
      << getIndex()
      << getParentIndex()
      << grumpyDirection[dynamic_cast<DcoNodeDesc*>(desc_)->getBranchedDir()]
      << CoinMessageEol;
#endif

    // end of grumpy message
    setStatus(AlpsNodeStatusFathomed);
    return;
  }
  if (subproblem_status!=BcpsSubproblemStatusOptimal and
      subproblem_status!=BcpsSubproblemStatusDualInfeasible) {
    // write the failed subproblem to mps file
    //model->solver()->writeMps("failed", "mps", 0.0);
    // log error message
    message_handler->message(DISCO_SOLVER_FAILED, *messages)
      << broker()->getProcRank()
      << getIndex()
      << CoinMessageEol;
  }
  if (subproblem_status==BcpsSubproblemStatusDualInfeasible) {
    // unbounded problem add cuts
    keepBounding = true;
    branch = false;
    generateVariables = false;
    generateConstraints = 1;
    return;
  }
  // subproblem is solved to optimality. Check feasibility of the solution.
  int numColsInf;
  int numRowsInf;
  double colInf;
  double rowInf;
  DcoSolution * sol = model->feasibleSolution(numColsInf, colInf,
                                              numRowsInf, rowInf);
  // Summary of following if else chain
  // if (both relaxed cols and rows are infeasible) {
  //   decide what to do
  // }
  // else if (relaxed cols are infeasible) {
  //   decide what to do
  // }
  // else if (relaxed rows are infeasible) {
  //   decide what to do
  // }
  // else {
  //   feasible solution found, store it. Fathom this node.
  // }
  if (numColsInf && numRowsInf) {
    // both relaxed columns and relaxed rows are infeasible
    //double cone_tol = model->dcoPar()->entry(DcoParams::coneTol);
    double gap;
    if (quality_<=-ALPS_OBJ_MAX ||
        !broker()->hasKnowledge(AlpsKnowledgeTypeSolution)) {
      // subproblem is not solved yet, or
      // no incumbent solution in broker
      gap = COIN_DBL_MAX;
    }
    else if (broker()->getIncumbentValue()<quality_) {
      // this node should be fathomed
      keepBounding = false;
      branch = false;
      generateConstraints = 0;
      generateVariables = false;
      return;
    }
    else {
      gap = (broker()->getIncumbentValue() - quality_)/fabs(quality_);
    }

    // Apply OA cut strategy, apply milp once OA is done.

    // OA cut generation strategy:
    // generate OA cuts for alpha many times after
    // alpha iterations generate OA cuts only if the last improvement was
    // larger than beta times current gap (i.e., we think we can get bound
    // larger than current incumbent and fathom, so we invest) and total number
    // of iteratios is less than gamma (we invested gamma iterations and could
    // not fathom, so stop).

    // get parameters
    int alpha = model->dcoPar()->entry(DcoParams::cutOaAlpha);
    double beta = model->dcoPar()->entry(DcoParams::cutOaBeta);
    int gamma = model->dcoPar()->entry(DcoParams::cutOaGamma);
    if (((bcpStats_.numBoundIter_<alpha)
         or (bcpStats_.lastImp_>beta*gap)) and
        (bcpStats_.numBoundIter_<gamma)) {
      // std::cout << "node " << index_
      //           << " iter " << bcpStats_.numBoundIter_
      //           << " total cuts "<< bcpStats_.numTotalCuts_
      //           << " last cuts "<< bcpStats_.numLastCuts_
      //           << " total imp "<< bcpStats_.totalImp_
      //           << " last imp "<< bcpStats_.lastImp_
      //           << " gap " << gap
      //           << std::endl;
      // cut
      keepBounding = true;
      branch = false;
      generateVariables = false;
      generateConstraints = 1;
    }
    else {
      // decided not to generate OA cuts, check MILP
      // invoke checking model->cutGenerationFrequency_
      int freq = model->cutGenerationFrequency();
      double delta = model->dcoPar()->entry(DcoParams::cutMilpDelta);
      if (bcpStats_.numMilpIter_ == 0 or
          (model->numConGenerators() > 1 and index_ % freq == 0
           and bcpStats_.lastImp_/bcpStats_.lastObjVal_ > delta)) {
        // generate MILP
        keepBounding = true;
        branch = false;
        generateVariables = false;
        generateConstraints = 2;
        bcpStats_.numMilpIter_++;
      }
      else {
        // no OA no MILP cuts
        // branch
        keepBounding = false;
        branch = true;
        generateVariables = false;
        generateConstraints = 0;
        bcpStats_.numMilpIter_ = 0;
      }
    }
  }
  else if (numColsInf) {
    // only relaxed columns are infeasible, relaxed rows are feasible
    // generate milp cuts or branch
    int freq = model->cutGenerationFrequency();
    double delta = model->dcoPar()->entry(DcoParams::cutMilpDelta);
    if (bcpStats_.numMilpIter_ == 0 or
        (model->numConGenerators() > 1 and index_%freq == 0
         and bcpStats_.lastImp_/bcpStats_.lastObjVal_ > delta)) {
      // generate MILP
      keepBounding = true;
      branch = false;
      generateVariables = false;
      generateConstraints = 2;
      bcpStats_.numMilpIter_++;
    }
    else {
      // no OA no MILP cuts
      // branch
      keepBounding = false;
      branch = true;
      generateVariables = false;
      generateConstraints = 0;
      bcpStats_.numMilpIter_ = 0;
    }
  }
  else if (numRowsInf) {
    // all relaxed cols are feasbile, only relaxed rows are infeasible
    // generate cuts using OA.
    keepBounding = true;
    branch = false;
    generateVariables = false;
    generateConstraints = 1;
  }
  else if (sol) {
    // all relaxed cols and rows are feasible
    sol->setDepth(depth_);
    //todo(aykut) solution has an index field I beleive that should be filled
    //in Alps level.
    sol->setIndex(broker()->getNumKnowledges(AlpsKnowledgeTypeSolution));
    // Store in Alps pool
    broker()->addKnowledge(AlpsKnowledgeTypeSolution,
                           sol,
                           model->objSense() * sol->getQuality());
    // update solver with the new incumbent if better
    double incum_val = broker()->getIncumbentValue();
    model->solver()->setDblParam(OsiDualObjectiveLimit,
                                 model->objSense()*incum_val);

#ifdef DISCO_DEBUG
    // set node status to fathomed.
    message_handler->message(0, "Dco", "Node is feasible, fathoming... ",
                             'G', DISCO_DLOG_BRANCH)
      << CoinMessageEol;
    // grumpy message
    model->dcoMessageHandler_->message(DISCO_GRUMPY_MESSAGE_MED,
                                       *model->dcoMessages_)
      << broker()->getProcRank()
      << broker()->timer().getTime()
      << grumpyMessage[DISCO_GRUMPY_INTEGER]
      << getIndex()
      << getParentIndex()
      << grumpyDirection[dynamic_cast<DcoNodeDesc*>(desc_)->getBranchedDir()]
      << model->objSense()*getQuality()
      << CoinMessageEol;
    // end of grumpy message
#endif

    // stop bounding, branching, cutting.
    keepBounding = false;
    branch = false;
    generateVariables = false;
    generateConstraints = 0;

    // set status to fathom
    setStatus(AlpsNodeStatusFathomed);
  }
  else {
    message_handler->message(DISCO_SHOULD_NOT_HAPPEN, *messages)
      << broker()->getProcRank()
      << getIndex()
      << CoinMessageEol;
  }
}

//todo(aykut) replace this with DcoModel::feasibleSolution????
void DcoTreeNode::checkRelaxedCols(int & numInf) {
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
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

// todo(aykut) this should go into the con generator. process the cuts given
// by cgl in disco con generator.
void DcoTreeNode::applyConstraints(BcpsConstraintPool const * conPool) {
  DcoModel * model = dynamic_cast<DcoModel*>(broker()->getModel());
  double scale_par = model->dcoPar()->entry(DcoParams::scaleConFactor);
  double density_par = model->dcoPar()->entry(DcoParams::denseConFactor);
  //double tailoff = model->dcoPar()->entry(DcoParams::tailOff);
  double cone_tol = model->dcoPar()->entry(DcoParams::coneTol);
  double const * sol = model->solver()->getColSolution();
  int num_cols = model->solver()->getNumCols();

  // Tranform constraints to Osi cut.
  int num_cuts = conPool->getNumConstraints();
  OsiRowCut const ** cuts_to_add = new OsiRowCut const * [num_cuts];
  int num_add = 0;
  std::vector<int> cuts_to_del;

  // iterate over cuts and
  //------------------------------------------
  // Remove following MILP cuts:
  //  (1) empty cuts
  //  (2) dense cuts
  //  (3) bad scaled cuts
  //  (4) weak cuts
  //  (5) almost parallel cuts
  //------------------------------------------
  for (int i=0; i<num_cuts; ++i) {
    DcoLinearConstraint * curr_con =
      dynamic_cast<DcoLinearConstraint*>(conPool->getConstraint(i));
    int length = curr_con->getSize();
    double const * elements = curr_con->getValues();
    int const * indices = curr_con->getIndices();
    double curr_con_lb = curr_con->getLbSoft();
    double curr_con_ub = curr_con->getUbSoft();

    // add all OA cuts
    if (curr_con->constraintType() == DcoConstraintTypeOA) {
      cuts_to_add[num_add++] = curr_con->createOsiRowCut(model);
      // update cut statistics
      model->conGenerators(curr_con->constraintType())->stats().addNumConsUsed(1);
      continue;
    }

    // (1) check whether cut is empty.
    if (length <= 0) {
      cuts_to_del.push_back(i);
      continue;
    }

    // (2) check whether cut is dense.
    if((int) density_par * length > num_cols) {
      // discard the cut
      cuts_to_del.push_back(i);
      continue;
    }

    // (3) check cut scaling
    double activity = 0.0;
    double maxElem = 0.0;
    double minElem = ALPS_DBL_MAX;
    double scaleFactor;

    // iterate over coefficients
    for (int k=0; k<length; ++k) {
      if (elements[k]==0.0) {
        // if a coef is exactly 0, ignore it.
        // what if all coef is zero?
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

    // (4) Check whether the cut is weak.
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
    if (violation < cone_tol) {
      // cut is weak, skip it.
      cuts_to_del.push_back(i);

#ifdef DISCO_DEBUG
      CoinMessageHandler * message_handler = model->dcoMessageHandler_;
      CoinMessages * messages = model->dcoMessages_;
      message_handler->message(DISCO_INEFFECTIVE_CUT, *messages)
        << broker()->getProcRank()
        << getIndex()
        << CoinMessageEol;
#endif

      continue;
    }

    // (5) Check almost parallel cuts
    // create dense cut
    double * dense_cut = new double[num_cols]();
    double cut_norm = 0.0;
    for (int i = 0; i < length; ++i) {
      dense_cut[indices[i]] = elements[i];
      cut_norm += elements[i]*elements[i];
    }
    cut_norm = sqrt(cut_norm);

    // {
    //   // check which one cuts the solution for problem r22c30k10i20
    //   double sol[] = {7.0, 2.0, 6.25, 8.0, 7.0, 3.85714285091261, 6.0,
    //                   3.0, 3.85714279912700, 8.0, 3.0, 7.4, 12, 5,
    //                   10.75, 8, 4, 6, 12, 3, 9, 9,
    //                   6, 5, 10, 4, 8.2, 7, 6, 2.4, 0.0};
    //   // check whether cut_i cuts the sol
    //   double activity = std::inner_product(sol, sol+30, dense_cut, 0.0);
    //   std::cout << curr_con_lb << " " << activity << " " << curr_con_ub << std::endl;
    //   if (activity < curr_con_lb or activity > curr_con_ub) {
    //     std::cout << "=================================================="
    //               << std::endl;
    //     std::cout << "size: " <<  curr_con->getSize() << std::endl;
    //     std::cout << "ind" << std::endl;
    //     for (int ii=0; ii<curr_con->getSize(); ++ii) {
    //       std::cout << curr_con->getIndices()[ii] << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "val" << std::endl;
    //     for (int ii=0; ii<curr_con->getSize(); ++ii) {
    //       std::cout << curr_con->getValues()[ii] << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "lb " << curr_con_lb << std::endl;
    //     std::cout << "ub " << curr_con_ub << std::endl;

    //     std::cout << "=================================================="
    //               << std::endl;
    //     //throw std::exception();
    //   }
    // }

    bool added = false;
    // Check whether cut is parallel to an existing constraint or cut.
    {
      // cut is stored at dense_cut
      // compute par := max _i { cut^T a_i / |cut| |a_i| } where a_i is
      // constraint. cut is discarded if par > 0.98
      //double par_factor = 0.0;
      // iterate over constraints
      CoinPackedMatrix const * mat = model->solver()->getMatrixByRow();
      double const * lb = model->solver()->getRowLower();
      double const * ub = model->solver()->getRowUpper();
      int num_rows = mat->getNumRows();
      // iterate over rows
      for (int k = 0; k < num_rows; ++k) {
        // get row
        //int row_size = mat->getVectorStarts()[k+1] - mat->getVectorStarts()[k];
        int row_size = mat->getVectorLast(k) - mat->getVectorFirst(k);
        int const * ind = mat->getIndices() + mat->getVectorFirst(k);
        double const * val = mat->getElements() + mat->getVectorFirst(k);
        double inn_prod = 0.0;
        double row_norm = 0.0;

        if (length != row_size) {
          // cut and rows are of different size
          // check next row
          continue;
        }

        for (int j = 0; j < row_size; ++j) {
          inn_prod += val[j]*dense_cut[ind[j]];
          row_norm += val[j]*val[j];
        }
        row_norm = sqrt(row_norm);
        // divide by norms
        inn_prod = inn_prod/(cut_norm*row_norm);
        // todo(aykut) parametrize the following.
        if (inn_prod > 0.99999) {
          // cut coeff are same, update bound if tighter

          // scale cut lower and upper bounds
          double scale = 0.0;
          scale = val[0]/dense_cut[ind[0]];
          if (curr_con_lb > -1e8 and lb[k] < scale * curr_con_lb) {
            if (index_ == 0) {
              added = true;
              model->solver()->setRowLower(k, scale*curr_con_lb-0.000001);
              break;
            }
            else {
              // keep the cut
              inn_prod = 0.0;
            }
          }
          if (curr_con_ub < 1e8 and ub[k] > scale * curr_con_ub) {
            // we can do this only at root, else just keep it as a cut
            if (index_ == 0) {
              added = true;
              model->solver()->setRowUpper(k, scale*curr_con_ub+0.000001);
              break;
            }
            else {
              // keep the cut
              inn_prod = 0.0;
            }
          }
        }
        if (inn_prod > 0.95 and !cuts_to_del.empty() and i != cuts_to_del.back()) {
          // almost parallel
          cuts_to_del.push_back(i);
        }
        if (!cuts_to_del.empty() and i == cuts_to_del.back()) {
          // already decided to delete the cut, no need to check other
          // constraints/cuts currently in the solver
          break;
        }
      }
    }
    // update cut statistics
    model->conGenerators(curr_con->constraintType())->stats().addNumConsUsed(1);
    if (not added) {
      cuts_to_add[num_add++] = curr_con->createOsiRowCut(model);
    }
    delete[] dense_cut;
  }

  // Add cuts to lp and adjust basis.
  CoinWarmStartBasis * ws = dynamic_cast<CoinWarmStartBasis*>
    (model->solver()->getWarmStart());

  int num_solver_rows = model->solver()->getNumRows();
  int num_solver_cols = model->solver()->getNumCols();

  // update statistics of the cut generator
  if (num_add > 0) {
    // update cut statistics
    bcpStats_.numTotalCuts_ += num_add;
    bcpStats_.numLastCuts_ = num_add;
    // apply cuts
    model->solver()->applyRowCuts(num_add, cuts_to_add);
    ws->resize(num_solver_rows + num_add, num_solver_cols);
    for (int i=0; i<num_add; ++i) {
      ws->setArtifStatus(num_solver_rows + i,
                         CoinWarmStartBasis::basic);
    }
    if (model->solver()->setWarmStart(ws) == false) {
#ifdef DISCO_DEBUG
      CoinMessageHandler * message_handler = model->dcoMessageHandler_;
      CoinMessages * messages = model->dcoMessages_;
      message_handler->message(DISCO_FAILED_WARM_START, *messages)
        << broker()->getProcRank()
        << getIndex()
        << CoinMessageEol;
#endif
    }
    for (int k=0; k<num_add; ++k) {
      delete cuts_to_add[k];
    }
  }

#ifdef DISCO_DEBUG
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  message_handler->message(DISCO_CUTS_ADDED, *messages)
    << broker()->getProcRank()
    << getIndex()
    << num_add
    << conPool->getNumConstraints()
    << CoinMessageEol;
#endif

  delete[] cuts_to_add;
  if (ws) delete ws;
}

/// Pack this into an encoded object.
AlpsReturnStatus DcoTreeNode::encode(AlpsEncoded * encoded) const {
  // get pointers for message logging
  assert(broker_);

  // return value
  AlpsReturnStatus status;
  status = AlpsTreeNode::encode(encoded);
  assert(status==AlpsReturnStatusOk);
  status = BcpsTreeNode::encode(encoded);
  assert(status==AlpsReturnStatusOk);

#ifdef DISCO_DEBUG
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  message_handler->message(DISCO_NODE_ENCODED, *messages)
    << broker()->getProcRank()
    << getIndex()
    << CoinMessageEol;
#endif

  return status;
}

/// Unpack into a new DcoTreeNode object and return a pointer to it.
AlpsKnowledge * DcoTreeNode::decode(AlpsEncoded & encoded) const {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  // return value
  AlpsReturnStatus status;
  // todo(aykut) we are decoing, how do we know the model to assign the new
  // node is same as the model_ of this? Is this due to fact that they are in
  // the same processor?
  AlpsNodeDesc * new_node_desc = new DcoNodeDesc();
  DcoTreeNode * new_node = new DcoTreeNode(new_node_desc);
  new_node->setBroker(broker_);
  new_node_desc = NULL;
  status = new_node->decodeToSelf(encoded);
  if (status != AlpsReturnStatusOk) {
    message_handler->message(DISCO_UNEXPECTED_DECODE_STATUS,
                             *messages)
      << __FILE__ << __LINE__ << CoinMessageEol;
  }
  return new_node;
}

/// Unpack into this from an encoded object.
AlpsReturnStatus DcoTreeNode::decodeToSelf(AlpsEncoded & encoded) {
  // get pointers related to message logging.
  assert(broker_);
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;

  // return value
  AlpsReturnStatus status;
  status = AlpsTreeNode::decodeToSelf(encoded);
  if (status != AlpsReturnStatusOk) {
    message_handler->message(DISCO_UNEXPECTED_DECODE_STATUS,
                             *messages)
      << __FILE__ << __LINE__ << CoinMessageEol;
  }
  status = BcpsTreeNode::decodeToSelf(encoded);
  if (status != AlpsReturnStatusOk) {
    message_handler->message(DISCO_UNEXPECTED_DECODE_STATUS,
                             *messages)
      << __FILE__ << __LINE__ << CoinMessageEol;
  }

#ifdef DISCO_DEBUG
  message_handler->message(DISCO_NODE_DECODED, *messages)
    << broker()->getProcRank()
    << getIndex()
    << CoinMessageEol;
#endif

  if (isPregnant()) {
    // delete branch object stored, we will re-create it. Ideally it should not
    // be sent in the first place.
    clearBranchObject();
    setStatus(AlpsNodeStatusEvaluated);
  }
  // initiate bcpStats_
  bcpStats_.numBoundIter_ = 0;
  bcpStats_.numTotalCuts_ = 0;
  bcpStats_.numLastCuts_ = 0;
  bcpStats_.numOaCuts_ = 0;
  bcpStats_.numMilpCuts_ = 0;
  bcpStats_.totalImp_ = 0.0;
  bcpStats_.lastImp_ = 0.0;
  bcpStats_.lastObjVal_ = 0.0;
  bcpStats_.startObjVal_ = 0.0;
  bcpStats_.numMilpIter_ = 0;
  bcpStats_.inactive_.clear();
  bcpStats_.generatorIndex_.clear();
  return status;
}
