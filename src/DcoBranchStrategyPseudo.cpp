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


#include "DcoBranchStrategyPseudo.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoTreeNode.hpp"
#include "DcoBranchObject.hpp"

DcoBranchStrategyPseudo::DcoBranchStrategyPseudo(DcoModel * model):
  BcpsBranchStrategy(model) {
  setType(DcoBranchingStrategyPseudoCost);
  // todo(aykut) think about parametrizing this.
  score_factor_ = 1.0/6.0;
  int num_relaxed = model->numRelaxedCols();
  down_num_ = new int[num_relaxed]();
  up_num_ = new int[num_relaxed]();
  down_derivative_ = new double[num_relaxed]();
  up_derivative_ = new double[num_relaxed]();
  // fill reverse map
  int const * relaxed_cols = model->relaxedCols();
  for (int i=0; i<num_relaxed; ++i) {
    rev_relaxed_[relaxed_cols[i]] = i;
  }
}

DcoBranchStrategyPseudo::~DcoBranchStrategyPseudo() {
  if (down_num_) {
    delete[] down_num_;
    down_num_ = NULL;
  }
  if (up_num_) {
    delete[] up_num_;
    up_num_ = NULL;
  }
  if (down_derivative_) {
    delete[] down_derivative_;
    down_derivative_ = NULL;
  }
  if (up_derivative_) {
    delete[] up_derivative_;
    up_derivative_ = NULL;
  }
}

int DcoBranchStrategyPseudo::createCandBranchObjects(BcpsTreeNode * node) {
  // get node
  DcoTreeNode * dco_node = dynamic_cast<DcoTreeNode*>(node);
  // update statistics
  update_statistics(dco_node);
  // get dco model and message stuff
  DcoModel * dco_model = dynamic_cast<DcoModel*>(model());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  // get number of relaxed columns
  // we assume all relaxed columns are integer variables.
  int num_relaxed = dco_model->numRelaxedCols();
  // get indices of relaxed object
  int const * relaxed = dco_model->relaxedCols();
  // store branch objects in bobjects
  std::vector<BcpsBranchObject*> bobjects;
  // iterate over relaxed columns and populate bobjects
  for (int i=0; i<num_relaxed; ++i) {
    int preferredDir;
    BcpsObject * curr_object = dco_model->getVariables()[relaxed[i]];
    double infeasibility = curr_object->infeasibility(dco_model, preferredDir);
    // check the amount of infeasibility
    if (infeasibility != 0.0) {
      double min = std::min(down_derivative_[i], up_derivative_[i]);
      double max = std::max(down_derivative_[i], up_derivative_[i]);
      // compute score
      double score = score_factor_*max + (1.0-score_factor_)*min;
      // create a branch object for this
      BcpsBranchObject * cb =
        curr_object->createBranchObject(dco_model, preferredDir);
      // set score
      cb->setScore(score);
      bobjects.push_back(cb);

      // debug stuff
      message_handler->message(DISCO_PSEUDO_REPORT, *messages)
        << dco_model->broker()->getProcRank()
        << relaxed[i]
        << score
        << CoinMessageEol;
    }
  }
  // add branch objects to branchObjects_
  setBranchObjects(bobjects);
  // bobjects are now owned by BcpsBranchStrategy, do not free them.
  bobjects.clear();
  // set the branch object member of the node
  dco_node->setBranchObject(new DcoBranchObject(bestBranchObject()));
  // compare branch objects and keep the best one at bestBranchObject_
  return 0;
}

int
DcoBranchStrategyPseudo::betterBranchObject(BcpsBranchObject const * current,
                                            BcpsBranchObject const * other) {
  // get model
  // DcoModel * dco_model = dynamic_cast<DcoModel*>(model());
  // CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  // CoinMessages * messages = dco_model->dcoMessages_;
  int res;
  if (current->score()>other->score()) {
    res = 1;
  }
  else {
    res = 0;
  }
  return res;
}

void DcoBranchStrategyPseudo::update_statistics(DcoTreeNode * node) {
  // return if this is the root node
  if (node->getParent()==NULL) {
    return;
  }
  // get dco model and message stuff
  DcoModel * dco_model = dynamic_cast<DcoModel*>(model());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  // get quality_ of this node, quality is sense*value
  double quality = node->getQuality();
  // get quality_ of the parent node
  double parent_quality = node->getParent()->getQuality();
  // is this node a down or up branch
  int dir = node->getDesc()->getBranchedDir();
  // index of the branched variable for the current node
  int branched_index = rev_relaxed_[node->getDesc()->getBranchedInd()];
  double branched_value = node->getDesc()->getBranchedVal();

  // update statistics
  double frac;
  if (dir==DcoNodeBranchDirectionDown) {
    frac = branched_value-floor(branched_value);
    double deriv = (quality-parent_quality) / frac;
    int n = down_num_[branched_index];
    double old = down_derivative_[branched_index];
    down_derivative_[branched_index] = (old*n + deriv)/(n+1);
    down_num_[branched_index]++;

    // debug stuff
    message_handler->message(DISCO_PSEUDO_DUP, *messages)
      << dco_model->broker()->getProcRank()
      << node->getDesc()->getBranchedInd()
      << old
      << down_derivative_[branched_index]
      << frac
      << CoinMessageEol;
  }
  else if (dir==DcoNodeBranchDirectionUp) {
    frac = ceil(branched_value)-branched_value;
    double deriv = (quality-parent_quality) / frac;
    int n = up_num_[branched_index];
    double old = up_derivative_[branched_index];
    up_derivative_[branched_index] = (old*n + deriv)/(n+1);
    up_num_[branched_index]++;

    // debug stuff
    message_handler->message(DISCO_PSEUDO_UUP, *messages)
      << dco_model->broker()->getProcRank()
      << node->getDesc()->getBranchedInd()
      << old
      << up_derivative_[branched_index]
      << frac
      << CoinMessageEol;
  }
  else {
    message_handler->message(9998, "Dco", "Invalid branching direction. ",
                             'E', 0)
      << CoinMessageEol;
  }
}
