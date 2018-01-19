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


#include "DcoBranchStrategyMaxInf.hpp"
#include "DcoModel.hpp"
#include "DcoTreeNode.hpp"
#include "DcoBranchObject.hpp"

DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf(DcoModel * model)
  : BcpsBranchStrategy(model) {
  setType(DcoBranchingStrategyMaxInfeasibility);
}

int DcoBranchStrategyMaxInf::createCandBranchObjects(BcpsTreeNode * node) {
  // get node
  DcoTreeNode * dco_node = dynamic_cast<DcoTreeNode*>(node);
  // get model
  DcoModel * dco_model = dynamic_cast<DcoModel*>(model());
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
      // create a branch object for this
      BcpsBranchObject * cb =
        curr_object->createBranchObject(dco_model, preferredDir);
      // set score
      cb->setScore(infeasibility);
      bobjects.push_back(cb);
    }
  }
  // add branch objects to branchObjects_
  setBranchObjects(bobjects);
  // bobjects are now owned by BcpsBranchStrategy, do not free them.
  std::vector<BcpsBranchObject*>::iterator it;
  bobjects.clear();
  // set the branch object member of the node
  dco_node->setBranchObject(new DcoBranchObject(bestBranchObject()));
  // compare branch objects and keep the best one at bestBranchObject_
  return 0;
}

/// Compare current to other, return 1 if current is better, 0 otherwise
int
DcoBranchStrategyMaxInf::betterBranchObject(BcpsBranchObject const * current,
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
