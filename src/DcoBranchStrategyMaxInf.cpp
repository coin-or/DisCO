#include "DcoBranchStrategyMaxInf.hpp"
#include "DcoModel.hpp"

DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf() {
}

DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf(DcoBranchStrategyMaxInf const & other) {
}

DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf(DcoModel * model):
  BcpsBranchStrategy(model) {
}

DcoBranchStrategyMaxInf::~DcoBranchStrategyMaxInf() {
}

BcpsBranchStrategy * DcoBranchStrategyMaxInf::clone() const {
  return new DcoBranchStrategyMaxInf(*this);
}

int DcoBranchStrategyMaxInf::createCandBranchObjects(int numPassesLeft,
                                                     double ub) {
  // get model
  DcoModel * model = dynamic_cast<DcoModel*>(model_);
  // get integer tolerance parameter
  double tolerance = model->dcoPar()->entry(DcoParams::integerTol);
  // get number of relaxed columns
  // we assume all relaxed columns are integer variables.
  int num_relaxed = model->numRelaxedCols();
  // get indices of relaxed object
  int const * relaxed = model->relaxedCols();
  // store branch objects in bobjects
  std::vector<BcpsBranchObject*> bobjects;
  // iterate over relaxed columns and populate bobjects
  for (int i=0; i<num_relaxed; ++i) {
    int preferredDir;
    BcpsObject * curr_object = model->getVariables()[relaxed[i]];
    double infeasibility = curr_object->infeasibility(model, preferredDir);
    // check the amount of infeasibility
    if (infeasibility != 0.0) {
      // create a branch object for this
      BcpsBranchObject * cb =
        curr_object->createBranchObject(model, preferredDir);
      bobjects.push_back(cb);
    }
  }
  // set number of branch objects
  //setNumBranchObjects(bobjects.size());
  // add branch objects to branchObjects_
  setBranchObjects(bobjects);
  // bobjects are now owned by BcpsBranchStrategy, do not free them.
  std::vector<BcpsBranchObject*>::iterator it;
  for (it=bobjects.begin(); it!=bobjects.end(); ++it) {
    *it = NULL;
  }
  bobjects.clear();
  // compare branch objects and keep the best one at bestBranchObject_
  return 0;
}

/** Compare branching object thisOne to bestSoFar. If thisOne is better than
    bestObject, return branching direction(1 or -1), otherwise return 0.  If
    bestSorFar is NULL, then always return branching direction(1 or -1).

    todo(aykut): this does not make sense to me. I think the return value
    can be made to represent something simpler. Now it represents a direction
    or update of the best so far.
*/
int DcoBranchStrategyMaxInf::betterBranchObject(BcpsBranchObject * b,
                                                BcpsBranchObject * bestSoFar) {
  // get model
  DcoModel * model = dynamic_cast<DcoModel*>(model_);
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  if (bestSoFar==NULL) {
    // best is NULL, set it to b.
    bestSoFar = b;
    return b->getDirection();
  }
  else {
    // compare b with bestSoFar
    if (b==NULL) {
      message_handler->message(9998, "Dco", "Can not compare, "
                               "NULL branch object is given. ",
                               'E', 0)
        << CoinMessageEol;
    }
    else {
      double value = b->getValue();
      // distance to upper
      double dt_upper = ceil(value) - value;
      // distance to lower
      double dt_lower = value - floor(value);
      // fractional distance value
      double d = (dt_upper<dt_lower) ? dt_upper : dt_lower;
      // branched value for the best
      double best_value = bestSoFar->getValue();
      double best_dt_upper = ceil(best_value) - best_value;
      double best_dt_lower = best_value - floor(best_value);
      double best_d = (best_dt_upper<best_dt_lower) ? best_dt_upper : best_dt_lower;
      if (d>best_d) {
        // update best
        bestSoFar = b;
        return b->getDirection();
      }
      else {
        // do not update best
        return 0;
      }
    }
  }
}
