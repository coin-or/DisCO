#include "DcoBranchStrategyStrong.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoTreeNode.hpp"


DcoBranchStrategyStrong::DcoBranchStrategyStrong(DcoModel * model)
  : BcpsBranchStrategy(model) {
  setType(static_cast<int>(DcoBranchingStrategyStrong));
}

// Assumes problem is not unbounded.
void DcoBranchStrategyStrong::updateScore(BcpsBranchObject * bobject,
                                            double orig_lb,
                                            double orig_ub,
                                            double orig_obj) const {
  // get dco model and message stuff
  DcoModel * dco_model = dynamic_cast<DcoModel*>(model());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  DcoBranchObject * dco_bobject = dynamic_cast<DcoBranchObject*>(bobject);
  // solve subproblem for the down branch
  dco_model->solver()->setColUpper(bobject->index(), dco_bobject->ubDownBranch());
  dco_model->solver()->solveFromHotStart();
  double down_obj = dco_model->solver()->getObjValue();
  // restore bound
  dco_model->solver()->setColUpper(bobject->index(), orig_ub);


  // solve subproblem for the up branch
  dco_model->solver()->setColLower(bobject->index(), dco_bobject->lbUpBranch());
  dco_model->solver()->solveFromHotStart();
  double up_obj = dco_model->solver()->getObjValue();
  // restore bound
  dco_model->solver()->setColLower(bobject->index(), orig_lb);
  // set score
  double down_diff = fabs(orig_obj-down_obj);
  double up_diff = fabs(orig_obj-up_obj);
  double score = down_diff>up_diff ? down_diff : up_diff;
  bobject->setScore(score);
}


int DcoBranchStrategyStrong::createCandBranchObjects(BcpsTreeNode * node) {
  // notes(aykut)
  // Considers all set of integers for now?
  // What happens in case of IPM solvers?
  // consider time limit
  // iter limit

  // get node
  DcoTreeNode * dco_node = dynamic_cast<DcoTreeNode*>(node);
  // get dco model and message stuff
  DcoModel * dco_model = dynamic_cast<DcoModel*>(model());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  // get number of relaxed columns
  // we assume all relaxed columns are integer variables.
  int num_relaxed = dco_model->numRelaxedCols();
  // get indices of relaxed object
  int const * relaxed = dco_model->relaxedCols();

  // create data to keep branching objects
  int cand_cap = dco_model->dcoPar()->entry(DcoParams::strongCandSize);
  cand_cap = CoinMax(CoinMin(cand_cap, num_relaxed), 1);
  BcpsBranchObject ** bobjects = new BcpsBranchObject*[cand_cap];
  int num_bobjects = 0;
  // pos is the index of the minimum score branch object
  int min_pos = -1;
  double min_score = ALPS_INFINITY;

  dco_model->solver()->markHotStart();
  dco_model->solver()->setIntParam(OsiMaxNumIterationHotStart, 50);

  double const obj_val = dco_model->solver()->getObjValue();
  double const * collb = dco_model->solver()->getColLower();
  double const * colub = dco_model->solver()->getColUpper();

  // iterate over integer cols, create branch objects, solve corresponding
  // relaxed problems and set scores.
  for (int i=0; i<num_relaxed; ++i) {
    // get corresponding variable
    int var_index = relaxed[i];
    BcpsObject * curr_var = dco_model->getVariables()[var_index];
    int preferDir;
    double infeasibility = curr_var->infeasibility(dco_model, preferDir);
    if (!infeasibility) {
      continue;
    }
    BcpsBranchObject * curr_branch_object = curr_var->createBranchObject(dco_model, preferDir);
    updateScore(curr_branch_object, collb[var_index], colub[var_index], obj_val);
    double curr_score = curr_branch_object->score();

    // if we have capacity add branch object
    // else check whether current performs better than the worst candidate
    // if it is add it to candidates.
    if (num_bobjects<cand_cap) {
      bobjects[num_bobjects] = curr_branch_object;
      if (curr_score<min_score) {
        min_score = curr_score;
        min_pos = num_bobjects;
      }
      num_bobjects++;
    }
    else if (curr_score>min_score) {
      delete bobjects[min_pos];
      bobjects[min_pos] = curr_branch_object;
      // find new minimum score candidate
      min_score = ALPS_INFINITY;
      for (int k=0; k<cand_cap; ++k) {
        if (bobjects[k]->score()<min_score) {
          min_score = bobjects[k]->score();
          min_pos = k;
        }
      }
    }
    else {
      // score is not enough to be a candidate
    }
  }
  if (num_bobjects==0) {
    std::cout << "All columns are feasible." << std::endl;
    throw std::exception();
  }
  dco_model->solver()->unmarkHotStart();

  // add branch objects to branchObjects_
  setBranchObjects(num_bobjects, bobjects);
  // bobjects are now owned by BcpsBranchStrategy, do not free them.
  //bobjects = NULL;
  // set the branch object member of the node
  dco_node->setBranchObject(new DcoBranchObject(bestBranchObject()));
  // compare branch objects and keep the best one at bestBranchObject_
  return 0;
}



int
DcoBranchStrategyStrong::betterBranchObject(BcpsBranchObject const * current,
                                            BcpsBranchObject const * other) {
  // get model
  DcoModel * dco_model = dynamic_cast<DcoModel*>(model());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  int res;
  if (current->score()>other->score()) {
    res = 1;
  }
  else {
    res = 0;
  }
  return res;
}
