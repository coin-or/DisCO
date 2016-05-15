#include "DcoBranchObject.hpp"

DcoBranchObject::DcoBranchObject(DcoModel * model, int colInd, int intScore,
                                 double dblScore, int direction, double value):
  BcpsBranchObject(model, colInd, intScore, dblScore, direction, value) {
  type_ = DcoBranchingObjectTypeInt;
  ubDownBranch_ = floor(value);
  lbUpBranch_ = ceil(value);
}


DcoBranchObject::~DcoBranchObject() {
}

BcpsBranchObject * DcoBranchObject::clone() const {
  DcoModel * dco_model = dynamic_cast<DcoModel*> (model());
  int index = getObjectIndex();
  double up_score = getUpScore();
  double down_score = getDownScore();
  int direction = getDirection();
  double value = getValue();
  BcpsBranchObject * nbo = new DcoBranchObject(dco_model, index, up_score,
                                               down_score, direction, value);
  return nbo;
}

double DcoBranchObject::branch(bool normalBranch) {
  DcoModel * model = dynamic_cast<DcoModel*> (model_);
  int col_index = model->relaxedCols()[objectIndex_];
  // Decrement number of branches left by 1.
  --numBranchesLeft_;
  if (direction_<0) {
    model->solver()->setColUpper(col_index, ubDownBranch_);
    direction_ = 1;
  }
  else {
    model->solver()->setColLower(col_index, lbUpBranch_);
    direction_ = -1;          // Swap direction
  }
  return 0.0;
}

bool DcoBranchObject::boundBranch() const {
  // return True if branching fixes the bound.
  // return True if variable we are branching is binary
  return false;
}
