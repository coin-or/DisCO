#include "DcoNodeDesc.hpp"

DcoNodeDesc::DcoNodeDesc() {
  // set if as left branch by default
  branchedDir_ = DcoNodeBranchDirectionLeft;
  branchedInd_ = -1;
  // set to 0.0 by default.
  branchedVal_ = 0.0;
  basis_ = NULL;
}

DcoNodeDesc::DcoNodeDesc(DcoModel * model): BcpsNodeDesc(model) {
  // set if as left branch by default
  branchedDir_ = DcoNodeBranchDirectionLeft;
  branchedInd_ = -1;
  // set to 0.0 by default.
  branchedVal_ = 0.0;
  basis_ = NULL;
}

DcoNodeDesc::~DcoNodeDesc() {
  if (basis_) {
    delete[] basis_;
  }
}

DcoNodeBranchDir DcoNodeDesc::branchedDir() const {
  return branchedDir_;
}

int DcoNodeDesc::branchedInd() const {
  return branchedInd_;
}

double DcoNodeDesc::branchedVal() const {
  return branchedVal_;
}

// set fields
void DcoNodeDesc::setBranchedDir(DcoNodeBranchDir dir) {
  branchedDir_ = dir;
}

void DcoNodeDesc::setBranchedInd(int ind) {
  branchedInd_ = ind;
}

void DcoNodeDesc::setBranchedVal(double val) {
  branchedVal_ = val;
}

void DcoNodeDesc::setBasis(CoinWarmStartBasis *& ws) {
  if (basis_) {
    delete basis_;
  }
  basis_= ws;
  ws = NULL;
}

/** Get warm start basis. */
CoinWarmStartBasis * DcoNodeDesc::getBasis() const {
  return basis_;
}
