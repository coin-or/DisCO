#include "DcoNodeDesc.hpp"

DcoNodeDesc::DcoNodeDesc() {
}

DcoNodeDesc::DcoNodeDesc(DcoModel * model): BcpsNodeDesc(model) {
}

DcoNodeDesc::~DcoNodeDesc() {
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
