#include "DcoConstraint.hpp"

#include <OsiRowCut.hpp>
#include <CoinHelperFunctions.hpp>

DcoConstraint::DcoConstraint(): cType_(DCO_LINEAR) {
}

DcoConstraint::DcoConstraint(double lbh, double ubh, double lbs,
			     double ubs, DcoConstraintType type):
  BcpsConstraint(lbh, ubh, lbs, ubs), cType_(type) {
}

DcoConstraint::DcoConstraint(OsiLorentzCone const * cone): cType_(DCO_CONIC) {
  cone_ = new OsiLorentzCone(*cone);
}

DcoConstraint::DcoConstraint(DcoConstraint const & other): cType_(other.type()) {
}

DcoConstraint::~DcoConstraint() {
}

DcoConstraintType DcoConstraint::type() const {
  return cType_;
}

// assumes linear constraint
int DcoConstraint::getSize() const {
  if (cType_==DCO_CONIC) {
    throw std::exception();
  }
  return size_;
}

/** Create a OsiRowCut based on this constraint. */
OsiRowCut * DcoConstraint::createOsiRowCut() const {
  double lower = CoinMax(getLbHard(), getLbSoft());
  double upper = CoinMin(getUbHard(), getUbSoft());
  OsiRowCut * cut = new OsiRowCut;
  if (!cut) {
    /* Out of memory. */
    // message_handler = ;
    // messages = ;
    // message_handler->message(DISCO_OUT_OF_MEMORY, *messages)
    //   << __FILE__ << __LINE__ << CoinMessageEol;
    throw CoinError("Out of memory", "createOsiRowCut", "DcoConstraint");
  }
  assert(size_ > 0);
  cut->setLb(lower);
  cut->setUb(upper);
  cut->setRow(size_, indices_, values_);
  return cut;
}
