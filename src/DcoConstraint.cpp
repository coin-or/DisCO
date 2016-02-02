#include "DcoConstraint.hpp"

DcoConstraint::DcoConstraint(): cType_(DCO_LINEAR) {
}

DcoConstraint::DcoConstraint(double lbh, double ubh, double lbs,
			     double ubs, DcoConstraintType type):
  BcpsConstraint(lbh, ubh, lbs, ubs), cType_(type) {
}

DcoConstraint::DcoConstraint(DcoConstraint const & other): cType_(other.type()) {
}

DcoConstraint::~DcoConstraint() {
}

DcoConstraintType DcoConstraint::type() const {
  return cType_;
}
