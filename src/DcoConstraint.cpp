#include "DcoConstraint.hpp"

#include "DcoModel.hpp"

DcoConstraint::DcoConstraint(double lb, double ub):
  BcpsConstraint(lb, ub, lb, ub) {
  type_ = DcoConstraintTypeNotSet;
}

DcoConstraint::~DcoConstraint() {
}
