#include "DcoConstraint.hpp"

#include "DcoModel.hpp"

DcoConstraint::DcoConstraint(double lb, double ub):
  BcpsConstraint(lb, ub, lb, ub) {
}

DcoConstraint::DcoConstraint(DcoConstraint const & other) {
}

DcoConstraint DcoConstraint::operator=(DcoConstraint const & rhs) {
}

DcoConstraint::~DcoConstraint() {
}
