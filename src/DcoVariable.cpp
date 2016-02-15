#include "DcoVariable.hpp"

DcoVariable::DcoVariable() {
}

DcoVariable::DcoVariable(double lbh, double ubh, double lbs, double ubs) {
}

DcoVariable::DcoVariable(double lbh, double ubh, double lbs, double ubs,
			 DcoIntegralityType it): BcpsVariable(lbh, ubh, lbs, ubs) {
  BcpsIntegral_t type;
  if (it==DcoIntegralityTypeCont) {
    type = 'C';
  }
  else {
    type = 'I';
  }
  setIntType(type);
}

DcoVariable::DcoVariable(DcoVariable const & other) {
}

DcoVariable::~DcoVariable() {
}
