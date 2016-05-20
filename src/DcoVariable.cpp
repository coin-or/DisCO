#include "DcoVariable.hpp"
#include "DcoModel.hpp"
#include "DcoBranchObject.hpp"

#include <cmath>

DcoVariable::DcoVariable() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

DcoVariable::DcoVariable(int index, double lbh, double ubh, double lbs, double ubs):
  BcpsVariable (lbh, ubh, lbs, ubs) {
  setObjectIndex(index);
}

DcoVariable::DcoVariable(int index, double lbh, double ubh, double lbs, double ubs,
                         DcoIntegralityType it): BcpsVariable(lbh, ubh, lbs, ubs) {
  BcpsIntegral_t type;
  if (it==DcoIntegralityTypeCont) {
    type = 'C';
  }
  else {
    type = 'I';
  }
  setIntType(type);
  setObjectIndex(index);
}

DcoVariable::DcoVariable(DcoVariable const & other) {
}

DcoVariable::~DcoVariable() {
}

BcpsObject * DcoVariable::clone() const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

double DcoVariable::infeasibility(BcpsModel * bcps_model,
                                  int & preferredDir) const {
  DcoModel * model = dynamic_cast<DcoModel*>(bcps_model);
  preferredDir = -1;
  // get integer tolerance parameter
  double tolerance = model->dcoPar()->entry(DcoParams::integerTol);
  double value = model->solver()->getColSolution()[getObjectIndex()];
  double dist_to_upper = ceil(value) - value;
  double dist_to_lower = value - floor(value);
  // return the minimum of distance to upper or lower
  double infeas;
  if (dist_to_upper>dist_to_lower) {
    preferredDir = -1;
    infeas = dist_to_lower;
  }
  else {
    preferredDir = 1;
    infeas = dist_to_upper;
  }
  if (infeas<tolerance) {
    infeas = 0.0;
  }
  return infeas;
}

BcpsBranchObject * DcoVariable::createBranchObject(BcpsModel * bcps_model,
                                                   int direction) const {
  DcoModel * model = dynamic_cast<DcoModel*>(bcps_model);
  int var_index = getObjectIndex();
  // get current value from solver
  double value = model->solver()->getColSolution()[var_index];
  // todo(aykut) we do not know the score of branch object since we do not know
  // the branching strategy.
  double score = 0.0;
  BcpsBranchObject * bo = new DcoBranchObject(model, var_index, score, value);
  //BcpsBranchObject(BcpsModel * model, int objectIndex, int upScore,
  //                 double downScore, int direction, double value)
  return bo;
}
