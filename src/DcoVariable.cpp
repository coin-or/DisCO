/*===========================================================================*
 * This file is part of the Discrete Conic Optimization (DisCO) Solver.      *
 *                                                                           *
 * DisCO is distributed under the Eclipse Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *          Aykut Bulut, Lehigh University                                   *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and      *
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#include "DcoVariable.hpp"
#include "DcoModel.hpp"
#include "DcoBranchObject.hpp"

#include <cmath>

DcoVariable::DcoVariable() {
  // objCoef_ = 0.0;
  // size_ = 0;
  // indices_ = NULL;
  // values_ = NULL;
}

DcoVariable::DcoVariable(int index, double lbh, double ubh, double lbs, double ubs):
  BcpsVariable (lbh, ubh, lbs, ubs) {
  setObjectIndex(index);
}

// DcoVariable::DcoVariable(int index, double lbh, double ubh, double lbs,
//                          double ubs, DcoIntegralityType it, double objCoef,
//                          int size, int const * indices, double const * values)
//   : BcpsVariable(lbh, ubh, lbs, ubs) {
//   BcpsIntegral_t type;
//   if (it==DcoIntegralityTypeCont) {
//     type = 'C';
//   }
//   else {
//     type = 'I';
//   }
//   setIntType(type);
//   setObjectIndex(index);
//   objCoef_ = objCoef;
//   size_ = size;
//   indices_ = new int[size];
//   std::copy(indices, indices+size, indices_);
//   values_ = new double[size];
//   std::copy(values, values+size, values_);
// }

DcoVariable::DcoVariable(DcoVariable const & other) {
}

DcoVariable::~DcoVariable() {
  // if (indices_) {
  //   delete[] indices_;
  // }
  // if (values_) {
  //   delete[] values_;
  // }
}

BcpsObject * DcoVariable::clone() const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

//todo(aykut) what if the variable is not integer/relaxed?
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
  BcpsBranchObject * bo = new DcoBranchObject(var_index, score, value);
  bo->setBroker(broker_);
  //BcpsBranchObject(BcpsModel * model, int objectIndex, int upScore,
  //                 double downScore, int direction, double value)
  return bo;
}

/// Encode this to an AlpsEncoded object.
AlpsReturnStatus DcoVariable::encode(AlpsEncoded * encoded) const {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  return AlpsReturnStatusOk;
}

/// Decode a given AlpsEncoded object to a new DcoVariable object and return
/// a pointer to it.
AlpsKnowledge * DcoVariable::decode(AlpsEncoded & encoded) const {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  DcoVariable * new_var = new DcoVariable();
  return new_var;
}

/// Decode a given AlpsEncoded object into self.
AlpsReturnStatus DcoVariable::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  return status;
}
