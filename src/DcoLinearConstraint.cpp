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


#include "DcoLinearConstraint.hpp"

#include <CoinHelperFunctions.hpp>
#include <OsiRowCut.hpp>

DcoLinearConstraint::DcoLinearConstraint() {
  size_ = 0;
  indices_ = NULL;
  values_ = NULL;
}

DcoLinearConstraint::DcoLinearConstraint(int size, int const * indices,
                                         double const * values, double lb,
                                         double ub):
  DcoConstraint(lb, ub) {
  size_ = size;
  indices_ = new int[size];
  std::copy(indices, indices+size, indices_);
  values_ = new double[size];
  std::copy(values, values+size, values_);
}

DcoLinearConstraint::DcoLinearConstraint(DcoLinearConstraint const & other):
  DcoConstraint(other) {
  size_ = other.getSize();
  // set indices
  indices_ = new int[size_];
  int const * other_indices = other.getIndices();
  std::copy(other_indices, other_indices+size_, indices_);
  // set values
  values_ = new double[size_];
  double const * other_values = other.getValues();
  std::copy(other_values, other_values+size_, values_);
}

DcoLinearConstraint &
DcoLinearConstraint::operator=(DcoLinearConstraint const & rhs) {
  size_ = rhs.getSize();
  // set indices
  indices_ = new int[size_];
  int const * rhs_indices = rhs.getIndices();
  std::copy(rhs_indices, rhs_indices+size_, indices_);
  // set values
  values_ = new double[size_];
  double const * rhs_values = rhs.getValues();
  std::copy(rhs_values, rhs_values+size_, values_);
  return *this;
}

DcoLinearConstraint::~DcoLinearConstraint() {
  if (indices_) {
    delete[] indices_;
  }
  if (values_) {
    delete[] values_;
  }
}

int DcoLinearConstraint::getSize() const {
  return size_;
}

int const * DcoLinearConstraint::getIndices() const {
  return indices_;
}

double const * DcoLinearConstraint::getValues() const {
  return values_;
}

/// Create a OsiRowCut based on this constraint.
OsiRowCut * DcoLinearConstraint::createOsiRowCut(DcoModel * model) const {
  double lower = CoinMax(getLbHard(), getLbSoft());
  double upper = CoinMin(getUbHard(), getUbSoft());
  OsiRowCut * cut = new OsiRowCut;
  if (!cut) {
    // Out of memory.
    model->dcoMessageHandler_->message(DISCO_OUT_OF_MEMORY,
                                            *(model->dcoMessages_))
      << __FILE__ << __LINE__ << CoinMessageEol;
    throw CoinError("Out of memory", "createOsiRowCut", "DcoConstraint");
  }
  assert(size_ > 0);
  cut->setLb(lower);
  cut->setUb(upper);
  cut->setRow(size_, indices_, values_);
  return cut;
}

double DcoLinearConstraint::infeasibility(BcpsModel * m,
                                          int & preferredWay) const {
  std::cerr << "Not implemented!" << std::endl;
  throw std::exception();
  return 0.0;
}


/// Encode this to an AlpsEncoded object.
AlpsReturnStatus DcoLinearConstraint::encode(AlpsEncoded * encoded) {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  return AlpsReturnStatusOk;
}

/// Decode a given AlpsEncoded object to an AlpsKnowledge object and return a
/// pointer to it.
AlpsKnowledge * DcoLinearConstraint::decode(AlpsEncoded & encoded) const {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  DcoLinearConstraint * con = new DcoLinearConstraint();
  // if (status!=AlpsReturnStatusOk) {
  //   std::cerr << "Unexpected decode status, "
  //             << "file: " <<  __FILE__
  //             << "line: " << __LINE__
  //             << std::endl;
  //   throw std::exception();
  // }
  return con;
}

// todo(aykut) this should be a pure virtual function in Alps level
// we can overload this function here due to cv-qualifier.
/// Decode a given AlpsEncoded object into self.
AlpsReturnStatus DcoLinearConstraint::decodeToSelf(AlpsEncoded & encoded) {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  return AlpsReturnStatusOk;
}
