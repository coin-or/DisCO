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


#include "DcoConicConstraint.hpp"
#include <numeric>

DcoConicConstraint::DcoConicConstraint() {
  coneType_ = DcoLorentzCone;
  coneSize_ = 0;
  members_ = NULL;
  numSupports_ = 0;
  supports_ = NULL;
  activeSupports_ = NULL;
}

/// Initializes constraint from given data.
DcoConicConstraint::DcoConicConstraint(DcoLorentzConeType type, int size,
                                       int const * members):
  DcoConstraint(0.0, COIN_DBL_MAX) {
  coneType_ = type;
  coneSize_ = size;
  members_ = new int[coneSize_];
  std::copy(members, members+size, members_);
  numSupports_ = 0;
  supports_ = NULL;
  activeSupports_ = NULL;
}

/// Copy constructor.
DcoConicConstraint::DcoConicConstraint(DcoConicConstraint const & other):
  DcoConstraint(other) {
  coneType_ = other.coneType();
  coneSize_ = other.coneSize();
  members_ = new int[coneSize_];
  int const * other_members = other.coneMembers();
  std::copy(other_members, other_members+coneSize_, members_);
  // copy support data
  numSupports_ = other.getNumSupports();
  DcoLinearConstraint const * const * other_supports = other.getSupports();
  supports_ = new DcoLinearConstraint*[numSupports_];
  for (int i=0; i<numSupports_; ++i) {
    supports_[i] = new DcoLinearConstraint(*other_supports[i]);
  }
  activeSupports_ = new int[coneSize_];
  int const * other_as = other.getActiveSupports();
  std::copy(other_as, other_as+numSupports_, activeSupports_);
}

/// Copy assignment operator.
DcoConicConstraint &
DcoConicConstraint::operator=(DcoConicConstraint const & rhs) {
  coneType_ = rhs.coneType();
  coneSize_ = rhs.coneSize();
  members_ = new int[coneSize_];
  int const * rhs_members = rhs.coneMembers();
  std::copy(rhs_members, rhs_members+coneSize_, members_);
  // copy support data
  numSupports_ = rhs.getNumSupports();
  DcoLinearConstraint const * const * rhs_supports = rhs.getSupports();
  supports_ = new DcoLinearConstraint*[numSupports_];
  for (int i=0; i<numSupports_; ++i) {
    supports_[i] = new DcoLinearConstraint(*rhs_supports[i]);
  }
  activeSupports_ = new int[coneSize_];
  int const * rhs_as = rhs.getActiveSupports();
  std::copy(rhs_as, rhs_as+numSupports_, activeSupports_);
  return *this;
}

/// Destructor.
DcoConicConstraint::~DcoConicConstraint() {
  if (members_) {
    delete[] members_;
  }
  if (supports_) {
    for (int i=0; i<numSupports_; ++i) {
      delete supports_[i];
    }
    delete[] supports_;
  }
  if (activeSupports_) {
    delete[] activeSupports_;
  }
}

double DcoConicConstraint::infeasibility(BcpsModel * m,
                                         int & preferredWay) const {
  DcoModel * model = dynamic_cast<DcoModel*>(m);
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  double infeasibility;
  // get solution stored in solver
  double const * sol = model->solver()->getColSolution();
  // get related portion of the solution
  double * par_sol = new double[coneSize_];
  for(int i=0; i<coneSize_; ++i) {
    par_sol[i] = sol[members_[i]];
  }
  // get cone tolerance
  double cone_tol = model->dcoPar()->entry(DcoParams::coneTol);
  if (coneType_==DcoLorentzCone) {
    // infeasibility is
    // |x_2n| - x_1, if |x_2n| - x_1 > coneTol
    //  0 otherwise
    double * p = par_sol;
    double norm = std::inner_product(p+1, p+coneSize_, p+1, 0.0);
    norm = sqrt(norm);
    infeasibility = norm - p[0];
  }
  else if (coneType_==DcoRotatedLorentzCone) {
    // infeasibility is
    // |x_3n|^2 - 2x_1x_2, if |x_3n|^2 - 2x_1x_2, > coneTol
    //  0 otherwise
    double * p = par_sol;
    double ss = std::inner_product(p+2, p+coneSize_, p+2, 0.0);
    infeasibility = ss - 2.0*p[0]*p[1];
  }
  else {
    // unknown cone type.
    message_handler->message(DISCO_UNKNOWN_CONETYPE, *messages)
      << __FILE__ << __LINE__ << CoinMessageEol;
    throw std::exception();
  }
  if (infeasibility<=cone_tol) {
    infeasibility = 0.0;
  }
  delete[] par_sol;
  return infeasibility;
}

/// Returns type of conic constraint.
DcoLorentzConeType DcoConicConstraint::coneType() const {
  return coneType_;
}

/// Return size of cone, i.e., number of variables in the cone.
int DcoConicConstraint::coneSize() const {
  return coneSize_;
}

/// Return array of cone members.
int const * DcoConicConstraint::coneMembers() const {
  return members_;
}

/// Return number of supports.
int DcoConicConstraint::getNumSupports() const {
  return numSupports_;
}

/// Return array of supports that approximates this conic constraint.
DcoLinearConstraint const * const * DcoConicConstraint::getSupports() const {
  return supports_;
}

/// Return array that gives active (binding) supports. Support i is active if
/// getSupports()[i] is 1, inactive if 0.
int const * DcoConicConstraint::getActiveSupports() const {
  return activeSupports_;
}

/// Encode this to an AlpsEncoded object.
AlpsReturnStatus DcoConicConstraint::encode(AlpsEncoded * encoded) {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  return AlpsReturnStatusOk;
}

/// Decode a given AlpsEncoded object to an AlpsKnowledge object and return a
/// pointer to it.
AlpsKnowledge * DcoConicConstraint::decode(AlpsEncoded & encoded) const {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  DcoConicConstraint * con = new DcoConicConstraint();
  return con;
}

// todo(aykut) this should be a pure virtual function in Alps level
// we can overload this function here due to cv-qualifier.
/// Decode a given AlpsEncoded object into self.
AlpsReturnStatus DcoConicConstraint::decodeToSelf(AlpsEncoded & encoded) {
  std::cerr << "Not implemented, "
            << "file: " <<  __FILE__
            << "line: " << __LINE__
            << std::endl;
  throw std::exception();
  return AlpsReturnStatusOk;
}
