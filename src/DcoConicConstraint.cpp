#include "DcoConicConstraint.hpp"

/// Initializes constraint from given data.
DcoConicConstraint::DcoConicConstraint(DcoLorentzConeType type, int size,
				       int const * members):
  DcoConstraint(0.0, DISCO_INFINITY) {
  coneType_ = type;
  size_ = size;
  members_ = new int[size_];
  std::copy(members, members+size, members_);
  numSupports_ = 0;
  supports_ = NULL;
  activeSupports_ = NULL;
}

/// Copy constructor.
DcoConicConstraint::DcoConicConstraint(DcoConicConstraint const & other):
  DcoConstraint(other) {
  coneType_ = other.getType();
  size_ = other.getSize();
  members_ = new int[size_];
  int const * other_members = other.getMembers();
  std::copy(other_members, other_members+size_, members_);
  // copy support data
  numSupports_ = other.getNumSupports();
  DcoLinearConstraint const * const * other_supports = other.getSupports();
  supports_ = new DcoLinearConstraint*[numSupports_];
  for (int i=0; i<numSupports_; ++i) {
    supports_[i] = new DcoLinearConstraint(*other_supports[i]);
  }
  activeSupports_ = new int[size_];
  int const * other_as = other.getActiveSupports();
  std::copy(other_as, other_as+numSupports_, activeSupports_);
}

/// Copy assignment operator.
DcoConicConstraint &
DcoConicConstraint::operator=(DcoConicConstraint const & rhs) {
  type_ = rhs.getType();
  size_ = rhs.getSize();
  members_ = new int[size_];
  int const * rhs_members = rhs.getMembers();
  std::copy(rhs_members, rhs_members+size_, members_);
  // copy support data
  numSupports_ = rhs.getNumSupports();
  DcoLinearConstraint const * const * rhs_supports = rhs.getSupports();
  supports_ = new DcoLinearConstraint*[numSupports_];
  for (int i=0; i<numSupports_; ++i) {
    supports_[i] = new DcoLinearConstraint(*rhs_supports[i]);
  }
  activeSupports_ = new int[size_];
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

/// Returns type of conic constraint.
DcoLorentzConeType DcoConicConstraint::getType() const {
  return coneType_;
}

/// Return size of cone, i.e., number of variables in the cone.
int DcoConicConstraint::getSize() const {
  return size_;
}

/// Return array of cone members.
int const * DcoConicConstraint::getMembers() const {
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
