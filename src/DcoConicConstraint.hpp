#ifndef DcoConicConstraint_hpp_
#define DcoConicConstraint_hpp_

#include "DcoConstraint.hpp"
#include "Dco.hpp"
#include "DcoLinearConstraint.hpp"

/*!
  DcoConicConstraint represents conic constraint. We consider Lorentz cones
  and rotated Lorentz cones only for now.
 */
class DcoConicConstraint: virtual public DcoConstraint {
  /// Cone type.
  DcoLorentzConeType coneType_;
  /// Cone size.
  int size_;
  /// Cone members.
  int * members_;
  /// Number of supports.
  int numSupports_;
  /// Array of linear constraints that approximates this conic constraint.
  DcoLinearConstraint ** supports_;
  /// Active (binding) supports. 1 for active and 0 otherwise.
  int * activeSupports_;
public:
  /// Initializes constraint from given data.
  DcoConicConstraint(DcoLorentzConeType type, int size, int const * members);
  /// Copy constructor.
  DcoConicConstraint(DcoConicConstraint const & other);
  /// Copy assignment operator.
  DcoConicConstraint & operator=(DcoConicConstraint const & rhs);
  /// Destructor.
  virtual ~DcoConicConstraint();
  /// Returns type of conic constraint.
  DcoLorentzConeType getType() const;
  /// Return size of cone, i.e., number of variables in the cone.
  int getSize() const;
  /// Return array of cone members.
  int const * getMembers() const;
  /// Return number of supports.
  int getNumSupports() const;
  /// Return array of supports that approximates this conic constraint.
  DcoLinearConstraint const * const * getSupports() const;
  /// Return array that gives active (binding) supports. Support i is active if
  /// getSupports()[i] is 1, inactive if 0.
  int const * getActiveSupports() const;
};

#endif
