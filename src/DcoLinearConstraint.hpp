#ifndef DcoLinearConstraint_hpp_
#define DcoLinearConstraint_hpp_

#include "DcoConstraint.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"

/*!
  DcoLinearConstraint represents linear constraint.
 */
class DcoLinearConstraint: virtual public DcoConstraint {
  //@{
  /// Number of variables with nonzero coefficient.
  int size_;
  /// Indices of non-zero coefficients.
  int * indices_;
  /// Values of non-zero coefficients.
  double * values_;
  //@}
public:
  DcoLinearConstraint(int size, int const * indices, double const * values,
		      double lb, double ub);
  DcoLinearConstraint(DcoLinearConstraint const & other);
  DcoLinearConstraint & operator=(DcoLinearConstraint const & rhs);
  virtual ~DcoLinearConstraint();
  int getSize() const;
  int const * getIndices() const;
  double const * getValues() const;
  virtual OsiRowCut * createOsiRowCut(DcoModel * model);
};

#endif
