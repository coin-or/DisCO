#ifndef DcoConstraint_hpp_
#define DcoConstraint_hpp_

#include <BcpsObject.h>

typedef enum {
  DCO_LINEAR=0,
  DCO_CONIC
} DcoConstraintType;

class DcoConstraint: public BcpsConstraint {
  DcoConstraintType const cType_;
public:
  DcoConstraint();
  DcoConstraint(double lbh, double ubh, double lbs, double ubs,
		DcoConstraintType type);
  DcoConstraint(DcoConstraint const & other);
  virtual ~DcoConstraint();
  DcoConstraintType type() const;
};

#endif
