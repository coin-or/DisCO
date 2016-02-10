#ifndef DcoConstraint_hpp_
#define DcoConstraint_hpp_

#include <OsiLorentzCone.hpp>
#include <BcpsObject.h>

typedef enum {
  DCO_LINEAR=0,
  DCO_CONIC
} DcoConstraintType;

class OsiRowCut;

/*!
  DcoConstraint inherits BcpsConstraint. BcpsConstraint inherits BcpsObject.
  BcpsObject inherits AlpsKnowledge.
  DcoConstraint -> BcpsConstraint -> BcpsObject -> AlpsKnowledge.

  <ul>
  <li> AlpsKnowledge<br>
       AlpsKnowledge is an abstract base class of any user-defined class that
       Alps has to know about in order to encode/decode. It has two fields
       (AlpsEncoded *) and (KnowledgeType) type_.

  <li> BcpsObject<br>
       BcpsObjects class represents a generic type for mathematical
       optimization problem object. All that is assumed about an object is
       that it has bounds and might have integrality constraint. It can be a
       variable or a constraint row (a mathematical formula of variables).

  <li> BcpsConstraint<br>
       Represents a constraint. Does not assume anything for a constraint
       other than having lower and upper bounds. This class does not have any
       fields other than the ones inherited.

  <li> DcoConstraint<br>
       Represents DisCO constraints. DisCO constraints comes in two types,
       (1) linear constraints (2) conic constraints. Conic constraints are
       of two types, Lorentz cones or rotated Lorentz Cones. We do not support
       generic conic constraints (|Ax-b| <= d^Tx-h) yet.

       If the constraint represented is a linear constraint we use fields
       \c size_, \c indices_ and \c values_. If it is a conic constraint
       field \c cone_ is used. Field \c cType_ keeps the type of the
       constraint. The irrelevant fields are ignored depending on the
       constraint type.
  </ul>

  todo(aykut) list:
  <ul>
  <li> Find a way to log messages. We need model (DcoModel) pointer for this.
  </ul>
 */

class DcoConstraint: public BcpsConstraint {
  DcoConstraintType const cType_;
  int size_;
  int * indices_;
  double * values_;
  OsiLorentzCone * cone_;
public:
  DcoConstraint();
  DcoConstraint(double lbh, double ubh, double lbs, double ubs,
		DcoConstraintType type);
  DcoConstraint(OsiLorentzCone const * cone);
  DcoConstraint(DcoConstraint const & other);
  virtual ~DcoConstraint();
  DcoConstraintType type() const;
  int getSize() const;
  /** Create a OsiRowCut based on this constraint. */
  OsiRowCut * createOsiRowCut() const;
};

#endif
