#ifndef DcoConstraint_hpp_
#define DcoConstraint_hpp_

#include <OsiLorentzCone.hpp>
#include <BcpsObject.h>
#include "Dco.hpp"

class OsiRowCut;
class DcoModel;
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

       It has the following fields, \c objectIndex_, \c repType_, \c intType_,
       \c validRegion_, \c status_, \c lbHard_, \c ubHard_, \c lbSoft_,
       \c ubSoft_, \c hashValue_, \c numInactive_, \c effectiveness_.



  <li> BcpsConstraint<br>
       Represents a constraint. Does not assume anything for a constraint
       other than having lower and upper bounds. This class does not have any
       fields other than the ones inherited.

  <li> DcoConstraint<br>
       Represents DisCO constraints. DisCO constraints comes in two types,
       (1) linear constraints (2) conic constraints. Conic constraints are
       of two types, Lorentz cones or rotated Lorentz Cones. We do not support
       generic conic constraints (|Ax-b| <= d^Tx-h) yet.

       This class is a base class for linear and conic constraints.
  </ul>

  todo(aykut) list:
  <ul>
  <li> Find a way to log messages. We need model (DcoModel) pointer for this.

  <li> createOsiRowCut() should be implemented in Bcps level.

  </ul>
 */

class DcoConstraint: public BcpsConstraint {
public:
  DcoConstraint() {}
  DcoConstraint(double lb, double ub);
  virtual ~DcoConstraint();
  /// Create an OsiRowCut based on this constraint. Returns NULL if this is a
  /// conic constraint.
  virtual OsiRowCut * createOsiRowCut(DcoModel * model) const = 0;
  /// return constraint type, linear or conic
  virtual DcoConstraintType constraintType() const = 0;
private:
};

#endif
