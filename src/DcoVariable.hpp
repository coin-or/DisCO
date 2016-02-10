#ifndef DcoVariable_hpp_
#define DcoVariable_hpp_

#include <BcpsObject.h>
#include "Dco.hpp"

/**
   DcoVariable inherits BcpsVariable. BcpsVariable inherits BcpsObject.
   BcpsObject inherits AlpsKnowledge.
   DcoVariable -> BcpsVariable -> BcpsObject -> AlpsKnowledge.

   AlpsKnowledge is an abstract base class of any user-defined class that Alps
   has to know about in order to encode/decode. It has two fields
   (AlpsEncoded *) and (KnowledgeType) type_.

   BcpsObject is a class for describing the objects that comprise a
   mathematical optimization problem. All that is assumed about an object is
   that it has bounds and might have integrality constraint.

   BcpsVariable is a class that represents a variable. It does not have any
   field other than ones inherited from BcpsObject.

   DcoModel:
   Columns are stored in cols_ inherited from BcpsModel. We keep
   number of integer variables at numIntegerCols_ and their indices at (int
   * intColIndices_). intColIndices_[0] gives the index of the first integer
   in the cols_ array.

   We will not keep objects_ separately. We will use constraints_ and
   variables_ inherited from BcpsModel.


  // indices of integer columns
  // columns are stored in cols_ inherited from BcpsModel
  ;  // size of numIntObjects_


   Blis:
   Integer variables are represented with DcoObjectInt. DcoObjectInt inherits
   BcpsObject class.

 */
class DcoVariable: public BcpsVariable {
public:
  DcoVariable();
  DcoVariable(double lbh, double ubh, double lbs, double ubs);
  DcoVariable(double lbh, double ubh, double lbs, double ubs,
	      DcoIntegralityType it);
  DcoVariable(DcoVariable const & other);
  virtual ~DcoVariable();
};

#endif
