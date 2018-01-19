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


#ifndef DcoVariable_hpp_
#define DcoVariable_hpp_

#include <BcpsObject.h>
#include "Dco.hpp"

/*!
   Represents a DisCO variable. DcoVariable inherits BcpsVariable. BcpsVariable
   inherits BcpsObject. BcpsObject inherits AlpsKnowledge.
   DcoVariable -> BcpsVariable -> BcpsObject -> AlpsKnowledge.

   # AlpsKnowledge
   AlpsKnowledge is an abstract base class of any user-defined class that Alps
   has to know about in order to encode/decode. It has two fields
   (AlpsEncoded *) encoded_ and (KnowledgeType) type_.

   # BcpsObject
   BcpsObject is a class for describing the objects that comprise a
   mathematical optimization problem. All that is assumed about an object is
   that it has bounds and might have integrality constraint.

   Bcps object has the following fields

   <ul>
   <li> objectIndex_
   <li> repType_
   <li> intType_
   <li> validRegion_
   <li> status_
   <li> lbHard_
   <li> ubHard_
   <li> lbSoft_
   <li> ubSoft_
   <li> hashValue_
   <li> numInactive_
   <li> effectiveness_
   </ul>


   BcpsVariable:
   BcpsVariable is a class that represents a variable. It does not have any
   field other than ones inherited from BcpsObject.

   In DcoModel columns are stored in cols_ inherited from BcpsModel. We keep
   number of integer variables at numIntegerCols_ and their indices at (int
   * intColIndices_). intColIndices_[0] gives the index of the first integer
   column in the cols_ array.

   DcoModel objects are kept at constraints_ and variables_ inherited from
   BcpsModel.

   In Blis, integer variables have their own class, BlisObjectInt.
   BlisObjectInt inherits BcpsObject class.
 */

class DcoVariable: virtual public BcpsVariable {
  // double objCoef_;
  // int size_;
  // int * indices_;
  // double * values_;
public:
  ///@name Cosntructors and Destructors
  //@{
  /// Default constructor.
  DcoVariable();
  /// Constructor with lower and upper bounds.
  DcoVariable(int index, double lbh, double ubh, double lbs, double ubs);
  // /// Constructor with bounds and integrality type.
  // DcoVariable(int index, double lbh, double ubh, double lbs, double ubs,
  //             DcoIntegralityType it, double objCoef, int size,
  //             int const * indices, double const * values);
  /// Copy constructor.
  DcoVariable(DcoVariable const & other);
  /// Destructor.
  virtual ~DcoVariable();
  //@}

  ///@ Get column information
  //@{
  /// Get objective coefficient
  // double objCoef() const { return objCoef_; }
  // /// Get size.
  // int size() const { return size_; }
  // /// Get indices.
  // int const * indices() const { return indices_; }
  // /// Get values.
  // double const * values() const { return values_; }
  //@}

  ///@name Virtual Functions inherited from BcpsObject
  //@{
  virtual BcpsObject * clone() const;
  /// Return the infeasibility of the variable. A variable can be infeasible
  /// only when it is integer.
  // todo(aykut) Why do we need model as input?
  virtual double infeasibility(BcpsModel * bcps_model, int & preferredDir) const;
  // virtual void feasibleRegion(BcpsModel *m) {}
  /// Create a branch object from this variable.
  virtual BcpsBranchObject * createBranchObject(BcpsModel * bcps_model, int way) const;

  // virtual BcpsBranchObject * preferredNewFeasible(BcpsModel *m) const;
  // virtual BcpsBranchObject * notPreferredNewFeasible(BcpsModel *m) const;
  // virtual void resetBounds(BcpsModel *m);
  // virtual bool boundBranch(BcpsModel *m) const;
  // virtual void floorCeiling(double &floorValue,
  //                        double &ceilingValue,
  //                        double value,
  //                        double tolerance) const;
  // virtual double upEstimate() const;
  // virtual double downEstimate() const;
  // virtual void printDesc();
  //@}

  ///@name Encode and Decode functions
  //@{
  /// Get encode from #AlpsKnowledge
  using AlpsKnowledge::encode;
  /// Encode this to an AlpsEncoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// Decode a given AlpsEncoded object to a new DcoVariable object and return
  /// a pointer to it.
  virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
  /// Decode a given AlpsEncoded object into self.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
  //@}
};

#endif
