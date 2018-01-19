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


#ifndef DcoBranchObject_hpp_
#define DcoBranchObject_hpp_

#include <BcpsBranchObject.h>
#include "DcoModel.hpp"

/*!
  Represents a DisCO branch object. DcoBranchObject inherits BcpsBranchObject.

  # DcoBranchObject

  DcoBranchObject represents a branch object for a simple branch on an
  integral variable.

  Has two fields, ubDownBranch_ and lbUpBranch_.

 */

class DcoBranchObject: virtual public BcpsBranchObject {
  /// upper bound of the down branch
  double ubDownBranch_;
  /// lower bound of the up branch
  double lbUpBranch_;
public:
  ///@name Constructor and Destructors.
  //@{
  /// Constructor.
  DcoBranchObject(int index, double score, double value);
  /// Copy constructor.
  DcoBranchObject(DcoBranchObject const & other);
  /// Helpful Copy constructor.
  DcoBranchObject(BcpsBranchObject const * other);
  /// Copy assignment operator
  DcoBranchObject & operator=(DcoBranchObject const & rhs);
  /// Destructor.
  virtual ~DcoBranchObject();
  //@}

  ///@name Other functions
  //@{
  /// Clone this to a new object and return pointer to it.
  //virtual BcpsBranchObject * clone() const;
  //@}

  ///@name Virtual functions inherited from BcpsBranchObject
  /// The number of branch arms created for this branch object.
  virtual int numBranches() const;
  /// The number of branch arms left to be evaluated.
  virtual int numBranchesLeft() const;
  /// Spit out a branch and, update this or superclass fields if necessary.
  virtual double branch(bool normalBranch = false);
  //@}

  ///@name Bound getting functions.
  //@{
  /// Get upper bound of the down branch.
  double ubDownBranch() const { return ubDownBranch_; }
  /// Get lower bound of the up branch.
  double lbUpBranch() const { return lbUpBranch_; }
  //@}

  ///@name Encode and Decode functions
  ///@{
  using AlpsKnowledge::encode;
  /// Encode the content of this into the given AlpsEncoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// Decode the given AlpsEncoded object into a new AlpsKnowledge object and
  /// return a pointer to it.
  virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
  /// Decode the given AlpsEncoded object into this.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
  ///@}

private:
  /// Disable default constructor.
  DcoBranchObject();
};

#endif
