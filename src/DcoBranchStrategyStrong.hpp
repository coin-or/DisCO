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


#ifndef DcoBranchStrategyStrong_h_
#define DcoBranchStrategyStrong_h_

#include <BcpsBranchObject.h>
#include <BcpsBranchStrategy.h>

#include "DcoModel.hpp"
#include "DcoBranchObject.hpp"

/*!
  Implements strong branching.
*/

class DcoBranchStrategyStrong : virtual public BcpsBranchStrategy {
  /// update score for the given branch object.
  void updateScore(BcpsBranchObject * bobject, double orig_lb,
                   double orig_ub, double orig_obj) const;
  // return integer infeasibility for the given value
  double infeas(double value) const;
 public:
  ///@name Constructor and Destructor.
  //@{
  /// Constructor.
  DcoBranchStrategyStrong(DcoModel * model);
  /// Destructor.
  virtual ~DcoBranchStrategyStrong() {}
  //@}

  ///@name Selecting and Creating branches.
  //@{
  /// Create a set of candidate branching objects from the given node.
  virtual int createCandBranchObjects(BcpsTreeNode * node);
  /// Compare current to other, return 1 if current is better, 0 otherwise
  virtual int betterBranchObject(BcpsBranchObject const * current,
                                 BcpsBranchObject const * other);
  //@}

 private:
  /// Disable default constructor.
  DcoBranchStrategyStrong();
  /// Disable copy constructor.
  DcoBranchStrategyStrong(DcoBranchStrategyStrong const & other);
  /// Disable copy assignment operator.
  DcoBranchStrategyStrong & operator=(DcoBranchStrategyStrong const & rhs);
};

#endif
