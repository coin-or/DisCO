/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2009, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef BlisBranchStrategyBilevel_h_
#define BlisBranchStrategyBilevel_h_

#include "BcpsBranchObject.h"
#include "BcpsBranchStrategy.h"
#include "BlisModel.h"

/** This class implements maximum infeasibility branching. */
class BlisBranchStrategyBilevel : public BcpsBranchStrategy {

 private:

    /** Illegal Assignment operator.*/
    BlisBranchStrategyBilevel& operator=(const BlisBranchStrategyBilevel& rhs);
    
 public:
    
    /** Bilevel Constructor. */
    BlisBranchStrategyBilevel() {}
  
    /** Bilevel Constructor. */
    BlisBranchStrategyBilevel(BlisModel *model) : BcpsBranchStrategy(model) {}
  
    /** Destructor. */
    virtual ~BlisBranchStrategyBilevel() {}
  
    /** Copy constructor. */
    BlisBranchStrategyBilevel(const BlisBranchStrategyBilevel &);
    
    /** Clone a brancing strategy. */
    virtual BcpsBranchStrategy * clone() const {
	return new BlisBranchStrategyBilevel(*this);
    }
  
    /** Create a set of candidate branching objects. */
    virtual int createCandBranchObjects(int numPassesLeft);
  
    /** Compare branching object thisOne to bestSoFar. If thisOne is better 
	than bestObject, return branching direction(1 or -1), otherwise
	return 0. 
	If bestSorFar is NULL, then always return branching direction(1 or -1).
    */
    virtual int betterBranchObject(BcpsBranchObject * thisOne,
				   BcpsBranchObject * bestSoFar);
};

#endif
