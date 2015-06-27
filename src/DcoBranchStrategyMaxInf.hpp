/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * ALPS is distributed under the Eclipse Public License as part of the       *
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
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef DcoBranchStrategyMaxInf_h_
#define DcoBranchStrategyMaxInf_h_

#include "BcpsBranchObject.h"
#include "BcpsBranchStrategy.h"
#include "DcoModel.hpp"

/** This class implements maximum infeasibility branching. */
class DcoBranchStrategyMaxInf : public BcpsBranchStrategy {

 private:

  /** Illegal Assignment operator.*/
    DcoBranchStrategyMaxInf& operator=(const DcoBranchStrategyMaxInf& rhs);
    
 public:
    
    /** MaxInf Constructor. */
    DcoBranchStrategyMaxInf() {
	type_ = static_cast<int>(DcoBranchingStrategyMaxInfeasibility);
    }

    /** MaxInf Constructor. */
    DcoBranchStrategyMaxInf(DcoModel *model) : BcpsBranchStrategy(model) {
	type_ = static_cast<int>(DcoBranchingStrategyMaxInfeasibility);
    }
    
    /** Destructor. */
    virtual ~DcoBranchStrategyMaxInf() {}
    
    /** Copy constructor. */
    DcoBranchStrategyMaxInf(const DcoBranchStrategyMaxInf &);
    
    /** Clone a brancing strategy. */
    virtual BcpsBranchStrategy * clone() const {
	return new DcoBranchStrategyMaxInf(*this);
    }
    
    /** Create a set of candidate branching objects. */
    virtual int createCandBranchObjects(int numPassesLeft, double ub);
    
    /** Compare branching object thisOne to bestSoFar. If thisOne is better 
	than bestObject, return branching direction(1 or -1), otherwise
	return 0. 
	If bestSorFar is NULL, then always return branching direction(1 or -1).
    */
    virtual int betterBranchObject(BcpsBranchObject * thisOne,
				   BcpsBranchObject * bestSoFar);
};

#endif
