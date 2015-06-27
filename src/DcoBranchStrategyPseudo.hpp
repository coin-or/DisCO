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


//#############################################################################
// NOTE: Borrow ideas from COIN/Cbc
//#############################################################################


#ifndef DcoBranchStrategyPseudo_h_
#define DcoBranchStrategyPseudo_h_

#include "BcpsBranchObject.h"
#include "BcpsBranchStrategy.h"
#include "DcoModel.hpp"


/** Dco branching strategy.
    This class implements pseudocost branching. */
class DcoBranchStrategyPseudo : public BcpsBranchStrategy {

 private:
    /** Illegal Assignment operator.*/
    DcoBranchStrategyPseudo& operator=(const DcoBranchStrategyPseudo& rhs);

    int relibility_;
    
 public:

    /** Default Constructor. */
    DcoBranchStrategyPseudo() : relibility_(1) {
	relibility_ = 1;
	type_ = static_cast<int>(DcoBranchingStrategyPseudoCost);
    }

    /** Useful Constructor. */
    DcoBranchStrategyPseudo(DcoModel *model, int rel)
	: BcpsBranchStrategy(model) {
        relibility_ = rel;
	type_ = static_cast<int>(DcoBranchingStrategyPseudoCost);
    }

    /** Destructor. */
    virtual ~DcoBranchStrategyPseudo() {}
    
    /** Copy constructor. */
    DcoBranchStrategyPseudo(const DcoBranchStrategyPseudo &);
    
    /** Set relibility. */
    void setRelibility(int rel) { relibility_ = rel; }    

    /** Clone a brancing strategy. */
    virtual BcpsBranchStrategy * clone() const {
	return new DcoBranchStrategyPseudo(*this);
    }
    
    /** Compare branching object thisOne to bestSoFar. If thisOne is better 
	than bestObject, return branching direction(1 or -1), otherwise
	return 0. 
	If bestSorFar is NULL, then always return branching direction(1 or -1).
    */
    virtual int betterBranchObject(BcpsBranchObject * thisOne,
				   BcpsBranchObject * bestSoFar);

    /** Create a set of candidate branching objects. */
    int createCandBranchObjects(int numPassesLeft, double ub);
};

#endif
