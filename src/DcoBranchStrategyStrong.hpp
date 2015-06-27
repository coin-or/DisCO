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


#ifndef DcoBranchStrategyStrong_h_
#define DcoBranchStrategyStrong_h_

#include "BcpsBranchObject.h"
#include "BcpsBranchStrategy.h"
#include "DcoModel.hpp"


//#############################################################################


typedef struct {
    int objectIndex;            // object index
    BcpsBranchObject * bObject; // the branching object
    int numIntInfUp;            // without odd ones
    int numObjInfUp;            // just odd ones
    bool finishedUp;            // true if solver finished
    int numIntInfDown;          // without odd ones
    int numObjInfDown;          // just odd ones
    bool finishedDown;          // true if solver finished
} DcoStrong;


//#############################################################################


/** This class implements strong branching. */
class DcoBranchStrategyStrong : public BcpsBranchStrategy {

 private:

    /** Illegal Assignment operator.*/
    DcoBranchStrategyStrong& operator=(const DcoBranchStrategyStrong& rhs);
    
 public:
    
    /** Strong Constructor. */
    DcoBranchStrategyStrong() {
	type_ = static_cast<int>(DcoBranchingStrategyStrong);
    }

    /** Strong Constructor. */
    DcoBranchStrategyStrong(DcoModel *model)
        : BcpsBranchStrategy(model) {
	type_ = static_cast<int>(DcoBranchingStrategyStrong);
    }
    
    /** Destructor. */
    virtual ~DcoBranchStrategyStrong() {}
    
    /** Copy constructor. */
    DcoBranchStrategyStrong(const DcoBranchStrategyStrong &);
    
    /** Clone a brancing strategy. */
    virtual BcpsBranchStrategy * clone() const {
	return new DcoBranchStrategyStrong(*this);
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
