/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
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
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#include "CoinHelperFunctions.hpp"

#include "DcoBranchObjectInt.hpp"


//#############################################################################

// Copy constructor
DcoBranchObjectInt::
DcoBranchObjectInt(const DcoBranchObjectInt & rhs)
    :
    BcpsBranchObject(rhs)
{
    downBranchBounds_[0] = rhs.downBranchBounds_[0];
    downBranchBounds_[1] = rhs.downBranchBounds_[1];
    upBranchBounds_[0] = rhs.upBranchBounds_[0];
    upBranchBounds_[1] = rhs.upBranchBounds_[1];
}

//#############################################################################

// Assignment operator
DcoBranchObjectInt &
DcoBranchObjectInt::operator=(const DcoBranchObjectInt& rhs)
{
    if (this != &rhs) {
	BcpsBranchObject::operator=(rhs);
	downBranchBounds_[0] = rhs.downBranchBounds_[0];
	downBranchBounds_[1] = rhs.downBranchBounds_[1];
	upBranchBounds_[0] = rhs.upBranchBounds_[0];
	upBranchBounds_[1] = rhs.upBranchBounds_[1];
    }
    return *this;
}

//#############################################################################

/*
  Perform a branch by adjusting the bounds of the specified variable.
  Note that each arm of the branch advances the object to the next arm by
  advancing the value of direction_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-direction branching.
  Returns change in guessed objective on next branch
*/
double
DcoBranchObjectInt::branch(bool normalBranch)
{
    DcoModel *model = dynamic_cast<DcoModel *>(model_);

    int iColumn = model->getIntColIndices()[objectIndex_];

    // Decrement number of branches left by 1.
    --numBranchesLeft_;

    if (direction_<0) {
#ifdef DISCO_DEBUG_MORE
	{
	    double olb, oub ;
	    olb = model->solver()->getColLower()[iColumn];
	    oub = model->solver()->getColUpper()[iColumn];
	    printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
		   iColumn,olb,oub,downBranchBounds_[0],downBranchBounds_[1]);
	}
#endif
	model->solver()->setColLower(iColumn, downBranchBounds_[0]);
	model->solver()->setColUpper(iColumn, downBranchBounds_[1]);
	direction_ = 1;
    }
    else {
#ifdef DISCO_DEBUG_MORE
	{
	    double olb, oub ;
	    olb = model->solver()->getColLower()[iColumn];
	    oub = model->solver()->getColUpper()[iColumn];
	    printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
		   iColumn,olb,oub,upBranchBounds_[0],upBranchBounds_[1]);
	}
#endif
	model->solver()->setColLower(iColumn, upBranchBounds_[0]);
	model->solver()->setColUpper(iColumn, upBranchBounds_[1]);
	direction_ = -1;	  // Swap direction
    }

    return 0.0;
}

//#############################################################################

// Print what would happen
void
DcoBranchObjectInt::print(bool normalBranch)
{
    DcoModel *model = dynamic_cast<DcoModel*>(model_);
    int iColumn = model->getIntColIndices()[objectIndex_];
    double olb, oub ;

    if (direction_ < 0) {
	olb = model->solver()->getColLower()[iColumn] ;
	oub = model->solver()->getColUpper()[iColumn] ;
	printf("DcoInteger would branch down on var %d: [%g,%g] => [%g,%g]\n",
	       iColumn,olb,oub,downBranchBounds_[0],downBranchBounds_[1]);

    }
    else {
	olb = model->solver()->getColLower()[iColumn] ;
	oub = model->solver()->getColUpper()[iColumn] ;
	printf("DcoInteger would branch up on var %d: [%g,%g] => [%g,%g]\n",
	       iColumn,olb,oub,upBranchBounds_[0],upBranchBounds_[1]);
    }
}

//#############################################################################
