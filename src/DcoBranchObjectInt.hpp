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
 *                                                                           *
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


//#############################################################################
// Borrow ideas from COIN/Cbc
//#############################################################################


#include "BcpsBranchObject.h"

#include "DcoModel.hpp"


//#############################################################################


class DcoBranchObjectInt : public BcpsBranchObject {

 protected:

    /** Down_[0]: the lower bound of down arm;
	Down_[1]: the upper bound of down arm; */
    double downBranchBounds_[2];

    /** Up_[0]: the lower bound of upper arm;
	Up_[1]: the upper bound of upper arm; */
    double upBranchBounds_[2];

 public:

    /** Default constructor. */
    DcoBranchObjectInt()
	:
	BcpsBranchObject()
	{
	    type_ = DcoBranchingObjectTypeInt;
	    downBranchBounds_[0] = 0.0;
	    downBranchBounds_[1] = 0.0;
	    upBranchBounds_[0] = 0.0;
	    upBranchBounds_[1] = 0.0;
	}

    /** Construct a branching object, which branching on variable varInd.
	\param varInd     the index of integer variable in object set
	\param direction  the direction of first branching: 1(up), -1(down)
	\param value      the fractional solution value of variable varInd
    */
    DcoBranchObjectInt(DcoModel * model,
			int varInd,
			int direction,
			double value)
	:
	BcpsBranchObject(model, varInd, direction, value)
	{
	    type_ = DcoBranchingObjectTypeInt;
	    int iColumn = model->getIntColIndices()[objectIndex_];
	    downBranchBounds_[0] = model->solver()->getColLower()[iColumn];
	    downBranchBounds_[1] = floor(value_);
	    upBranchBounds_[0] = ceil(value_);
	    upBranchBounds_[1] = model->getColUpper()[iColumn];
	}

    /** Construct a branching object, which branching on variable varInd.
	\param varInd     the index of integer variable in object set
	\param intScore   the integer score/goodness
	\param dblScore   the double score/goodness
	\param direction  the direction of first branching: 1(up), -1(down)
	\param value      the fractional solution value of variable varInd
    */
    DcoBranchObjectInt(DcoModel * model,
			int varInd,
			int intScore,
			double dblScore,
			int direction,
			double value)
	:
	BcpsBranchObject(model, varInd, intScore, dblScore, direction, value)
	{
	    type_ = DcoBranchingObjectTypeInt;
	    int iColumn = model->getIntColIndices()[objectIndex_];
	    downBranchBounds_[0] = model->solver()->getColLower()[iColumn];
	    downBranchBounds_[1] = floor(value_);
	    upBranchBounds_[0] = ceil(value_);
	    upBranchBounds_[1] = model->getColUpper()[iColumn];
	}

    /** Create a degenerate branching object.
	Specifies a `one-direction branch'. Calling branch() for this
	object will always result in lowerValue <= x <= upperValue.
	Used to fix a variable when lowerValue = upperValue.
    */
    DcoBranchObjectInt(DcoModel * model,
			int varInd,
			int direction,
			double lowerValue,
			double upperValue)
	:
	BcpsBranchObject(model, varInd, direction, lowerValue)
	{
	    type_ = DcoBranchingObjectTypeInt;
	    numBranchesLeft_ = 1;
	    downBranchBounds_[0] = lowerValue;
	    downBranchBounds_[1] = upperValue;
	    upBranchBounds_[0] = lowerValue;
	    upBranchBounds_[1] = upperValue;
	}

    /** Copy constructor. */
    DcoBranchObjectInt(const DcoBranchObjectInt &);

    /** Assignment operator. */
    DcoBranchObjectInt & operator = (const DcoBranchObjectInt& rhs);

    /** Clone. */
    virtual BcpsBranchObject * clone() const {
	return (new DcoBranchObjectInt(*this));
    }

    /** Destructor. */
    virtual ~DcoBranchObjectInt() {}

    /** Set the bounds for the variable according to the current arm
	of the branch and advances the object state to the next arm.
	Returns change in guessed objective on next branch. */
    virtual double branch(bool normalBranch = false);

    /** \brief Print something about branch - only if log level high. */
    virtual void print(bool normalBranch);

    /** Get down arm bounds. */
    const double *getDown() const { return downBranchBounds_; }

    /** Get upper arm bounds. */
    const double *getUp() const { return upBranchBounds_; }

 protected:

    /** Pack Disco portion to an encoded object. */
    AlpsReturnStatus encodeDco(AlpsEncoded *encoded) const {
	assert(encoded);
	AlpsReturnStatus status = AlpsReturnStatusOk;
	int j;
	// TODO: N-way.
	for (j = 0; j < 2; ++j) {
	    encoded->writeRep(downBranchBounds_[j]);
	}
	for (j = 0; j < 2; ++j) {
	    encoded->writeRep(upBranchBounds_[j]);
	}

	return status;
    }

    /** Unpack Disco portion from an encoded object. */
    AlpsReturnStatus decodeDco(AlpsEncoded &encoded) {
	AlpsReturnStatus status = AlpsReturnStatusOk;
	int j;
	// TODO: N-way.
	for (j = 0; j < 2; ++j) {
	    encoded.readRep(downBranchBounds_[j]);
	}
	for (j = 0; j < 2; ++j) {
	    encoded.readRep(upBranchBounds_[j]);
	}

	return status;
    }

 public:

    /** Pack to an encoded object. */
    virtual AlpsReturnStatus encode(AlpsEncoded *encoded) const {
	AlpsReturnStatus status = AlpsReturnStatusOk;

	status = encodeBcps(encoded);
	status = encodeDco(encoded);

	return status;
    }

    /** Unpack a branching object from an encoded object. */
    virtual AlpsReturnStatus decode(AlpsEncoded &encoded) {

	AlpsReturnStatus status = AlpsReturnStatusOk;

	status = decodeBcps(encoded);
	status = decodeDco(encoded);

	return status;
    }

};
