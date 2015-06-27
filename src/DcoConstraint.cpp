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

#include "OsiRowCut.hpp"

#include "DcoConstraint.hpp"
#include "DcoModel.hpp"

//#############################################################################

/* Default constructor. */
DcoConstraint::DcoConstraint() 
    :size_(0), indices_(NULL), values_(NULL) {}

//#############################################################################

/* Useful constructor. */
DcoConstraint::DcoConstraint(int s, const int *ind, const double *val)
{
    size_ = s;
    indices_ = new int [s];
    values_ = new double [s];
    memcpy(indices_, ind, s * sizeof(int));
    memcpy(values_, val, s * sizeof(double));
}

//#############################################################################

/* Useful constructor. */
DcoConstraint::DcoConstraint(double lbh, double ubh, double lbs, double ubs) 
    :
    BcpsConstraint(lbh, ubh, lbs, ubs),
    size_(0), indices_(NULL), values_(NULL) {}

//#############################################################################

/* Useful constructor. */
DcoConstraint::DcoConstraint(double lbh, double ubh, double lbs, double ubs,
			       int s, const int *ind, const double *val)
    : 
    BcpsConstraint(lbh, ubh, lbs, ubs)
{
    size_ = s;
    indices_ = new int [s];
    values_ = new double [s];
    memcpy(indices_, ind, s * sizeof(int));
    memcpy(values_, val, s * sizeof(double));
}

//#############################################################################

/** Destructor. */
DcoConstraint::~DcoConstraint()
{ 
    delete [] indices_; indices_ = NULL;
    delete [] values_; values_ = NULL;
}

//#############################################################################

/** Copy constructor. */
DcoConstraint::DcoConstraint(const DcoConstraint & rhs) 
    : BcpsConstraint(rhs) 
{
    size_ = rhs.size_;
    
    if (size_ < 0) {
	std::cout << "ERROR: size_ = " << size_ << std::endl;
	assert(size_);
    }
    if (size_ > 0) {
	indices_ = new int [size_];
	values_ = new double [size_];
	memcpy(indices_, rhs.indices_, size_ * sizeof(int));
	memcpy(values_, rhs.values_, size_ * sizeof(double));
    }
    else {
	indices_ = NULL;
	values_ = NULL;
    }
}

//#############################################################################
   
/** Pack Dco part into an encoded object. */
AlpsReturnStatus 
DcoConstraint::encodeDco(AlpsEncoded *encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    if (size_ < 0) {
	std::cout << "ERROR: encodeDco: size_=" << size_<<std::endl;
	assert(size_ > 0);
    }

    //std::cout << "encodeDco: constraint length/size_=" << size_<<std::endl;

    encoded->writeRep(indices_, size_);
    encoded->writeRep(values_, size_);
    return status;
}    

//#############################################################################

/** Unpack Dco part from a encode object. */
AlpsReturnStatus 
DcoConstraint::decodeDco(AlpsEncoded &encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    encoded.readRep(indices_, size_);
    if (size_ < 0) {
	std::cout << "ERROR: decodeDco: con1, size_=" << size_<<std::endl;
	assert(size_ > 0);
    }
    //std::cout << "----- decodeDco: con1, size_=" << size_<<std::endl;

    encoded.readRep(values_, size_);
    if (size_ < 0) {
	std::cout << "ERROR: decodeDco: con2, size_=" << size_<<std::endl;
	assert(size_ > 0);
    }

    //std::cout << "----- decodeDco: con2, size_=" << size_<<std::endl;
    return status;
}

//#############################################################################

AlpsReturnStatus 
DcoConstraint::encode(AlpsEncoded *encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    status = encodeBcpsObject(encoded);
    status = encodeDco(encoded);
    return status;
}

//#############################################################################

/** Decode a constraint from an encoded object. */
AlpsKnowledge* 
DcoConstraint::decode(AlpsEncoded& encoded) const 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    DcoConstraint* con = new DcoConstraint();    
    
    // Unpack Bcps object part.
    status = con->decodeBcpsObject(encoded);
    if (status) {
	throw CoinError("Failed to decode Bcps part",
			"decode", 
			"DcoObject");
    }
    
    // Unpack Dco part.
    status = con->decodeDco(encoded);
    if (status) {
	throw CoinError("Failed to decode Dco part", 
			"decode", 
			"DcoObject");
    }
    
    return con;
}

//#############################################################################

/** Compute hash value. */
void 
DcoConstraint::hashing(BcpsModel *model)
{
    assert(model != NULL);
    DcoModel *m = dynamic_cast<DcoModel *>(model);
    
    int k, ind;
    const double * randoms = m->getConRandoms();

    hashValue_ = 0.0;    
    for (k = 0; k < size_; ++k) {
	ind = indices_[k];
	hashValue_ += randoms[ind] * ind;
    }
#ifdef DISCO_DEBUG_MORE
    std::cout << "hashValue_=" << hashValue_ << std::endl;
#endif
}

//#############################################################################

OsiRowCut *
DcoConstraint::createOsiRowCut()
{
    double lower = CoinMax(getLbHard(), getLbSoft());
    double upper = CoinMin(getUbHard(), getUbSoft());
    
    OsiRowCut * cut = new OsiRowCut;
    if (!cut) {
        /* Out of memory. */
	throw CoinError("Out of Memory", "Dco_constraintToOsiCut", "NONE");
    }
    
    assert(size_ > 0);
    
    cut->setLb(lower);
    cut->setUb(upper);
    cut->setRow(size_, indices_, values_);
    
    return cut;
}

//#############################################################################

double
DcoConstraint::violation(const double *lpSolution)
{
    int k, varInd;
    double activity = 0.0;
    double rowLower = CoinMax(lbHard_, lbSoft_);
    double rowUpper = CoinMin(ubHard_, ubSoft_);
    double violation = -ALPS_DBL_MAX; // Any negative number is OK

    for (k = 0; k < size_; ++k) {
	varInd = indices_[k];
	activity += values_[varInd] * lpSolution[varInd];
    }
    
    if (rowLower > -ALPS_INFINITY) {
	violation = rowLower - activity;
    }
    if (rowUpper < ALPS_INFINITY) {
	violation = CoinMax(violation, activity - rowUpper);
    }
                    
    return violation;
}

//#############################################################################
