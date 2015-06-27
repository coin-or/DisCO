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

#ifndef DcoConstraint_h_
#define DcoConstraint_h_

#include "BcpsObject.h"

class OsiRowCut;

//#############################################################################

class DcoConstraint : public BcpsConstraint {
    
 protected:
    /** Number of nonzero coefficients */
    int size_;
    /** Variable indices */
    int *indices_;
    /** Value of nonzero coefficients */
    double *values_;

 public:

    /** Default constructor. */
    DcoConstraint();

    /** Useful constructor. */
    DcoConstraint(int s, const int *ind, const double *val);

    /** Useful constructor. */
    DcoConstraint(double lbh, double ubh, double lbs, double ubs);
    
    /** Useful constructor. */
    DcoConstraint(double lbh, double ubh, double lbs, double ubs,
                   int size, const int *ind, const double *val);
    /** Destructor. */
    virtual ~DcoConstraint();

    /** Copy constructor. */
    DcoConstraint(const DcoConstraint & rhs);
    
    /** Return data  */
    /**@{*/
    int getSize() const       { return size_; }
    int* getIndices() const   { return indices_; }
    double* getValues() const { return values_; }    
    /**@}*/
    
    /** Set data  */
    /**@{*/
    void setData(int s, const int *ind, const double *val) {
	if (size_ < s) {
	    delete [] indices_;
	    delete [] values_;
	    indices_ = new int [s];
	    values_ = new double [s];
	}
	size_ = s;
	memcpy(indices_, ind, sizeof(int) * s);
	memcpy(values_, val, sizeof(double) * s);
    }
    /**@}*/

 protected:

    /** Pack Dco part into an encoded object. */
    AlpsReturnStatus encodeDco(AlpsEncoded *encoded);

    /** Unpack Dco part from a encode object. */
    AlpsReturnStatus decodeDco(AlpsEncoded &encoded);
	    
 public:

    /** Create a OsiRowCut based on this constraint. */
    OsiRowCut *createOsiRowCut();

    /** Compute a hash key. */
    virtual void hashing(BcpsModel *model=NULL);

    /** Check if violates a given lp solution. */
    double violation(const double *lpSolution);

    using AlpsKnowledge::encode ;
    /** Pack into a encode object. */
    virtual AlpsReturnStatus encode(AlpsEncoded *encoded);
    
    /** Decode a constraint from an encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded& encoded) const;
};

//#############################################################################

#endif
