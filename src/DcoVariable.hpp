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

#ifndef DcoVariable_h_
#define DcoVariable_h_

#include "BcpsObject.h"

//#############################################################################

class DcoVariable : public BcpsVariable {

 private:

    double objCoef_;
    int size_;
    int *indices_;
    double *values_;
    
 public:

    DcoVariable() : objCoef_(0.0), size_(0), indices_(NULL), values_(NULL) {}
    
    DcoVariable(double obj, int s, const int *ind, const double *val)
	{
            objCoef_ = obj;
	    size_ = s;
	    indices_ = new int [s];
	    values_ = new double [s];
	    memcpy(indices_, ind, s * sizeof(int));
	    memcpy(values_, val, s * sizeof(double));
	}
    
    DcoVariable(double lbh, double ubh, double lbs, double ubs) 
	:
	BcpsVariable(lbh, ubh, lbs, ubs),
        objCoef_(0.0),
	size_(0), indices_(NULL), values_(NULL)
	{}

    DcoVariable(double lbh, double ubh, double lbs, double ubs,
                 double obj, int s, const int *ind, const double *val)
        : 
        BcpsVariable(lbh, ubh, lbs, ubs)
        {
            objCoef_ = obj;
            size_ = s;
            indices_ = new int [s];
            values_ = new double [s];
            memcpy(indices_, ind, s * sizeof(int));
            memcpy(values_, val, s * sizeof(double));
        }
    
    virtual ~DcoVariable(){ 
        delete [] indices_; indices_ = NULL;
        delete [] values_; values_ = NULL;
    }
    
    /** Return data  */
    /**@{*/
    double getObjCoef()    { return objCoef_; }    
    int getSize() const     { return size_; }
    int* getIndices() const { return indices_; }
    double* getValues()     { return values_; }    
    /**@}*/
    
    /** Set data  */
    /**@{*/
    void setData(int s, const int *ind, const double *val) {
	if (size_ < s) {
	    delete [] indices_; indices_ = NULL;
	    delete [] values_; values_ = NULL;
	    indices_ = new int [s];
	    values_ = new double [s];
	}
	size_ = s;
	memcpy(indices_, ind, sizeof(int) * s);
	memcpy(values_, val, sizeof(double) * s);
    }
    void setObjCoef(double coef)    { objCoef_ = coef; }    
   /**@}*/

 protected:

   /** Pack Dco part into an encoded object. */
    AlpsReturnStatus encodeDco(AlpsEncoded *encoded) {
	AlpsReturnStatus status = AlpsReturnStatusOk;

        //std::cout << "****** encodeDco var: size_ = " << size_ << std::endl;

	encoded->writeRep(objCoef_);
	encoded->writeRep(indices_, size_);
	encoded->writeRep(values_, size_);
        
	return status;
    }    

    /** Unpack Dco part from a encode object. */
    AlpsReturnStatus decodeDco(AlpsEncoded &encoded) {
	AlpsReturnStatus status = AlpsReturnStatusOk;

	encoded.readRep(objCoef_);
	encoded.readRep(indices_, size_);
	encoded.readRep(values_, size_);
        
        //std::cout << "****** decodeDco var: size_ = " << size_ << std::endl;

	return status;
    }
    
 public:
    
    using AlpsKnowledge::encode ;
    /** Pack to a encode object. */
    virtual AlpsReturnStatus encode(AlpsEncoded *encoded){
	AlpsReturnStatus status;

	status = encodeBcpsObject(encoded);
	status = encodeDco(encoded);
	
	return status;
    }

    /** Decode a variable from an encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded &encoded) const {
	AlpsReturnStatus status = AlpsReturnStatusOk;
	DcoVariable * var = new DcoVariable();    
	
	// Unpack Bcps part.
	status = var->decodeBcpsObject(encoded);
	if (status) {
	    throw CoinError("Failed to decode Bcps part of var",
			    "decode",
			    "DcoObject");
	}
	
	// Unpack Dco part.
	status = var->decodeDco(encoded);
	if (status) {
	    throw CoinError("Failed to decode Dco part of var", 
			    "decode", 
			    "DcoObject");
	}
	return var;
    }
    
};

//#############################################################################

#endif /* End of file */

