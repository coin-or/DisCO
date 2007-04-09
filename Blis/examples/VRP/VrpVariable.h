/*===========================================================================*
 * This file is part of a solver for the Vehicle Routing Problem             *
 * developed using the BiCePS Linear Integer Solver (BLIS).                  *
 *                                                                           *
 * This solver is distributed under the Common Public License as part of     * 
 * the COIN-OR repository (http://www.coin-or.org).                          *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2007 Yan Xu and Ted Ralphs.                                 *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef VrpVariable_h_
#define VrpVariable_h_

//#############################################################################

#include "BlisVariable.h"

//#############################################################################

/** Variable class for VRP. */
class VrpVariable : public BlisVariable 
{    
private:

   /* The endpoints of the edge */
   int ends[2];
   double objCoef_;
   int size_;
   int *indices_;
   double *values_;
   
public:

    VrpVariable() : objCoef_(0.0), size_(0), indices_(NULL), values_(NULL) {
       ends[0] = 0;
       ends[1] = 0;
    }

    VrpVariable(int v1, int v2, int cost) {
       ends[0] = (v1 < v2) : v1 ? v2;
       ends[1] = (v1 < v2) : v2 ? v1;
       size_ = 2;
       indices_ = new int[size_];
       values_ = new int[size_];
       indices_[0] = ends[0];
       indices_[1] = ends[1];
       values_[0] = values_[1] = 1.0;
       intType_ = 'B';
       lbHard_ = 0.0;
       ubHard_ = 1.0;
    }

    virtual ~VrpVariable(){ 
       if (size_ > 0) {
	  delete [] indices_; indices_ = NULL;
	  delete [] values_; values_ = NULL;
       }
    
};

//#############################################################################

#endif
