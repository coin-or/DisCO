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

#ifndef VrpCutGenerator_h_
#define VrpCutGenerator_h_

//#############################################################################

#include "BlisConGenerator.h"
#include "VrpModel.h"
#include "VrpNetwork.h"
#include "VrpMacros.h"

class VrpModel;

//#############################################################################

class VrpCutGenerator : public BlisConGenerator  
{
private:

   VrpModel *model_;
   char **coef_list;
   int *ref_;
   double *cutVal_;
   char *cutList_;
   char *inSet_;
   VrpNetwork *n_;
   
public:
    
    /** Default construtor. */
   VrpCutGenerator();
   
   VrpCutGenerator(int edgenum, int vertnum); 
   
   /** Destructor. */
   virtual ~VrpCutGenerator() {
      delete [] ref_; ref_ = 0;
      delete [] cutVal_; cutVal_ = 0;
      delete [] cutList_; cutList_ = 0;
      delete [] inSet_; inSet_ = 0;
      delete [] n_; n_ = 0;
   }
   
   /** Generate cons for the client model.
       
   Evaluate the state of the client model and decide whether to 
   generate cons. The generated cons are inserted into and returned 
   in the collection of cons \p cs.
   
   If \p fullScan is true, the generator is obliged to call the CGL
   \c generateCuts routine.  Otherwise, it is free to make a local 
   decision. The current implementation uses \c strategy_ to decide.
   
   The routine returns true if reoptimisation is needed (because the 
   state of the solver interface has been modified).
   */

   virtual bool generateCons(OsiCuts &cs, bool fullScan);
  
   CoinPackedVector *getSolution();

   void createNet(CoinPackedVector *vec);
   
   bool connectivityCuts(OsiCuts &cs);
   
   bool addCut(OsiCuts &cs, char *coef, int rhs, int type);
};

//#############################################################################

#endif
