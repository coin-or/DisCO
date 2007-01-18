/*===========================================================================*
 * This file is part of the Abstract Library for Parallel Search (ALPS).     *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "BlisParams.h"

using std::make_pair;

void 
BlisParams::createKeywordList() {

  //--------------------------------------------------------
  // Create the list of keywords for parameter file reading
  //--------------------------------------------------------

  //--------------------------------------------------------
  // CharPar
  //--------------------------------------------------------
 
   keys_.push_back(make_pair(std::string("Blis_cutDuringRampup"),
			     AlpsParameter(AlpsCharPar, cutDuringRampup)));
   
   keys_.push_back(make_pair(std::string("Blis_presolve"),
			     AlpsParameter(AlpsCharPar, presolve)));
    
   keys_.push_back(make_pair(std::string("Blis_sharePseudocostRampup"),
			     AlpsParameter(AlpsCharPar,
					   sharePseudocostRampup)));
    
   keys_.push_back(make_pair(std::string("Blis_sharePseudocostSearch"),
			     AlpsParameter(AlpsCharPar,
					   sharePseudocostSearch)));
    
  //--------------------------------------------------------
  // BoolArrayPar
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  // Int Parameters
  //--------------------------------------------------------

  keys_.push_back(make_pair(std::string("Blis_branchStrategy"),
			    AlpsParameter(AlpsIntPar, branchStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutClique"),
			    AlpsParameter(AlpsIntPar, cutClique)));

  keys_.push_back(make_pair(std::string("Blis_cutGomory"),
			    AlpsParameter(AlpsIntPar, cutGomory)));

  keys_.push_back(make_pair(std::string("Blis_cutFlowCover"),
			    AlpsParameter(AlpsIntPar, cutFlowCover)));

  keys_.push_back(make_pair(std::string("Blis_cutKnapsack"),
			    AlpsParameter(AlpsIntPar, cutKnapsack)));
  
  keys_.push_back(make_pair(std::string("Blis_cutMir"),
			    AlpsParameter(AlpsIntPar, cutMir)));
  
  keys_.push_back(make_pair(std::string("Blis_cutOddHole"),
			    AlpsParameter(AlpsIntPar, cutOddHole)));
  
  keys_.push_back(make_pair(std::string("Blis_cutProbing"),
			    AlpsParameter(AlpsIntPar, cutProbing)));

  keys_.push_back(make_pair(std::string("Blis_cutStrategy"),
			    AlpsParameter(AlpsIntPar, cutStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutTwoMir"),
			    AlpsParameter(AlpsIntPar, cutTwoMir)));

  keys_.push_back(make_pair(std::string("Blis_difference"),
			    AlpsParameter(AlpsIntPar, difference)));
  
  keys_.push_back(make_pair(std::string("Blis_heurRound"),
			    AlpsParameter(AlpsIntPar, heurRound)));
  
  keys_.push_back(make_pair(std::string("Blis_heurStrategy"),
                            AlpsParameter(AlpsIntPar, heurStrategy)));
  
  keys_.push_back(make_pair(std::string("Blis_lookAhead"),
			    AlpsParameter(AlpsIntPar, lookAhead)));
  
  keys_.push_back(make_pair(std::string("Blis_pseudoRelibility"),
			    AlpsParameter(AlpsIntPar, pseudoRelibility)));

  keys_.push_back(make_pair(std::string("Blis_rampUpBranchStrategy"),
			    AlpsParameter(AlpsIntPar, rampUpBranchStrategy)));

  keys_.push_back(make_pair(std::string("Blis_strongCandSize"),
			    AlpsParameter(AlpsIntPar, strongCandSize)));
  
  //--------------------------------------------------------
  // Double Parameters.
  //--------------------------------------------------------
  
  keys_.push_back(make_pair(std::string("Blis_cutFactor"),
			    AlpsParameter(AlpsDoublePar, cutFactor)));
  
  keys_.push_back(make_pair(std::string("Blis_cutoffInc"),
			    AlpsParameter(AlpsDoublePar, cutoffInc)));
  
  keys_.push_back(make_pair(std::string("Blis_denseConFactor"),
			    AlpsParameter(AlpsDoublePar, denseConFactor)));

  keys_.push_back(make_pair(std::string("Blis_integerTol"),
			    AlpsParameter(AlpsDoublePar, integerTol)));
  
  keys_.push_back(make_pair(std::string("Blis_optimalRelGap"),
			    AlpsParameter(AlpsDoublePar, optimalRelGap)));
  
  keys_.push_back(make_pair(std::string("Blis_optimalAbsGap"),
			    AlpsParameter(AlpsDoublePar, optimalAbsGap)));
  
  keys_.push_back(make_pair(std::string("Blis_pseudoWeight"),
			    AlpsParameter(AlpsDoublePar, pseudoWeight)));
  
  keys_.push_back(make_pair(std::string("Blis_scaleConFactor"),
			    AlpsParameter(AlpsDoublePar, scaleConFactor)));
  
  //--------------------------------------------------------
  // String Parameters.
  //--------------------------------------------------------

}

//#############################################################################

void 
BlisParams::setDefaultEntries() {

  //-------------------------------------------------------------
  // Char Parameters.
  //-------------------------------------------------------------

  setEntry(cutDuringRampup, false);
  setEntry(presolve, false);
  setEntry(sharePseudocostRampup, true);
  setEntry(sharePseudocostSearch, false);

  //-------------------------------------------------------------
  // Int Parameters.
  //-------------------------------------------------------------

  setEntry(branchStrategy, 1);
  setEntry(cutClique, BLIS_NOT_SET);
  setEntry(cutGomory, BLIS_NOT_SET);
  setEntry(cutFlowCover, BLIS_NOT_SET);
  setEntry(cutKnapsack, BLIS_NOT_SET);
  setEntry(cutMir, BLIS_NOT_SET);
  setEntry(cutOddHole, BLIS_NOT_SET);
  setEntry(cutProbing, BLIS_NOT_SET);
  setEntry(cutStrategy, BLIS_AUTO);
  setEntry(cutTwoMir, BLIS_NOT_SET);
  setEntry(difference, -1);
  setEntry(heurRound, BLIS_NOT_SET);
  setEntry(heurStrategy, -1);
  setEntry(lookAhead, 4);
  setEntry(pseudoRelibility, 8);
  setEntry(rampUpBranchStrategy, BLIS_BS_PSEUDOCOST);
  setEntry(strongCandSize, 10);
  
  //-------------------------------------------------------------
  // Double Parameters
  //-------------------------------------------------------------

  setEntry(cutFactor, 4.0);
  setEntry(cutoffInc, 1.0e-6);
  setEntry(denseConFactor, 5.0);
  setEntry(integerTol, 1.0e-5);
  setEntry(optimalRelGap, 1.0e-6);
  setEntry(optimalAbsGap, 1.0e-4);
  setEntry(pseudoWeight, 1.0);
  setEntry(scaleConFactor, 1000000.0);

  //-------------------------------------------------------------
  // String Parameters
  //-------------------------------------------------------------
  
}
