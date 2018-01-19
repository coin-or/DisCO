/*===========================================================================*
 * This file is part of the Discrete Conic Optimization (DisCO) Solver.      *
 *                                                                           *
 * DisCO is distributed under the Eclipse Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *          Aykut Bulut, Lehigh University                                   *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and      *
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#include "DcoParams.hpp"

using std::make_pair;

void DcoParams::createKeywordList() {
  //--------------------------------------------------------
  // Create the list of keywords for parameter file reading
  //--------------------------------------------------------
  //--------------------------------------------------------
  // CharPar
  //--------------------------------------------------------
  keys_.push_back(make_pair(std::string("Dco_cutRampUp"),
                            AlpsParameter(AlpsBoolPar, cutRampUp)));
  keys_.push_back(make_pair(std::string("Dco_presolve"),
                            AlpsParameter(AlpsBoolPar, presolve)));
  keys_.push_back(make_pair(std::string("Dco_shareConstraints"),
                            AlpsParameter(AlpsBoolPar,
                                          shareConstraints)));
  keys_.push_back(make_pair(std::string("Dco_shareVariables"),
                            AlpsParameter(AlpsBoolPar,
                                          shareVariables)));
  keys_.push_back(make_pair(std::string("Dco_sharePseudocostRampUp"),
                            AlpsParameter(AlpsBoolPar,
                                          sharePseudocostRampUp)));
  keys_.push_back(make_pair(std::string("Dco_sharePseudocostSearch"),
                            AlpsParameter(AlpsBoolPar,
                                          sharePseudocostSearch)));
  //--------------------------------------------------------
  // BoolArrayPar
  //--------------------------------------------------------
  //--------------------------------------------------------
  // Int Parameters
  //--------------------------------------------------------
  keys_.push_back(make_pair(std::string("Dco_branchStrategy"),
                            AlpsParameter(AlpsIntPar, branchStrategy)));
  keys_.push_back(make_pair(std::string("Dco_branchStrategyRampUp"),
                            AlpsParameter(AlpsIntPar, branchStrategyRampUp)));
  keys_.push_back(make_pair(std::string("Dco_cutDisable"),
                            AlpsParameter(AlpsIntPar, cutDisable)));
  keys_.push_back(make_pair(std::string("Dco_cutStrategy"),
                            AlpsParameter(AlpsIntPar, cutStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutGenerationFrequency"),
                            AlpsParameter(AlpsIntPar,
                                          cutGenerationFrequency)));
  keys_.push_back(make_pair(std::string("Dco_cutCliqueStrategy"),
                            AlpsParameter(AlpsIntPar, cutCliqueStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutGomoryStrategy"),
                            AlpsParameter(AlpsIntPar, cutGomoryStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutFlowCoverStrategy"),
                            AlpsParameter(AlpsIntPar, cutFlowCoverStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutKnapsackStrategy"),
                            AlpsParameter(AlpsIntPar, cutKnapsackStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutMirStrategy"),
                            AlpsParameter(AlpsIntPar, cutMirStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutOddHoleStrategy"),
                            AlpsParameter(AlpsIntPar, cutOddHoleStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutProbingStrategy"),
                            AlpsParameter(AlpsIntPar, cutProbingStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutTwoMirStrategy"),
                            AlpsParameter(AlpsIntPar, cutTwoMirStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutIpmStrategy"),
                            AlpsParameter(AlpsIntPar, cutIpmStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutIpmIntStrategy"),
                            AlpsParameter(AlpsIntPar, cutIpmIntStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutOaStrategy"),
                            AlpsParameter(AlpsIntPar, cutOaStrategy)));
  keys_.push_back(make_pair(std::string("Dco_cutOaAlpha"),
                            AlpsParameter(AlpsIntPar, cutOaAlpha)));
  keys_.push_back(make_pair(std::string("Dco_cutOaGamma"),
                            AlpsParameter(AlpsIntPar, cutOaGamma)));
  keys_.push_back(make_pair(std::string("Dco_cutOaSlackLimit"),
                            AlpsParameter(AlpsIntPar, cutOaSlackLimit)));
  /// MILP Auto cut generation strategy parameters
  keys_.push_back(make_pair(std::string("Dco_cutMilpGamma"),
                            AlpsParameter(AlpsIntPar, cutMilpGamma)));

  keys_.push_back(make_pair(std::string("Dco_cutCliqueFreq"),
                            AlpsParameter(AlpsIntPar, cutCliqueFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutGomoryFreq"),
                            AlpsParameter(AlpsIntPar, cutGomoryFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutFlowCoverFreq"),
                            AlpsParameter(AlpsIntPar, cutFlowCoverFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutKnapsackFreq"),
                            AlpsParameter(AlpsIntPar, cutKnapsackFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutMirFreq"),
                            AlpsParameter(AlpsIntPar, cutMirFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutOddHoleFreq"),
                            AlpsParameter(AlpsIntPar, cutOddHoleFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutProbingFreq"),
                            AlpsParameter(AlpsIntPar, cutProbingFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutTwoMirFreq"),
                            AlpsParameter(AlpsIntPar, cutTwoMirFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutIpmFreq"),
                            AlpsParameter(AlpsIntPar, cutIpmFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutIpmIntFreq"),
                            AlpsParameter(AlpsIntPar, cutIpmIntFreq)));
  keys_.push_back(make_pair(std::string("Dco_cutOaFreq"),
                            AlpsParameter(AlpsIntPar, cutOaFreq)));
  keys_.push_back(make_pair(std::string("Dco_difference"),
                            AlpsParameter(AlpsIntPar, difference)));
  keys_.push_back(make_pair(std::string("Dco_heurStrategy"),
                            AlpsParameter(AlpsIntPar, heurStrategy)));
  keys_.push_back(make_pair(std::string("Dco_heurCallFrequencyy"),
                            AlpsParameter(AlpsIntPar, heurCallFrequency)));
  keys_.push_back(make_pair(std::string("Dco_heurRoundStrategy"),
                            AlpsParameter(AlpsIntPar, heurRoundStrategy)));
  keys_.push_back(make_pair(std::string("Dco_heurRoundFreq"),
                            AlpsParameter(AlpsIntPar, heurRoundFreq)));
  keys_.push_back(make_pair(std::string("Dco_lookAhead"),
                            AlpsParameter(AlpsIntPar, lookAhead)));
  keys_.push_back(make_pair(std::string("Dco_pseudoReliability"),
                            AlpsParameter(AlpsIntPar, pseudoReliability)));
  keys_.push_back(make_pair(std::string("Dco_sharePcostDepth"),
                            AlpsParameter(AlpsIntPar, sharePcostDepth)));
  keys_.push_back(make_pair(std::string("Dco_sharePcostFrequency"),
                            AlpsParameter(AlpsIntPar, sharePcostFrequency)));
  keys_.push_back(make_pair(std::string("Dco_strongCandSize"),
                            AlpsParameter(AlpsIntPar, strongCandSize)));
  // conic cut related
  // keys_.push_back(make_pair(std::string("Dco_conicCutStrategy"),
  //                           AlpsParameter(AlpsIntPar, conicCutStrategy)));
  // keys_.push_back(make_pair(std::string("Dco_conicCutGenerationFrequency"),
  //                           AlpsParameter(AlpsIntPar, conicCutGenerationFrequency)));
  // keys_.push_back(make_pair(std::string("Dco_conicCutMirStrategy"),
  //                           AlpsParameter(AlpsIntPar, conicCutMirStrategy)));
  // keys_.push_back(make_pair(std::string("Dco_conicCutGD1Strategy"),
  //                           AlpsParameter(AlpsIntPar, conicCutGD1Strategy)));
  // keys_.push_back(make_pair(std::string("Dco_conicCutGD2Strategy"),
  //                           AlpsParameter(AlpsIntPar, conicCutGD2Strategy)));
  // keys_.push_back(make_pair(std::string("Dco_conicCutMirFreq"),
  //                           AlpsParameter(AlpsIntPar, conicCutMirFreq)));
  // keys_.push_back(make_pair(std::string("Dco_conicCutGD1Freq"),
  //                           AlpsParameter(AlpsIntPar, conicCutGD1Freq)));
  // keys_.push_back(make_pair(std::string("Dco_conicCutGD2Freq"),
  //                           AlpsParameter(AlpsIntPar, conicCutGD2Freq)));
  keys_.push_back(make_pair(std::string("Dco_logLevel"),
                            AlpsParameter(AlpsIntPar, logLevel)));
  keys_.push_back(make_pair(std::string("Dco_presolveNumPass"),
                            AlpsParameter(AlpsIntPar, presolveNumPass)));
  keys_.push_back(make_pair(std::string("Dco_approxNumPass"),
                            AlpsParameter(AlpsIntPar, approxNumPass)));
  //--------------------------------------------------------
  // Double Parameters.
  //--------------------------------------------------------
  keys_.push_back(make_pair(std::string("Dco_cutFactor"),
                            AlpsParameter(AlpsDoublePar, cutFactor)));
  keys_.push_back(make_pair(std::string("Dco_cutoff"),
                            AlpsParameter(AlpsDoublePar, cutoff)));
  keys_.push_back(make_pair(std::string("Dco_objTol"),
                            AlpsParameter(AlpsDoublePar, objTol)));
  keys_.push_back(make_pair(std::string("Dco_denseConFactor"),
                            AlpsParameter(AlpsDoublePar, denseConFactor)));
  keys_.push_back(make_pair(std::string("Dco_integerTol"),
                            AlpsParameter(AlpsDoublePar, integerTol)));
  keys_.push_back(make_pair(std::string("Dco_coneTol"),
                            AlpsParameter(AlpsDoublePar, coneTol)));
  keys_.push_back(make_pair(std::string("Dco_objSense"),
                            AlpsParameter(AlpsDoublePar, objSense)));
  keys_.push_back(make_pair(std::string("Dco_optimalRelGap"),
                            AlpsParameter(AlpsDoublePar, optimalRelGap)));
  keys_.push_back(make_pair(std::string("Dco_optimalAbsGap"),
                            AlpsParameter(AlpsDoublePar, optimalAbsGap)));
  keys_.push_back(make_pair(std::string("Dco_pseudoWeight"),
                            AlpsParameter(AlpsDoublePar, pseudoWeight)));
  keys_.push_back(make_pair(std::string("Dco_scaleConFactor"),
                            AlpsParameter(AlpsDoublePar, scaleConFactor)));
  keys_.push_back(make_pair(std::string("Dco_tailOff"),
                            AlpsParameter(AlpsDoublePar, tailOff)));
  keys_.push_back(make_pair(std::string("Dco_presolveTolerance"),
                            AlpsParameter(AlpsDoublePar, presolveTolerance)));
  keys_.push_back(make_pair(std::string("Dco_approxFactor"),
                            AlpsParameter(AlpsDoublePar, approxFactor)));
  keys_.push_back(make_pair(std::string("Dco_cutOaBeta"),
                            AlpsParameter(AlpsDoublePar, cutOaBeta)));
  keys_.push_back(make_pair(std::string("Dco_cutOaSlack1"),
                            AlpsParameter(AlpsDoublePar, cutOaSlack1)));
  keys_.push_back(make_pair(std::string("Dco_cutOaSlack2"),
                            AlpsParameter(AlpsDoublePar, cutOaSlack2)));
  keys_.push_back(make_pair(std::string("Dco_cutMilpDelta"),
                            AlpsParameter(AlpsDoublePar, cutMilpDelta)));
//--------------------------------------------------------
  // String Parameters.
  //--------------------------------------------------------
}

//#############################################################################

void DcoParams::setDefaultEntries() {
  //-------------------------------------------------------------
  // Char Parameters.
  //-------------------------------------------------------------
  setEntry(cutRampUp, true);
  setEntry(presolve, false);
  setEntry(shareConstraints, false);
  setEntry(shareVariables, false);
  setEntry(sharePseudocostRampUp, true);
  setEntry(sharePseudocostSearch, false);
  // presolve parameters
  setEntry(presolveKeepIntegers, true);
  setEntry(presolveTransform, true);
  //-------------------------------------------------------------
  // Int Parameters.
  //-------------------------------------------------------------
  setEntry(branchStrategy, DcoBranchingStrategyPseudoCost);
  setEntry(branchStrategyRampUp, DcoBranchingStrategyPseudoCost);
  setEntry(cutStrategy, DcoCutStrategyNotSet);
  setEntry(cutGenerationFrequency, 1);
  setEntry(cutDisable, 20);
  setEntry(cutCliqueStrategy, DcoCutStrategyAuto);
  setEntry(cutGomoryStrategy, DcoCutStrategyAuto);
  setEntry(cutFlowCoverStrategy, DcoCutStrategyAuto);
  setEntry(cutKnapsackStrategy, DcoCutStrategyAuto);
  setEntry(cutMirStrategy, DcoCutStrategyAuto);
  setEntry(cutOddHoleStrategy, DcoCutStrategyAuto);
  setEntry(cutProbingStrategy, DcoCutStrategyNone);
  setEntry(cutTwoMirStrategy, DcoCutStrategyAuto);
  setEntry(cutIpmStrategy, DcoCutStrategyNotSet);
  setEntry(cutIpmIntStrategy, DcoCutStrategyNotSet);
  setEntry(cutOaStrategy, DcoCutStrategyNotSet);
  /// OA cut strategy parameters
  setEntry(cutOaAlpha, 1);
  setEntry(cutOaGamma, 50);
  setEntry(cutOaSlackLimit, 3);
  /// MILP cut generation strategy parameters
  setEntry(cutMilpGamma, 20);

  setEntry(cutCliqueFreq, 100);
  setEntry(cutGomoryFreq, 100);
  setEntry(cutFlowCoverFreq, 100);
  setEntry(cutKnapsackFreq, 100);
  setEntry(cutMirFreq, 100);
  setEntry(cutOddHoleFreq, 100);
  setEntry(cutProbingFreq, 100);
  setEntry(cutTwoMirFreq, 100);
  setEntry(cutIpmFreq, 1);
  setEntry(cutIpmIntFreq, 1);
  setEntry(cutOaFreq, 1);
  setEntry(difference, -1);
  setEntry(heurStrategy, DcoHeurStrategyPeriodic);
  setEntry(heurCallFrequency, 1);
  setEntry(heurRoundStrategy, DcoHeurStrategyPeriodic);
  setEntry(heurRoundFreq, 100);
  setEntry(lookAhead, 4);
  setEntry(pseudoReliability, 8);
  setEntry(sharePcostDepth, 30);
  setEntry(sharePcostFrequency, 100);
  setEntry(strongCandSize, 1000);
  setEntry(logLevel, 2);
  setEntry(presolveNumPass, 5);
  setEntry(approxNumPass, 400);
  //-------------------------------------------------------------
  // Double Parameters
  //-------------------------------------------------------------
  setEntry(cutFactor, 4.0);
  setEntry(cutoff, ALPS_INC_MAX);
  setEntry(objTol, 1.0e-6);
  setEntry(denseConFactor, 5.0);
  setEntry(integerTol, 1.0e-5);
  setEntry(coneTol, 1.0e-5);
  setEntry(objSense, 1.0);
  setEntry(optimalRelGap, 1.0e-4);
  setEntry(optimalAbsGap, 1.0e-6);
  setEntry(pseudoWeight, 0.8);
  setEntry(scaleConFactor, 1000.0);
  setEntry(tailOff, 1e-8);
  setEntry(presolveTolerance, 0.0);
  // approximation factor, used in OA
  setEntry(approxFactor, 1.0);
  // threshold for cut activity used in approximateCones()
  setEntry(cutOaSlack1, 0.0001);
  // threshold for cut activity used in bounding loop
  setEntry(cutOaSlack2, 0.005);
  setEntry(cutOaBeta, 0.001);
  /// MILP cut generation strategy parameters
  setEntry(cutMilpDelta, 0.0001);
  //-------------------------------------------------------------
  // String Parameters
  //-------------------------------------------------------------
}
