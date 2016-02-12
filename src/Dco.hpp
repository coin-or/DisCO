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

#ifndef Dco_hpp_
#define Dco_hpp_

#include "AlpsConfig.h"
#include "BcpsConfig.h"
#include "DcoConfig.hpp"

//#############################################################################

#define DISCO_INFINITY 1e30

enum DcoLorentzConeType {
  DcoLorentzCone = 0,
  DcoRotatedLorentzCone
};

enum DcoLpStatus{
   DcoLpStatusOptimal,
   DcoLpStatusAbandoned,
   DcoLpStatusPrimalInfeasible,
   DcoLpStatusDualInfeasible,
   DcoLpStatusPrimalObjLim,
   DcoLpStatusDualObjLim,
   DcoLpStatusIterLim,
   DcoLpStatusUnknown
};

//#############################################################################

enum DcoReturnStatus {
   DcoReturnStatusOk = 0,
   DcoReturnStatusErrLp,
   DcoReturnStatusInfeasible,
   DcoReturnStatusUnbounded,
   DcoReturnStatusOverObjLim,
   DcoReturnStatusFeasible,
   DcoReturnStatusBranch,
   DcoReturnStatusUnknown
};

#if 0
#define DISCO_ERR_LP         100
#define DISCO_INF            200
#define DISCO_UNBOUND        201
#define DISCO_OPTIMAL          0
#define DISCO_UNKNOWN        202
#endif

//#############################################################################

enum DcoCutStrategy{
   DcoCutStrategyNotSet = -1,
   DcoCutStrategyNone = 0,
   DcoCutStrategyRoot,
   DcoCutStrategyAuto,
   DcoCutStrategyPeriodic
};

enum DcoConicCutStrategy{
   DcoConicCutStrategyNotSet = -1,
   DcoConicCutStrategyNone = 0,
   DcoConicCutStrategyRoot,
   DcoConicCutStrategyAuto,
   DcoConicCutStrategyPeriodic
};

enum DcoHeurStrategy{
   DcoHeurStrategyNotSet = -1,
   DcoHeurStrategyNone = 0,
   DcoHeurStrategyRoot,
   DcoHeurStrategyAuto,
   DcoHeurStrategyPeriodic,
   DcoHeurStrategyBeforeRoot // Before solving first relaxation
};

#if 0
#define DISCO_NOT_SET       -555
#define DISCO_ROOT            -2
#define DISCO_AUTO            -1
#define DISCO_NONE             0
#endif

//#############################################################################

enum DcoHotStartStrategy{
   DcoHotStartBranchIncorrect,
   DcoHotStartBranchCorrect
};

//#############################################################################

enum DcoBranchingStrategy{
   DcoBranchingStrategyMaxInfeasibility,
   DcoBranchingStrategyPseudoCost,
   DcoBranchingStrategyReliability,
   DcoBranchingStrategyStrong,
   DcoBranchingStrategyBilevel
};

//#############################################################################

enum DcoSolutionType {
    DcoSolutionTypeBounding,
    DcoSolutionTypeBranching,
    DcoSolutionTypeDiving,
    DcoSolutionTypeHeuristic,
    DcoSolutionTypeStrong
};

//#############################################################################

/** Branching object type. */
enum DcoBranchingObjectType {
    DcoBranchingObjectTypeNone = 0,
    DcoBranchingObjectTypeInt,
    DcoBranchingObjectTypeSos,
    DcoBranchingObjectTypeBilevel
};

/** Node branch direction, is it a left node or right */
enum DcoNodeBranchDir {
    DcoNodeBranchDirectionLeft = 0,
    DcoNodeBranchDirectionRight
};

/** Integral type */
enum DcoIntegralityType {
    DcoIntegralityTypeCont = 0,
    DcoIntegralityTypeInt
};


//#############################################################################

#define DISCO_CUT_DISABLE            20

#define DISCO_HEUR_ROUND_DISABLE     1000000

#define DISCO_PSEUDO                 21

//#############################################################################

#endif
