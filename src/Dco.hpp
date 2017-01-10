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
 * Copyright (C) 2001-2017, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef Dco_hpp_
#define Dco_hpp_

#include "AlpsConfig.h"
#include "BcpsConfig.h"
#include "DcoConfig.hpp"

/*!
  /mainpage

  # Bcps ideas for future

  We keep relaxed cols/rows (the whole object or integrality of the object) in
  BcpsModel in list BcpsObject ** relaxed_. This way Bcps can check the
  feasibility of its relaxed objects through virtual functions that will be
  impelemnted in Disco level.

  Can subproblems have different relaxed objects? Yes they can. But make sure
  all relaxed objects are in the relaxed_. If subproblems have different
  relaxed objects we will get a more depocposition like algorithm. Moreover if
  not all integrality is relaxed we might need a discrete solver (a DcoModel
  instance?) to solve the problem.

  Subproblems might have different solvers? A subprobllem might be an LP, SOCO,
  MILP or even MISOCO. How will DcoTreeNode::bound know about the solver to be
  used?

  When subproblem is solved with the solver we go back and check the
  feasibility of the objects that are relaxed in the subproblem (objects in
  relaxed_). This is done through infeasible() virtual function defined in
  BcpsObject.
  After this check we have a list of objects that are infeasible. At this
  point we need to decide what to do with them. Options are (1) generate
  cuts (using generateConstraints()), (2) lift the subproblem (using
  generateVariables()), (3) branch.
  (1) We can generate cuts when infeasible objects are continuous or
  integer. Generate a separating hyperplane that cuts the subproblem solution
  from feasible region.
  (2) Branching when the object is continuous. This is similar to branch
  on variables. Branching procedure should create new subproblems where
  infeasible object has new upper and lower bounds.

  feasibility checking should be implemented in Bcps level. Itertate over
  cols/rows and check their feasiblity. Store infeasible cols (BcpsVariable)
  and rows (BcpsConstraint). This function checks the feasibility of the
  solution stored in the solver at the time of the call.

  BcpsModel::feasibleSolution() should be able to check cols or rows only.
  In a typical branch and bound we need to check feasibility of cols only.
  In DisCO we may want to check both or check cols only.

  # DisCO ideas

  # Style guides

  <ul>

    <li> If a function does not fit in the same line of its decleration in the
    header file then it is not supposed to be there. Move it to .cpp file.

    <li> All functions defined in .hpp are already inline since bla.bla version
    of gcc. This makes inline key-words in the header files redundants.  I will
    just remove them, since I can not bear redundant stuff floating around.

    <li> Define pointers as "Type * name". * is not next to Type or variable
    name. * is separate since it is neither part of the Type nor variable
    name. This is the way recommended by Bjarne Stroustrup and it makes sense.

    <il> We use const specifiers as suggested by Bjarne Stroustrup in his book.
    Check the code to see how I declare const objects, check the book to see
    why.

  </ul>

*/

#define DISCO_INFINITY 1e30

enum DcoConstraintType {
  DcoConstraintTypeLinear = 0,
  DcoConstraintTypeConic
};

enum DcoLorentzConeType {
  DcoLorentzCone = 0,
  DcoRotatedLorentzCone
};

enum DcoSubproblemStatus{
  DcoSubproblemStatusOptimal,
  DcoSubproblemStatusAbandoned,
  DcoSubproblemStatusPrimalInfeasible,
  DcoSubproblemStatusDualInfeasible,
  DcoSubproblemStatusPrimalObjLim,
  DcoSubproblemStatusDualObjLim,
  DcoSubproblemStatusIterLim,
  DcoSubproblemStatusUnknown
};

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


enum DcoCutStrategy {
  DcoCutStrategyNotSet = -1,
  DcoCutStrategyNone = 0,
  DcoCutStrategyRoot,
  DcoCutStrategyAuto,
  DcoCutStrategyPeriodic
};

enum DcoConicCutStrategy {
  DcoConicCutStrategyNotSet = -1,
  DcoConicCutStrategyNone = 0,
  DcoConicCutStrategyRoot,
  DcoConicCutStrategyAuto,
  DcoConicCutStrategyPeriodic
};

enum DcoHeurStrategy {
  DcoHeurStrategyNotSet = -1,
  DcoHeurStrategyNone = 0,
  DcoHeurStrategyRoot,
  DcoHeurStrategyAuto,
  DcoHeurStrategyPeriodic,
  DcoHeurStrategyBeforeRoot // Before solving first relaxation
};

enum DcoHeurType {
  DcoHeurTypeNotSet = -1,
  DcoHeurTypeRounding
};

enum DcoHotStartStrategy{
  DcoHotStartBranchIncorrect,
  DcoHotStartBranchCorrect
};

enum DcoBranchingStrategy{
  DcoBranchingStrategyMaxInfeasibility,
  DcoBranchingStrategyPseudoCost,
  DcoBranchingStrategyReliability,
  DcoBranchingStrategyStrong,
  DcoBranchingStrategyBilevel
};

enum DcoSolutionType {
  DcoSolutionTypeBounding,
  DcoSolutionTypeBranching,
  DcoSolutionTypeDiving,
  DcoSolutionTypeHeuristic,
  DcoSolutionTypeStrong
};

/** Branching object type. */
enum DcoBranchingObjectType {
  DcoBranchingObjectTypeNone = 0,
  DcoBranchingObjectTypeInt,
  DcoBranchingObjectTypeSos,
  DcoBranchingObjectTypeBilevel
};

/** Node branch direction, is it a left node or right */
enum DcoNodeBranchDir {
  DcoNodeBranchDirectionDown = 0,
  DcoNodeBranchDirectionUp
};

/** Integral type */
enum DcoIntegralityType {
  DcoIntegralityTypeCont = 0,
  DcoIntegralityTypeInt
};

#endif
