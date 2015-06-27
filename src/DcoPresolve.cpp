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

#include "DcoPresolve.hpp"

//#############################################################################

OsiConicSolverInterface *
DcoPresolve::preprocess(OsiConicSolverInterface & origModel,
                         double feasibilityTolerance,
                         bool keepIntegers,
                         int numberPasses,
                         const char * prohibited)
{
  OsiSolverInterface * om = dynamic_cast<OsiSolverInterface*>(&origModel);
  OsiSolverInterface * si = presolvedModel(*om,
					   feasibilityTolerance,
					   keepIntegers,
					   numberPasses,
					   prohibited);
  OsiConicSolverInterface * si_conic = dynamic_cast<OsiConicSolverInterface*>(si);
  return si_conic;
}

//#############################################################################

void 
DcoPresolve::postprocess(bool updateStatus)
{
    postsolve(updateStatus);
}

//#############################################################################
