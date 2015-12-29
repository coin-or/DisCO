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


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "CoinTime.hpp"
#include "OsiConicSolverInterface.hpp"
#include "CglProbing.hpp"

#include "BcpsObjectPool.h"

#include "DcoHelp.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoConGenerator.hpp"
#include "DcoConstraint.hpp"

//#############################################################################

// Copy constructor
DcoConGenerator::DcoConGenerator(const DcoConGenerator & rhs): DcoConGeneratorBase(rhs) {
  setGenerator(rhs.generator());
}

//#############################################################################
// Assignment operator
DcoConGenerator &
DcoConGenerator::operator=( const DcoConGenerator& rhs) {
  if (this != &rhs) {
    DcoConGeneratorBase::operator=(rhs);
    setGenerator(rhs.generator());
  }
  return *this;
}

//#############################################################################

// Generate constraints for the model data contained in si.
// The generated constraints are inserted into and returned in the
// constraint pool.
// Default implementation use Cgl cut generators.
bool
DcoConGenerator::generateConstraints(BcpsConstraintPool &conPool)
{
  bool status = false;
#if defined(__OA__)
  OsiSolverInterface * solver = model_->solver();
#else
  OsiConicSolverInterface * solver = model_->solver();
#endif

#if defined(DISCO_DEBUG_MORE)
  // todo(aykut) no such function in DcoModel
  // std::cout << "model_->getNodeCount() = " << model_->getNodeCount()
  //	    << std::endl;
#endif

  //--------------------------------------------------
  // Start to generate constraints...
  //--------------------------------------------------

  assert(generator_ != NULL);

  int j;
  OsiCuts newOsiCuts;
  CglProbing* generator = dynamic_cast<CglProbing *>(generator_);

  if (generator) {
    // It is CglProbing - return tight column bounds
    CglTreeInfo info;
    generator->generateCutsAndModify(*solver, newOsiCuts, &info);
    const double * tightLower = generator->tightLower();
    const double * lower = solver->getColLower();
    const double * tightUpper = generator->tightUpper();
    const double * upper = solver->getColUpper();
    const double * solution = solver->getColSolution();

    int numberColumns = solver->getNumCols();
    double primalTolerance = 1.0e-8;
    for (j = 0; j < numberColumns; ++j) {
      if ( (tightUpper[j] == tightLower[j]) &&
	   (upper[j] > lower[j]) ) {
	// Fix column j
	solver->setColLower(j, tightLower[j]);
	solver->setColUpper(j, tightUpper[j]);
	if ( (tightLower[j] > solution[j] + primalTolerance) ||
	     (tightUpper[j] < solution[j] - primalTolerance) ) {
	  status = true;
	}
      }
    }
  }
  else {
    // Other Cgl cut generators
    generator_->generateCuts(*solver, newOsiCuts);
  }

  //--------------------------------------------------
  // Create blis constraints and remove zero length row cuts.
  //--------------------------------------------------

  int numNewConstraints = newOsiCuts.sizeRowCuts();
  for (j = 0; j < numNewConstraints; ++j) {
    OsiRowCut & rCut = newOsiCuts.rowCut(j);
    int len = rCut.row().getNumElements();

#ifdef DISCO_DEBUG_MORE
    std::cout << "Cut " << j<<": length = " << len << std::endl;
#endif
    if (len > 0) {
      // Create DcoConstraints from OsiCuts.
      DcoConstraint *blisCon = DcoOsiCutToConstraint(&rCut);
      conPool.addConstraint(blisCon);
    }
    else if (len == 0) {
      // Empty cuts
#ifdef DISCO_DEBUG
      std::cout << "WARNING: Empty cut from " << name_ << std::endl;
#endif
    }
    else {
#ifdef DISCO_DEBUG
      std::cout << "ERROR: Cut length = " << len << std::endl;
#endif
      // Error
      assert(0);
    }
  }

  // Adjust cut strategy.
  if ( (strategy_ == DcoCutStrategyAuto) &&
       (noConsCalls_ > DISCO_CUT_DISABLE) ) {
    strategy_ = DcoCutStrategyNone;
  }

  return status;
}

//#############################################################################
