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

#ifndef DcoPresolve_H_
#define DcoPresolve_H_

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiPresolve.hpp"
#include <OsiConicSolverInterface.hpp>

//#############################################################################

/** A interface to Osi/Coin Presolve. */
class DcoPresolve : public OsiPresolve
{
private:

    CoinPresolveMatrix *preMatrix_;
    CoinPostsolveMatrix *postMatrix_;

public:

    /** Default constructor (empty object) */
    DcoPresolve() :
        preMatrix_(0),
        postMatrix_(0) {}

    /** Virtual destructor */
    virtual ~DcoPresolve() {
        delete preMatrix_;
        delete postMatrix_;
    }

    /** Presolve */
    virtual OsiConicSolverInterface *preprocess(OsiConicSolverInterface & origModel,
                                           double feasibilityTolerance=0.0,
                                           bool keepIntegers=true,
                                           int numberPasses=5,
                                           const char * prohibited=NULL);

    /** Postsolve */
    virtual void postprocess(bool updateStatus=true);
};

#endif

//#############################################################################
