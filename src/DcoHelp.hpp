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

#ifndef DcoHelp_h_
#define DcoHelp_h_

#include "AlpsEncoded.h"

#include "Dco.hpp"

class CoinWarmStartBasis;
class OsiRowCut;
class DcoConstraint;
class DcoModel;

//#############################################################################

/** Convert a OsiRowCut to a Dco Contraint. */
DcoConstraint * DcoOsiCutToConstraint(const OsiRowCut *rowCut);

/** Strong branching on a variable colInd. */
DcoReturnStatus DcoStrongBranch(DcoModel *model, double objValue, int colInd, double x,
                                  const double *saveLower, const double *saveUpper,
		                          bool &downKeep, bool &downFinished, double &downDeg,
		                          bool &upKeep, bool &upFinished, double &upDeg);

/** Pack coin warm start into an encoded object. */
int DcoEncodeWarmStart(AlpsEncoded *encoded, const CoinWarmStartBasis *ws);

/** Unpack coin warm start from an encoded object. */
CoinWarmStartBasis *DcoDecodeWarmStart(AlpsEncoded &encoded,
					AlpsReturnStatus *rc);

/** Compute and return a hash value of an Osi row cut. */
double DcoHashingOsiRowCut(const OsiRowCut *rowCut, 
			    const DcoModel *model);

/** Check if a row cut parallel with another row cut. */
bool DcoParallelCutCut(OsiRowCut * rowCut1,
			OsiRowCut * rowCut2,
			double threshold = 1.0);

/** Check if a row cut parallel with a constraint. */
bool DcoParallelCutCon(OsiRowCut * rowCut,
			DcoConstraint * con,
			double threshold = 1.0);

/** Check if a row cut parallel with a constraint. */
bool DcoParallelConCon(DcoConstraint * con1,
			DcoConstraint * con2,
			double threshold = 1.0);


#endif
