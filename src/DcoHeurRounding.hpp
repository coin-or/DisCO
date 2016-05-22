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
// This file is modified from COIN/Cbc/CbcHeuristic.hpp
//#############################################################################


#ifndef DcoHeurRounding_hpp_
#define DcoHeurRounding_hpp_

// STL headers

// Disco headers
#include "DcoHeuristic.hpp"

/*!
  Implements simple rounding heuristic described in Achterberg's dissretation.
*/

/** Heuristic base class */
class DcoHeurRounding {
public:
  ///@name Constructors and Destructor.
  //@{
  /// Useful constructor.
  DcoHeurRounding(DcoModel * model, char const * name,
                  DcoHeurStrategy strategy, int frequency);
  /// Destructor.
  virtual ~DcoHeurRounding() { }
  //@}

  ///@name Finding solutions.
  //@{
  /// returns a solution if found, NULL otherwise.
  virtual DcoSolution * searchSolution();
  //@}

private:
  /// Disable default constructor.
  DcoHeurRounding();
  /// Disable copy constructor.
  DcoHeurRounding(const DcoHeurRounding & other);
  /// Disable copy assignment operator
  DcoHeurRounding & operator=(const DcoHeurRounding & rhs);
};

#endif
