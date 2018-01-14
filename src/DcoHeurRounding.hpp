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
 * Ted Ralphs. All Rights Reserved.                                          *
 *===========================================================================*/


#ifndef DcoHeurRounding_hpp_
#define DcoHeurRounding_hpp_

// Disco headers
#include "DcoHeuristic.hpp"

/*!
  Implements simple rounding heuristic described in Achterberg's dissretation.

  # Ideas:
  When rounding solutions we can round integer leading variables up.

*/

/** Heuristic base class */
class DcoHeurRounding: virtual public DcoHeuristic {
  void bound_fix(int * down_fix, int * up_fix);
  void bound_fix2(char sense, int row_index, int * down_fix, int * up_fix);
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
  virtual DcoSolution * searchSolution2();
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
