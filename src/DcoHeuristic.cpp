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


#include "DcoHeuristic.hpp"

void DcoHeurStats::reset() {
  numCalls_ = 0;
  numNoSolCalls_ = 0;
  time_ = 0.0;
  numSolutions_ = 0;
}

DcoHeuristic::DcoHeuristic(DcoModel * model, char const * name,
                           DcoHeurStrategy strategy,
                           int frequency) {
  model_ = model;
  name_ = name;
  strategy_ = strategy;
  frequency_ = frequency;
  type_ = DcoHeurTypeNotSet;
  stats_.reset();
}
