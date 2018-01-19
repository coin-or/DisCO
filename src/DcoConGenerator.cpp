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


#include "DcoConGenerator.hpp"

void DcoConGeneratorStats::reset() {
  numConsGenerated_ = 0;
  numConsUsed_ = 0;
  time_ = 0.0;
  numCalls_ = 0;
  numNoConsCalls_ = 0;
}

DcoConGenerator::DcoConGenerator(DcoModel * model,
                                 DcoConstraintType type,
                                 char const * name,
                                 DcoCutStrategy strategy,
                                 int frequency):
  name_(name), type_(type), model_(model), strategy_(strategy),
  frequency_(frequency) {
  stats_.reset();
}

DcoConGenerator::~DcoConGenerator() {
  model_ = NULL;
}
