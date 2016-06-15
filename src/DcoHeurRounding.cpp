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

#include <CoinMessageHandler.hpp>
//#include <CoinMessage.hpp>

#include "DcoHeurRounding.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"

DcoHeurRounding::DcoHeurRounding(DcoModel * model, char const * name,
                                 DcoHeurStrategy strategy, int frequency)
  : DcoHeuristic(model, name, strategy, frequency) {
  setType(DcoHeurTypeRounding);
}

DcoSolution * DcoHeurRounding::searchSolution() {
  /// todo(aykut) disable heuristic search for now.
  return NULL:
  if (strategy() == DcoHeurStrategyNone) {
    // This heuristic has been disabled.
    return NULL;
  }
  DcoModel * dcom = model();
  // get required pointers for log messages
  CoinMessageHandler * message_handler = dcom->dcoMessageHandler_;
  CoinMessages * messages = dcom->dcoMessages_;

  // call bound_fix function to compute find up_fix and down_fix
  int num_rel_cols = dcom->numRelaxedCols();
  // if down_fix[i] is 0 variable i can be fixed to its lower bound.
  int * down_fix = new int[num_rel_cols]();
  int * up_fix = new int[num_rel_cols]();
  // == compute up and down fixes
  bound_fix(down_fix, up_fix);

  // check up_fix and down_fix
  return NULL;
}

void DcoHeurRounding::bound_fix(int * down_fix, int * up_fix) {
  DcoModel * dcom = model();
  // get required pointers for log messages
  CoinMessageHandler * message_handler = dcom->dcoMessageHandler_;
  CoinMessages * messages = dcom->dcoMessages_;

  int num_rel_cols = dcom->numRelaxedCols();
  int num_rows = dcom->solver()->getNumRows();
  char const * row_sense = dcom->solver()->getRowSense();

  double infinity = dcom->solver()->getInfinity();
  // iterate over rows and update up and down fixed.
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]=='R') {
      if (dcom->solver()->getColUpper()[i]>=infinity
          and dcom->solver()->getColLower()[i]<=-infinity) {
        // both upper and lower bound are not finite,
        // do nothing
        continue;
      }
      // upper bound is infinity
      if (dcom->solver()->getColUpper()[i]>=infinity) {
        bound_fix2('G', i, down_fix, up_fix);
        continue;
      }
      // lower bound is negative infinity
      if (dcom->solver()->getColLower()[i]<=-infinity) {
        bound_fix2('L', i, down_fix, up_fix);
        continue;
      }
      // bounds are not infinity if we are here
      // note that once we have R rows with finite bounds
      // bound fixing procedure is same as E rows.
      bound_fix2('E', i, down_fix, up_fix);
    }
    else if (row_sense[i]=='N') {
      // row is not restricted, do nothing
      continue;
    }
    else if (row_sense[i]=='E' or row_sense[i]=='L' or row_sense[i]=='G') {
      bound_fix2(row_sense[i], i, down_fix, up_fix);
    }
    else {
      // unknown row sense
      std::stringstream error;
      error << "Unknown row sense "
            << row_sense[i];
      message_handler->message(9998, "Dco", error.str().c_str(),
                               'E', 0)
        << CoinMessageEol;
    }
  }
}

void DcoHeurRounding::bound_fix2(char sense, int row_index, int * down_fix, int * up_fix) {
  //char row_sense = dcom->solver()->getRowSense()[row_index];

  CoinPackedMatrix const * matrix = model()->solver()->getMatrixByRow();
  int const * indices = matrix->getIndices();
  double const * values = matrix->getElements();
  int const * lengths = matrix->getVectorLengths();
  int const * starts = matrix->getVectorStarts();
  double tailoff = model()->dcoPar()->entry(DcoParams::tailOff);
  // iterate over varaibles of row row_index
  for (int i=starts[row_index]; i<starts[row_index]+lengths[row_index]; ++i) {
    if (-tailoff<=values[i] and values[i]<=tailoff) {
      // coeffciient is very close to 0. log warning message
      std::stringstream warning;
      warning << "Coefficient of variable "
              << indices[i]
              << " in row "
              << row_index
              << " is "
              << values[i]
              << ", very close to 0.";
      model()->dcoMessageHandler_->message(3000, "Dco", warning.str().c_str(),
                                           'W', 0)
        << CoinMessageEol;
    }
    if (sense=='E') {
      up_fix[indices[i]]++;
      down_fix[indices[i]]++;
    }
    else if (sense=='L') {
      if (values[i]>0) {
        up_fix[indices[i]]++;
      }
      else {
        down_fix[indices[i]]++;
      }
    }
    else if (sense=='G') {
      if (values[i]>0) {
        down_fix[indices[i]]++;
      }
      else {
        up_fix[indices[i]]++;
      }
    }
    else {
      // unknown row sense
      std::stringstream error;
      error << "Unexpected row sense "
            << sense;
      model()->dcoMessageHandler_->message(9998, "Dco", error.str().c_str(),
                                           'E', 0)
        << CoinMessageEol;
    }
  }
}
