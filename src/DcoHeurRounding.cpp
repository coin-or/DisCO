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


#include <numeric>

#include <CoinMessageHandler.hpp>
//#include <CoinMessage.hpp>

#include "DcoHeurRounding.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoModel.hpp"
#include "DcoSolution.hpp"
#include "DcoConicConstraint.hpp"


DcoHeurRounding::DcoHeurRounding(DcoModel * model, char const * name,
                                 DcoHeurStrategy strategy, int frequency)
  : DcoHeuristic(model, name, strategy, frequency) {
  setType(DcoHeurTypeRounding);
}

DcoSolution * DcoHeurRounding::searchSolution() {
  if (strategy() == DcoHeurStrategyNone) {
    // This heuristic has been disabled.
    return NULL;
  }
  else if (strategy() == DcoHeurStrategyAuto) {
    // todo(aykut) the disable threshold should be a parameter in DcoParams.
    int disable_threshold = 1000000;
    if (stats().numNoSolCalls() > disable_threshold) {
      return NULL;
    }
  }
  DcoSolution * dco_sol = NULL;
  // get pointers for message logging
  CoinMessageHandler * message_handler = model()->dcoMessageHandler_;
  CoinMessages * messages = model()->dcoMessages_;

  // Get a copy of original matrix (and by row for rounding);
  CoinPackedMatrix const * matrix = model()->solver()->getMatrixByCol();
  CoinPackedMatrix const * matrixByRow = model()->solver()->getMatrixByRow();

#if defined(__OA__)
  OsiSolverInterface * solver = model()->solver();
#else
  OsiConicSolverInterface * solver = model()->solver();
#endif
  double const * colLower = solver->getColLower();
  double const * colUpper = solver->getColUpper();
  double const * rowLower = solver->getRowLower();
  double const * rowUpper = solver->getRowUpper();
  double const * obj = solver->getObjCoefficients();
  double integerTol = model()->dcoPar()->entry(DcoParams::integerTol);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  int numRows = matrix->getNumRows();
  int numIntegers = model()->numIntegerCols();
  int const * integerCols = model()->integerCols();
  double direction = solver->getObjSense();
  double sol_quality = direction * solver->getObjValue();
  //double newSolutionValue = direction * solver->getObjValue();

  // Column copy
  double const * element = matrix->getElements();
  int const * row = matrix->getIndices();
  int const * columnStart = matrix->getVectorStarts();
  int const * columnLength = matrix->getVectorLengths();
  // Row copy
  double const * elementByRow = matrixByRow->getElements();
  int const * column = matrixByRow->getIndices();
  int const * rowStart = matrixByRow->getVectorStarts();
  int const * rowLength = matrixByRow->getVectorLengths();

  // Get solution array for heuristic solution
  int numCols = solver->getNumCols();
  double * sol = new double [numCols];
  std::copy(solver->getColSolution(), solver->getColSolution()+numCols, sol);

  double * rowActivity = new double[numRows]();

  for (int i=0; i<numCols; i++) {
    double value = sol[i];
    if (value) {
      for (int j=columnStart[i]; j<columnStart[i]+columnLength[i]; j++) {
        int iRow = row[j];
        rowActivity[iRow] += value*element[j];
      }
    }
  }
  // check was feasible - if not adjust (cleaning may move)
  for (int i=0; i< numRows; i++) {
    if(rowActivity[i] < rowLower[i]) {
      rowActivity[i] = rowLower[i];
    }
    else if(rowActivity[i] > rowUpper[i]) {
      rowActivity[i] = rowUpper[i];
    }
  }
  for (int i=0; i<numIntegers; i++) {
    int iColumn = integerCols[i];
    double value = sol[iColumn];
    if (fabs(floor(value + 0.5) - value) > integerTol) {
      double below = floor(value);
      double newValue = sol[iColumn];
      double cost = direction * obj[iColumn];
      double move;
      if (cost > 0.0) {
        // try up
        move = 1.0 - (value - below);
      }
      else if (cost < 0.0) {
        // try down
        move = below - value;
      }
      else {
        // won't be able to move unless we can grab another variable
        // just for now go down
        move = below-value;
      }
      newValue += move;
      sol[iColumn] = newValue;
      sol_quality += move * cost;
      for (int j=columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn]; j++) {
        int iRow = row[j];
        rowActivity[iRow] += move * element[j];
      }
    }
  }

  double penalty = 0.0;
  // see if feasible
  for (int i=0; i< numRows; i++) {
    double value = rowActivity[i];
    double thisInfeasibility = 0.0;
    if (value < rowLower[i] - primalTolerance) {
      thisInfeasibility = value - rowLower[i];
    }
    else if (value > rowUpper[i] + primalTolerance) {
      thisInfeasibility = value - rowUpper[i];
    }
    if (thisInfeasibility) {
      // See if there are any slacks I can use to fix up
      // maybe put in coding for multiple slacks?
      double bestCost = 1.0e50;
      int iBest = -1;
      double addCost = 0.0;
      double newValue = 0.0;
      double changeRowActivity = 0.0;
      double absInfeasibility = fabs(thisInfeasibility);
      for (int k=rowStart[i]; k<rowStart[i]+rowLength[i]; k++) {
        int iColumn = column[k];
        if (columnLength[iColumn] == 1) {
          double currentValue = sol[iColumn];
          double elementValue = elementByRow[k];
          double lowerValue = colLower[iColumn];
          double upperValue = colUpper[iColumn];
          double gap = rowUpper[i] - rowLower[i];
          double absElement = fabs(elementValue);
          if (thisInfeasibility * elementValue > 0.0) {
            // we want to reduce
            if ((currentValue - lowerValue) * absElement >=
                absInfeasibility) {

              // possible - check if integer
              double distance = absInfeasibility / absElement;
              double thisCost =
                -direction * obj[iColumn] * distance;
              if (solver->isInteger(iColumn)) {
                distance = ceil(distance - primalTolerance);
                assert (currentValue - distance >=
                        lowerValue - primalTolerance);
                if (absInfeasibility - distance * absElement
                    < -gap - primalTolerance)
                  thisCost = 1.0e100; // no good
                else
                  thisCost =
                    -direction*obj[iColumn]*distance;
              }
              if (thisCost < bestCost) {
                bestCost = thisCost;
                iBest = iColumn;
                addCost = thisCost;
                newValue = currentValue - distance;
                changeRowActivity = -distance * elementValue;
              }
            }
          }
          else {
            // we want to increase
            if ((upperValue - currentValue) * absElement >=
                absInfeasibility) {
              // possible - check if integer
              double distance = absInfeasibility / absElement;
              double thisCost =
                direction * obj[iColumn] * distance;
              if (solver->isInteger(iColumn)) {
                distance = ceil(distance - 1.0e-7);
                assert (currentValue - distance <=
                        upperValue + primalTolerance);
                if (absInfeasibility - distance * absElement
                    < -gap - primalTolerance)
                  thisCost = 1.0e100; // no good
                else
                  thisCost =
                    direction*obj[iColumn]*distance;
              }
              if (thisCost < bestCost) {
                bestCost = thisCost;
                iBest = iColumn;
                addCost = thisCost;
                newValue = currentValue + distance;
                changeRowActivity = distance * elementValue;
              }
            }
          }
        }
      }
      if (iBest >= 0) {
        /*printf("Infeasibility of %g on row %d cost %g\n",
          thisInfeasibility,i,addCost);*/
        sol[iBest] = newValue;
        thisInfeasibility = 0.0;
        sol_quality += addCost;
        rowActivity[i] += changeRowActivity;
      }
      penalty += fabs(thisInfeasibility);
    }
  }

  // Could also set SOS (using random) and repeat
  if (!penalty) {
    // Got a feasible solution. Try to improve.
    //seed_++;
    //CoinSeedRandom(seed_);
    // Random number between 0 and 1.

    double randomNumber = CoinDrand48();
    int start[2];
    int end[2];
    int iRandom = (int) (randomNumber * ((double) numIntegers));
    start[0] = iRandom;
    end[0] = numIntegers;
    start[1] = 0;
    end[1] = iRandom;
    // todo(aykut) number of passes is a parameter, use it. set default value
    // to 2.
    for (int iPass = 0; iPass<2; iPass++) {
      for (int i=start[iPass]; i<end[iPass]; i++) {
        int iColumn = integerCols[i];
#ifdef DISCO_DEBUG
        double value = sol[iColumn];
        assert(fabs(floor(value + 0.5) - value) < integerTol);
#endif
        double cost = direction * obj[iColumn];
        double move = 0.0;
        if (cost > 0.0)
          move = -1.0;
        else if (cost < 0.0)
          move = 1.0;
        while (move) {
          bool good = true;
          double newValue = sol[iColumn] + move;
          if (newValue < colLower[iColumn] - primalTolerance||
              newValue > colUpper[iColumn] + primalTolerance) {
            move = 0.0;
          }
          else {
            // see if we can move
            for (int j=columnStart[iColumn];
                 j<columnStart[iColumn]+columnLength[iColumn]; j++) {
              int iRow = row[j];
              double newActivity =
                rowActivity[iRow] + move*element[j];
              if (newActivity < rowLower[iRow] - primalTolerance
                  ||
                  newActivity > rowUpper[iRow]+primalTolerance) {
                good = false;
                break;
              }
            }
            if (good) {
              sol[iColumn] = newValue;
              sol_quality += move * cost;
              for (int j=columnStart[iColumn];
                   j<columnStart[iColumn]+columnLength[iColumn]; j++) {
                int iRow = row[j];
                rowActivity[iRow] += move*element[j];
              }
            }
            else {
              move=0.0;
            }
          }
        }
      }
    }
    // paranoid check
    std::fill_n(rowActivity, numRows, 0.0);
    for (int i=0; i<numCols; i++) {
      double value = sol[i];
      if (value) {
        for (int j=columnStart[i]; j<columnStart[i]+columnLength[i]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += value * element[j];
        }
      }
    }


    // check was approximately feasible
    bool feasible = true;
    for (int i=0; i<numRows; i++) {
      if (rowActivity[i] < rowLower[i]) {
        if (rowActivity[i] < rowLower[i] - 1000.0*primalTolerance)
          feasible = false;
      }
      else if (rowActivity[i] > rowUpper[i]) {
        if (rowActivity[i] > rowUpper[i] + 1000.0*primalTolerance)
          feasible = false;
      }
    }
    if (feasible!=false) {
      // check whether the solution is conic feasible
      double cone_tol = model()->dcoPar()->entry(DcoParams::coneTol);
      //cone_tol = 1e-10;
      int num_linear_rows = model()->getNumCoreLinearConstraints();
      int num_conic_rows = model()->getNumCoreConicConstraints();
      for (int i=num_linear_rows; i<num_linear_rows+num_conic_rows; ++i) {
        DcoConicConstraint * con =
          dynamic_cast<DcoConicConstraint*> (model()->getConstraints()[i]);
        int const * members = con->coneMembers();
        DcoLorentzConeType type = con->coneType();
        int size = con->coneSize();
        double * values = new double[size];
        for (int i=0; i<size; ++i) {
          values[i] = sol[members[i]];
        }
        double term1 = 0.0;
        double term2 = 0.0;
        if (type==DcoLorentzCone) {
          term1 = values[0];
          term2 = std::inner_product(values+1, values+size, values+1, 0.0);
          term2 = sqrt(term2);
        }
        else if (type==DcoRotatedLorentzCone) {
          term1 = 2.0*sol[members[0]]*sol[members[1]];
          term2 = std::inner_product(values+2, values+size, values+2, 0.0);
        }
        else {
          message_handler->message(DISCO_UNKNOWN_CONETYPE,
                                   *messages)
            << type << CoinMessageEol;
        }
        if (term1-term2<-cone_tol) {
          delete[] values;
          feasible = false;
          break;
        }
        delete[] values;
      }
    }

    if (feasible) {
      // new solution found, store solution.
      dco_sol = new DcoSolution(numCols, sol, sol_quality);
      dco_sol->setBroker(model()->broker_);
      stats().addNumSolutions();
    }
    else {
      // update statistics
      stats().addNoSolCalls();
      // Can easily happen
      //printf("Debug DcoHeurRound giving bad solution\n");
    }
  }
  delete [] sol;
  delete [] rowActivity;
  return dco_sol;
}


DcoSolution * DcoHeurRounding::searchSolution2() {
  /// todo(aykut) disable heuristic search for now.
  return NULL;
  if (strategy() == DcoHeurStrategyNone) {
    // This heuristic has been disabled.
    return NULL;
  }
  DcoModel * dcom = model();
  // get required pointers for log messages
  //CoinMessageHandler * message_handler = dcom->dcoMessageHandler_;
  //CoinMessages * messages = dcom->dcoMessages_;

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
  //CoinMessages * messages = dcom->dcoMessages_;

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
