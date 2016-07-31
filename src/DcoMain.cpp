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

#include <iostream>

#include "CoinError.hpp"
#include "CoinTime.hpp"


#if defined(__OA__)
  #include "OsiClpSolverInterface.hpp"
  typedef OsiClpSolverInterface SOLVER;
#else
  #include "OsiConicSolverInterface.hpp"
  #if defined(__OSI_MOSEK__)
    // use mosek as solver
    #include <OsiMosekSolverInterface.hpp>
    //#define IPM_SOLVER OsiMosekSolverInterface
    typedef OsiMosekSolverInterface SOLVER;
  #elif defined(__OSI_CPLEX__)
    // use cplex as solver
    #include <OsiCplexSolverInterface.hpp>
    typedef OsiCplexSolverInterface SOLVER;
  #elif defined(__COLA__)
    // use COLA as solver
    #include <ColaModel.hpp>
    typedef OsiCplexSolverInterface SOLVER;
  #else
    // use ipopt as solver
    #include <OsiIpoptSolverInterface.hpp>
    typedef OsiIpoptSolverInterface SOLVER;
  #endif
#endif

#include "CglFlowCover.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"

#include "DcoConfig.hpp"
#include "DcoModel.hpp"

#if  COIN_HAS_MPI
  #include "AlpsKnowledgeBrokerMPI.h"
#else
  #include "AlpsKnowledgeBrokerSerial.h"
#endif

// NOTE: gcc compiler doesn't recognize COIN_HAS_CLP, COIN_HAS_MPI

//#############################################################################

int main(int argc, char *argv[]) {
  try {
    //Set up lp solver
    // OsiConicSolverInterface * solver = new ColaModel();
    // solver.getModelPtr()->setDualBound(1.0e10);
    // solver.messageHandler()->setLogLevel(0);

    // if both mosek and cplex is available, choose mosek.
    // if both not available then choose cola
    OsiSolverInterface * solver = new SOLVER();
#if defined(__OA__)
    // set Clpc specific option to get unboundedness directions (if problem is unbounded)
    dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMoreSpecialOptions(0);
    // set clp specific options for no solver output
    dynamic_cast<OsiClpSolverInterface*>(solver)->setHintParam(OsiDoReducePrint,false,OsiHintDo, 0);
    //lpSolver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
#endif
    solver->setHintParam(OsiDoReducePrint,true,OsiHintDo, 0);
    // Create DisCO model
    DcoModel model;
    model.setSolver(solver);

#ifdef  COIN_HAS_MPI
    AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
    AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif
    // Search for best solution
    broker.search(&model);
    // Report the best solution found and its ojective value
    broker.printBestSolution();
    delete solver;
  }
  catch(CoinError& er) {
    std::cerr << "\nDISCO ERROR: \"" << er.message()
	      << "\""<< std::endl
	      << "             from function \"" << er.methodName()
	      << "\""<< std::endl
	      << "             from class \"" << er.className()
	      << "\"" << std::endl;
  }
  catch(...) {
    std::cerr << "Something went wrong!" << std::endl;
  }
  return 0;
}
//#############################################################################
