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


#include "DcoModel.hpp"
#include "DcoSolution.hpp"
#include "DcoTreeNode.hpp"
#include "DcoBranchObject.hpp"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

#if defined(__OA__)
  #include <OsiClpSolverInterface.hpp>
  typedef OsiClpSolverInterface LINEAR_SOLVER;
// get SOCO solver
#else
  #include <OsiConicSolverInterface.hpp>
#if defined(__OSI_MOSEK__)
  // use mosek
  #include <OsiMosekSolverInterface.hpp>
  typedef OsiMosekSolverInterface SOCO_SOLVER;
#elif defined(__OSI_CPLEX__)
  // use cplex
  #include <OsiCplexSolverInterface.hpp>
  typedef OsiCplexSolverInterface SOCO_SOLVER;
#elif defined(__OSI_IPOPT__)
  // use ipopt
  #include <OsiIpoptSolverInterface.hpp>
  typedef OsiIpoptSolverInterface SOCO_SOLVER;
#elif defined(__COLA__)
  // use cola
  #include <ColaModel.hpp>
  typedef ColaModel SOCO_SOLVER;
#endif
#endif

int main(int argc, char *argv[]) {

  // print host info and sleep
  // this is for debugging
  // {
  //   int i = 0;
  //   char hostname[256];
  //   gethostname(hostname, sizeof(hostname));
  //   printf("PID %d on %s ready for attach\n", getpid(), hostname);
  //   fflush(stdout);
  //   while (0 == i)
  //     sleep(5);
  // }

#if defined(__OA__)
  OsiSolverInterface * solver = new LINEAR_SOLVER();
  solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
  // clp specific options for getting unboundedness directions
  dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMoreSpecialOptions(0);
  dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setLogLevel(0);
#else
  OsiConicSolverInterface * solver = new SOCO_SOLVER();
  solver->setHintParam(OsiDoReducePrint, true, OsiHintTry);
#endif
  // Create DisCO model
  DcoModel model;
  model.setSolver(solver);
#ifdef  COIN_HAS_MPI
  AlpsKnowledgeBrokerMPI broker(argc, argv, model);
  //broker.passInMessageHandler(model.dcoMessageHandler_);
  // Register model, solution, and tree node
  broker.registerClass(AlpsKnowledgeTypeModel, new DcoModel);
  broker.registerClass(AlpsKnowledgeTypeSolution, new DcoSolution);
  broker.registerClass(AlpsKnowledgeTypeNode, new DcoTreeNode);
  broker.registerClass(AlpsKnowledgeTypeNodeDesc, new DcoNodeDesc);
  broker.registerClass(999, new DcoBranchObject(-1, 0.0, 0.0));
#else
  AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif

  // Search for best solution
  broker.search(&model);
  model.reportFeasibility();
  // Report the best solution found and its ojective value
  broker.printBestSolution();

  delete solver;
  return 0;
}
//#############################################################################
