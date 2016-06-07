#include <OsiClpSolverInterface.hpp>

#include "DcoModel.hpp"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

int main(int argc, char *argv[]) {
#if defined(__OA__)
  // todo(aykut) what about setting solver log level
  OsiSolverInterface * solver = new OsiClpSolverInterface();
  // for unboundedness directions set option
  dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMoreSpecialOptions(0);
  solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
  dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setLogLevel(0);
#else
#if defined(__OSI_MOSEK__)
  OsiConicSolverInterface * solver = new OsiMosekSolverInterface();
#else
#if defined(__OSI_CPLEX__)
  OsiConicSolverInterface * solver = new OsiCplexSolverInterface();
#else
#if defined(__COLA__)
  OsiConicSolverInterface * solver = new ColaModel();
#endif
#endif
#endif
#endif
  // Create DisCO model
  DcoModel model;
  model.setSolver(solver);
#ifdef  COIN_HAS_MPI
  AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
  AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif
  // Search for best solution
  //todo(aykut) log levels are set at setupSelf. Maybe they should be set
  //earlier.
  // since broker starts printing logs before setupSelf is even called.
  broker.search(&model);
  // Report the best solution found and its ojective value
  broker.printBestSolution();
  model.reportFeasibility();
  delete solver;
  return 0;
}
//#############################################################################
