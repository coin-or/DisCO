#include <OsiClpSolverInterface.hpp>

#include "DcoModel.hpp"
#include "DcoSolution.hpp"
#include "DcoTreeNode.hpp"
#include "DcoBranchObject.hpp"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
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
  OsiSolverInterface * solver = new OsiClpSolverInterface();
  // for unboundedness directions set option
  dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setMoreSpecialOptions(0);
  solver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
  dynamic_cast<OsiClpSolverInterface*>(solver)->getModelPtr()->setLogLevel(0);
#elif defined(__OSI_MOSEK__)
  OsiConicSolverInterface * solver = new OsiMosekSolverInterface();
#elif defined(__OSI_CPLEX__)
  OsiConicSolverInterface * solver = new OsiCplexSolverInterface();
#elif defined(__COLA__)
  OsiConicSolverInterface * solver = new ColaModel();
#endif

  // Create DisCO model
  DcoModel model;
  model.setSolver(solver);
#ifdef  COIN_HAS_MPI
  AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
  AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif
  //broker.passInMessageHandler(model.dcoMessageHandler_);
  // Register model, solution, and tree node
  broker.registerClass(AlpsKnowledgeTypeModel, new DcoModel);
  broker.registerClass(AlpsKnowledgeTypeSolution, new DcoSolution);
  broker.registerClass(AlpsKnowledgeTypeNode, new DcoTreeNode);
  broker.registerClass(AlpsKnowledgeTypeNodeDesc, new DcoNodeDesc);
  broker.registerClass(999, new DcoBranchObject(-1, 0.0, 0.0));

  // Search for best solution
  broker.search(&model);
  // Report the best solution found and its ojective value
  broker.printBestSolution();

  model.reportFeasibility();
  delete solver;
  return 0;
}
//#############################################################################
