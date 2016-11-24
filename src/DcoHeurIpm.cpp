#include <numeric>

#include <CoinMessageHandler.hpp>
//#include <CoinMessage.hpp>

#include "DcoHeurIpm.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoModel.hpp"
#include "DcoSolution.hpp"
#include "DcoConicConstraint.hpp"

#if defined(__MOSEK_EXIST__)
  // use mosek
  #include <OsiMosekSolverInterface.hpp>
  typedef OsiMosekSolverInterface SOCO_SOLVER;
#elif defined(__CPLEX_EXIST__)
  // use cplex
  #include <OsiCplexSolverInterface.hpp>
  typedef OsiCplexSolverInterface SOCO_SOLVER;
#endif

DcoHeurIpm::DcoHeurIpm(DcoModel * model, char const * name,
                       DcoHeurStrategy strategy, int frequency)
  : DcoHeuristic(model, name, strategy, frequency) {
  setType(DcoHeurTypeIpm);
}

// this heur should only be called in case of __OA__
// For this heuristic to work, we need an IPM solver
// do we assume all integer variables are feasible?
DcoSolution * DcoHeurIpm::searchSolution() {
#if defined(__OA__)
  DcoSolution * dco_sol = NULL;
  // get pointers for message logging
  CoinMessageHandler * message_handler = model()->dcoMessageHandler_;
  CoinMessages * messages = model()->dcoMessages_;
  // update statistics
  // check integrality, if integer infeasible then return null
  int numColsInf;
  int numRowsInf;
  double colInf;
  double rowInf;
  DcoSolution * sol = model()->feasibleSolution(numColsInf, colInf,
                                              numRowsInf, rowInf);
  // Debug message
  // IPM heuristic is called, solution is not integral
  // IPM heuristic is called and solution with quality ... is found.
  if (numColsInf==0) {
#if defined(__CPLEX_EXIST__) || (__MOSEK_EXIST__)
    // load problem to solver, model()->solver() is a linear solver
    // since __OA__
    stats().addCalls();
    std::cout << "Integer conic infeasible solution!" << std::endl;
    OsiConicSolverInterface * solver = new SOCO_SOLVER();
    OsiSolverInterface * si = model()->solver();
    solver->loadProblem(*(si->getMatrixByRow()),
                        si->getColLower(),
                        si->getColUpper(),
                        si->getObjCoefficients(),
                        si->getRowLower(),
                        si->getRowUpper());
    // add cones
    int const * coneStart = model()->coneStart();
    int const * coneMembers = model()->coneMembers();
    for (int i=0; i<model()->getNumCoreConicConstraints(); ++i) {
      // do not relax conic constraints, add them to the conic solver
      OsiLorentzConeType osi_type;
      if (model()->coneType()[i]==1) {
        osi_type = OSI_QUAD;
      }
      else if (model()->coneType()[i]==2) {
        osi_type = OSI_RQUAD;
      }
      // for (int j=0; j<(coneStart[i+1]-coneStart[i]); ++j) {
      //   std::cout << coneMembers[coneStart[i]+j] << " ";
      // }
      // std::cout << std::endl;
      solver->addConicConstraint(osi_type, coneStart[i+1]-coneStart[i],
                                  coneMembers+coneStart[i]);
    }
    // fix integer solutions
    int const * integer_cols = model()->integerCols();
    double const * sol = model()->solver()->getColSolution();
    for (int i=0; i<model()->numIntegerCols(); ++i) {
      solver->setColBounds(integer_cols[i], sol[integer_cols[i]],
                           sol[integer_cols[i]]);
    }
    // solve
    solver->initialSolve();
    int num_cols = solver->getNumCols();
    if (solver->isProvenOptimal()) {
      std::cout << "Feasible solution found!"
                << si->getObjValue()
                << " "
                << solver->getObjValue()
                << std::endl;
      dco_sol = new DcoSolution(num_cols, solver->getColSolution(),
                                solver->getObjValue());
      stats().addNumSolutions();
    }
    else {
      stats().addNoSolCalls();
    }
#endif
  }
#endif
  return dco_sol;
}
