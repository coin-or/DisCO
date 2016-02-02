#ifndef DcoModel_hpp_
#define DcoModel_hpp_

#include <BcpsModel.h>
#include <OsiConicSolverInterface.hpp>
#include <OsiLorentzCone.hpp>
#include <BcpsBranchStrategy.h>


#include "DcoParams.hpp"

class DcoModel: public BcpsModel {
  //------------------------------------------------------
  // LP SOLVER.
  //------------------------------------------------------
#if defined(__OA__)
  OsiSolverInterface * solver_;
#else
  OsiConicSolverInterface * solver_;
#endif
  //------------------------------------------------------
  // PROBLEM DATA
  //------------------------------------------------------
  /** Problem matrix */
  CoinPackedMatrix * matrix_;
  /** Variable and constraint bounds. */
  //@{
  double *colLB_;
  double *colUB_;
  double *rowLB_;
  double *rowUB_;
  //@}
  /** Number of columns/rows/elements */
  //@{
  int numCols_;
  int numRows_;
  //@}
  /** Objective function. */
  //@{
  double objSense_;
  double *objCoef_;
  //@}
  /** Cone data. */
  //@{
  OsiLorentzCone * cone_;
  int numCones_;
  //@}
  /** Column types. */
  //@{
  int numIntObjects_;
  int * intColIndices_;  // size of numIntObjects_
  //@}
  //------------------------------------------------------
  // SOLUTION.
  //------------------------------------------------------
  int numSolutions_;
  /** Incumbent objective value. */
  double incObjValue_;
  /** Incumbent */
  double * incumbentSol_;
  /** Cutoff in lp solver. */
  double cutoff_;
  /** Cutoff increment. */
  double cutoffInc_;
  //------------------------------------------------------
  // SEARCHING.
  //------------------------------------------------------
  int * intObjIndices_; // size of numCols_
  char * colType_;
  /** Starting var/con bounds for processing each node */
  //@{
  double * startColLB_;
  double * startColUB_;
  double * startRowLB_;
  double * startRowUB_;
  //@}
  /** Variable selection function. */
  BcpsBranchStrategy * branchStrategy_;
  BcpsBranchStrategy * rampUpBranchStrategy_;
  /** Number of objects. */
  int numObjects_;
  /** The set of objects. */
  BcpsObject ** objects_;
  /** The objects that can be shared. */
  char * sharedObjectMark_;
  /** Active node. */
  AlpsTreeNode * activeNode_;
  //------------------------------------------------------
  // PARAMETERS, STATISTICS, and MESSAGE
  //------------------------------------------------------
  /** Dco parameters. */
  DcoParams * DcoPar_;
  /** Message handler. */
  CoinMessageHandler * dcoMessageHandler_;
  /** Dco messages. */
  CoinMessages dcoMessages_;
  /** Number of processed nodes. */
  int numNodes_;
  /** Number of lp(Simplex) iterations. */
  int numIterations_;
  /** Average number of lp iterations to solve a subproblem. */
  int aveIterations_;
  /** The number of passes during bounding procedure.*/
  int boundingPass_;
  /** Current relative optimal gap. */
  double currRelGap_;
  /** Current absolute optimal gap. */
  double currAbsGap_;
public:
  DcoModel();
  virtual ~DcoModel();
#if defined(__OA__)
  void setSolver(OsiSolverInterface * solver);
#else
  void setSolver(OsiConicSolverInterface * solver);
#endif
};

#endif
