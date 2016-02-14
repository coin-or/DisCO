#ifndef DcoModel_hpp_
#define DcoModel_hpp_

#include <BcpsModel.h>
#include <OsiConicSolverInterface.hpp>
#include <OsiLorentzCone.hpp>
#include <BcpsBranchStrategy.h>


#include "DcoParams.hpp"
#include "DcoConstraint.hpp"


/*
  todo(aykut): Bcps ideas
  We keep relaxed objects (the whole object or integrality of the object) in
  list BcpsObject ** relaxed_.

  Can subproblems have different relaxed objects? Yes they can. But make sure
  all relaxed objects are in the relaxed_.

  Subproblems might have different solvers? A subprobllem might be an LP, SOCO,
  MILP or even MISOCO. How will DcoTreeNode::bound know about the solver to be
  used?

  When subproblem is solved with the solver we go back and check the
  feasibility of the objects that are relaxed in the subproblem (objects in
  relaxed_). This is done through infeasible() virtual function defined in
  BcpsObject.
  After this check we have a list of objects that are infeasible. At this
  point we need to decide what to do with them. Options are (1) generate
  cuts (using generateConstraints()), (2) lift the subproblem (using
  generateVariables()), (3) branch.
  (1) We can generate cuts when infeasible objects are continuous or
  integer. Generate a separating hyperplane that cuts the subproblem solution
  from feasible region.
  (2) Branching when the object is continuous. This is similar to branch
  on variables. Branching procedure should create new subproblems where
  infeasible object has new upper and lower bounds.
*/

/**
   Represents a discrete conic optimization problem.
   Relax same set of rows/columns for the all the subproblems.
 */

class DcoModel: public BcpsModel {
  /// Subproblem solver.
#if defined(__OA__)
  OsiSolverInterface * solver_;
#else
  OsiConicSolverInterface * solver_;
#endif
  ///@name Variable and constraint bounds.
  //@{
  /// Column lower bound.
  double * colLB_;
  /// Column upper bound.
  double * colUB_;
  /// Row lower bound.
  double * rowLB_;
  /// Row upper bound.
  double * rowUB_;
  //@}
  ///@name Number of columns and rows
  //@{
  /// Number of columns.
  int numCols_;
  /// Number of rows (constraints), linear + conic.
  int numRows_;
  /// Number of linear rows.
  int numLinearRows_;
  /// Number of conic rows.
  int numConicRows_;
  //@}
  /// Problem matrix (linear constraints only).
  CoinPackedMatrix * matrix_;
  ///@name Objective function
  //@{
  double objSense_;
  double * objCoef_;
  //@}
  ///@name Column types
  //@{
  /// Number of integer columns in the problem.
  int numIntegerCols_;
  /// Indices of integer columns. Columns are stored in cols_ inherited from
  /// BcpsModel. Size of numIntegerCols_.
  int * integerCols_;
  //@}
  //------------------------------------------------------
  // SOLUTION.
  //------------------------------------------------------
  ///@name Solution related
  //@{
  /// Number of solutions.
  int numSolutions_;
  /// Incumbent objective value.
  double incObjValue_;
  /// Incumbent solution.
  double * incumbentSol_;
  /// Cutoff in lp solver.
  double cutoff_;
  /// Cutoff increment.
  double cutoffInc_;
  /// Current relative optimal gap.
  double currRelGap_;
  /// Current absolute optimal gap.
  double currAbsGap_;
  //@}
  ///@name Variable selection function.
  //@{
  /// Branchs strategy.
  BcpsBranchStrategy * branchStrategy_;
  /// Ramp up branch strategy.
  BcpsBranchStrategy * rampUpBranchStrategy_;
  //@}
  /// Active node.
  AlpsTreeNode * activeNode_;
  ///@name Dco parameters.
  //@{
  /// DisCO parameter.
  DcoParams * dcoPar_;
  //@}
  ///@name Statistics
  //@{
  /// Number of processed nodes.
  int numNodes_;
  /// Number of lp(Simplex) iterations.
  int numIterations_;
  /// Average number of lp iterations to solve a subproblem.
  int aveIterations_;
  //@}
  ///@name Relaxed objects data.
  //@{
  /// Number of relaxed columns
  int numRelaxedCols_;
  /// Array of indices to relaxed columns.
  int * relaxedCols_;
  /// Number of relaxed rows
  int numRelaxedRows_;
  /// Array of indices to relaxed rows.
  int * relaxedRows_;
  //@}
  // Private Functions
  ///@name Read Helpers
  //@{
  /// Add variables to the model. Helps readInstance function.
  void readAddVariables(CoinMpsIO * reader);
  /// Add linear constraints to the model. Helps readInstance function.
  void readAddLinearConstraints(CoinMpsIO * reader);
  /// Add conic constraints to the model. Helps readInstance function.
  void readAddConicConstraints(CoinMpsIO * reader);
  //@}
public:
  ///@name Message printing
  //@{
  /// DisCO message handler.
  CoinMessageHandler * dcoMessageHandler_;
  /// DisCO messages.
  CoinMessages * dcoMessages_;
  //@}
public:
  ///@name Constructors and Destructors
  //@{
  /// Default constructor.
  DcoModel();
  /// Destructor.
  virtual ~DcoModel();
  //@}
  ///@name Solver related
  //@{
  /// Set solver
#if defined(__OA__)
  void setSolver(OsiSolverInterface * solver);
#else
  void setSolver(OsiConicSolverInterface * solver);
#endif
  /// Get solver
#if defined(__OA__)
  OsiSolverInterface * solver() {return solver_;}
#else
  OsiConicSolverInterface * solver() {return solver_;}
#endif
  //@}
  ///@name Other functions
  //@{
  /// Approximate cones.
  void approximateCones();
  //@}
  ///@name Querry problem data
  //@{
  /// Get number of core variables.
  int getNumCoreVariables() const {return numCols_;}
  /// Get number of core linear constraints.
  int getNumCoreLinearConstraints() const {return numLinearRows_;}
  /// Get number of core conic constraints.
  int getNumCoreConicConstraints() const {return numConicRows_;}
  /// Get column lower bounds.
  double * colLB() {return colLB_;}
  /// Get column upper bounds.
  double * colUB() {return colUB_;}
  /// Get row lower bounds.
  double * rowLB() {return rowLB_;}
  /// Get row upper bounds.
  double * rowUB() {return rowUB_;}
  //@}
  ///@name Querry relaxed problem objects
  //@{
  /// Get number of relaxed columns.
  int numRelaxedCols() const {return numRelaxedCols_;}
  /// Get array of indices to relaxed columns.
  int const * relaxedCols() const {return relaxedCols_;}
  /// Get number of relaxed rows
  int numRelaxedRows() const {return numRelaxedRows_;}
  /// Get array of indices to relaxed rows.
  int const * relaxedRows() const {return relaxedRows_;}
  //@}
  ///@name Virtual functions from AlpsModel
  //@{
  /// Read in the problem instance. Currently linear Mps files and Mosek
  /// style conic mps files.
  virtual void readInstance(char const * dataFile);
  /// Do necessary work to make model ready for use, such as classify
  /// variable and constraint types.
  virtual bool setupSelf();
  /// Preprocessing the model. Default does nothing. We do not have any
  /// preprocessing for now.
  virtual void preprocess();
  /// Postprocessing the model. Default does nothing. We do not have any
  /// postprocessing for now.
  virtual void postprocess();
  /// Create the root node.
  virtual AlpsTreeNode * createRoot();
  //@}
};

#endif
