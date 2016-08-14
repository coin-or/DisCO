#ifndef DcoModel_hpp_
#define DcoModel_hpp_

#include <BcpsModel.h>
#include <OsiConicSolverInterface.hpp>
#include <OsiLorentzCone.hpp>
#include <BcpsBranchStrategy.h>


#include "DcoParams.hpp"
#include "DcoConstraint.hpp"

class DcoConGenerator;
class DcoSolution;
class DcoHeuristic;

class CglCutGenerator;
class CglConicCutGenerator;

/**
   Represents a discrete conic optimization problem (master problem).
   Some set of rows/columns will be relaxed in this problem to get subproblems
   represented by the branch and bound tree nodes.

   # Fields of DcoModel
   relaxedCols_ keeps the set of relaxed columns. relaxedRows_ keeps the set
   of relaxed rows. We relax integer columns only, their integrality
   constraint is relaxed. When OA algorithm is used we relax rows corresponding
   to conic constraints.

   In DcoModel columns are stored in cols_ inherited from BcpsModel. We keep
   number of integer variables at numIntegerCols_ and their indices at (int
   * intColIndices_). intColIndices_[0] gives the index of the first integer
   column in the cols_ array.

   DcoModel objects are kept at constraints_ and variables_ inherited from
   BcpsModel.

   In Blis (MILP solver built on top of Bcps), integer variables have their own
   class, BlisObjectInt.  BlisObjectInt inherits BcpsObject class.

   # Cut generation

   Pointers to cut generators are stored in conGenerators_. Type of generators
   are DcoConGenerator. DcoConGenerator is an abstract base class (ABC) for
   constraint generators.  Two different classes implements this ABC,
   DcoLinearConGenerator and DcoConicConGenerator. DcoConGenerator has a single
   pure virtual function, generateConstraints(conPool). This function generates
   constraints and add them to the given pool.

   DcoLinearConGenerator implements linear cuts. CglCutGenerator is given as an
   input to constructor.

   DcoConicConGenerator implements generating supports for conic
   constraints. Constructor takes a CglConicCutGenerator as an input.

   # Heuristics

   # setupSelf()

   In serial code this function is called after readInstance() method.

   In parallel code it is called after readInstance() in the master
   processor. In other processors it is called after the received encoded
   object (AlpsEncoded instance) is decoded to self. This decode part covers
   the following fields, variables_, constraints_, dcoPar_ and objSense_.

   setupSelf() should genrate all fields of DcoModel from the 4 data fields
   mentioned above (variables_, constraints_, dcoPar_ and objSense_).

   This function also loads the problem defined by self to the solver. It
   sets/creates branching strategy, cut generator and heuristics object.

*/

class DcoModel: public BcpsModel {
  /// Subproblem solver.
#if defined(__OA__)
  OsiSolverInterface * solver_;
#else
  OsiConicSolverInterface * solver_;
#endif

  ///==========================================================================
  /// Fields that will be set by ::readInstance() and sent to other processors
  /// through network. ::setupSelf() will use these fields to set the rest.
  ///==========================================================================
  ///@name Variable and constraint bounds.
  //@{
  // colLB_ and colUB_ used when installing subproblems in
  // each node. Having it a class member we do not need to allocate and delete
  // memory every time a subproblem is installed to the solver in a node.
  /// Column lower bound, corresponds to the last subproblem installed.
  double * colLB_;
  /// Column upper bound, corresponds to the last subproblem installed.
  double * colUB_;
  /// Row lower bound.
  double * rowLB_;
  /// Row upper bound.
  double * rowUB_;
  //@}

  ///@name Constraint matrix (linear part) and conic constraints. These can be
  /// freed at the end of ::setupSelf()
  //@{
  /// Constraint matrix.
  CoinPackedMatrix * matrix_;
  /// We keep cones in basic form for now, it is easier to send/receive
  int * coneStart_;
  int * coneMembers_;
  int * coneType_;
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
  int * isInteger_;
  //@}
  ///==========================================================================

  ///==========================================================================
  /// Fields that will be set by ::setupSelf(). constraints_ and variables_
  /// inherited from BcpsModel will also set by ::setupSelf().
  ///==========================================================================
  ///@name Variable selection function.
  //@{
  /// Branchs strategy.
  BcpsBranchStrategy * branchStrategy_;
  /// Ramp up branch strategy.
  BcpsBranchStrategy * rampUpBranchStrategy_;
  //@}

  ///@name Dco parameters.
  //@{
  /// DisCO parameter.
  DcoParams * dcoPar_;
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

  ///@name Heuristics
  //@{
  DcoHeurStrategy heurStrategy_;
  int heurFrequency_;
  std::vector<DcoHeuristic*> heuristics_;
  //@}

  ///@name Cut generator related.
  //@{
  /// global cut strategy, it will be set with respect to specific cut
  /// strategies. It will be set to the most allowing one, ie. if we have
  /// strategies with root and periodic calls, it will be set to periodic.
  DcoCutStrategy cutStrategy_;
  /// Cut generation frequency, it will be set with respect to specific cut
  /// strategies. It will be set to the most frequent one, ie. if we have
  /// strategies with frequencies 10 and 20, it will be set to 10.
  int cutGenerationFrequency_;
  /// Constraint generators.
  std::vector<DcoConGenerator*> conGenerators_;
  //@}
  ///==========================================================================


  // Private Functions
  ///@name Read Helpers
  //@{
  /// Add variables to the model. Helps readInstance function.
  void setupAddVariables();
  /// Add linear constraints to the model. Helps readInstance function.
  void setupAddLinearConstraints();
  /// Add conic constraints to the model. Helps readInstance function.
  void setupAddConicConstraints();
  //@}

  ///@name Setup Helpers
  //@{
  /// Set log levels, Alps, Bcps and Disco
  void setMessageLevel();
  /// Set branching strategy from parameters.
  void setBranchingStrategy();
  /// Add constraint generators with respect to parameters.
  void addConstraintGenerators();
  /// Add heuristics
  void addHeuristics();
  //@}

  /// write parameters to oustream
  void writeParameters(std::ostream& outstream) const;

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
  /// return to branch strategy.
  BcpsBranchStrategy * branchStrategy() {return branchStrategy_;}
  /// return Dco Parameter
  DcoParams const * dcoPar() const {return dcoPar_;}
  /// get upper bound of the objective value for minimization
  double bestQuality();
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
  /// Get objective sense, 1 for min, -1 for max
  double objSense() const { return objSense_; }
  /// Get number of integer variables.
  int numIntegerCols() const { return numIntegerCols_; }
  /// Get indices of integer variables. Size of numIntegerCols().
  int const * integerCols() const { return integerCols_; }
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

  ///@name Constraint Generation related.
  //@{
  /// Add constraint generator using linear Cgl.
  void addConGenerator(CglCutGenerator * cgl_gen, char const * name,
                       DcoCutStrategy dco_strategy, int frequency);
  /// Add constraint generator using conic Cgl.
  void addConGenerator(CglConicCutGenerator * cgl_gen, char const * name,
                       DcoCutStrategy dco_strategy, int frequency);
  /// Add constraint generator.
  void addConGenerator(DcoConGenerator * dco_gen);
  /// Get the number of constraint generators.
  int numConGenerators() const { return conGenerators_.size(); }
  /// Get a specific constraint generator.
  DcoConGenerator * conGenerators(int i) const { return conGenerators_[i]; }
  /// Get global cut strategy. It will be set using specific cut strategies, to
  /// the most allowing one. If we have strategies with root and periodic
  /// calls, it will be set to periodic.
  DcoCutStrategy cutStrategy() const {return cutStrategy_;}
  /// Set global cut strategy. It will be set using specific cut strategies, to
  /// the most allowing one. If we have strategies with root and periodic
  /// calls, it will be set to periodic.
  void setCutStrategy(DcoCutStrategy strategy) {cutStrategy_ = strategy;}
  //@}

  ///@name Heuristics related
  //@{
  // get number of heuristics
  int numHeuristics() const { return heuristics_.size(); }
  // get a constant specific heuristic, for reading statistics.
  DcoHeuristic const * heuristics(int i) const { return heuristics_[i]; }
  // get a specific heuristic, for solution search
  DcoHeuristic * heuristics(int i) { return heuristics_[i]; }
  //@}


  /// Check feasiblity of subproblem solution, store number of infeasible
  /// columns and rows.
  virtual DcoSolution * feasibleSolution(int & numInfColumns,
                                         int & numInfRows);

  ///@name Virtual functions from AlpsModel
  //@{
  /// Read in the problem instance. Currently linear Mps files and Mosek
  /// style conic mps files.
  virtual void readInstance(char const * dataFile);
  /// Reads in parameters.
  /// This function is called from AlpsKnowledgeBrokerSerial::initializeSearch
  /// It reads and stores the parameters in alpsPar_ inherited from AlpsModel.
  virtual void readParameters(int const argnum, char const * const * arglist);
  /// Do necessary work to make model ready for use, such as classify
  /// variable and constraint types.
  /// Called from AlpsKnowledgeBrokerSerial::initializeSearch. Called
  /// after readParameters and preprocess.
  virtual bool setupSelf();
  /// Preprocessing the model. Default does nothing. We do not have any
  /// preprocessing for now.
  /// Called from AlpsKnowledgeBrokerSerial::initializeSearch. Called
  /// after readParameters and before setupSelf.
  virtual void preprocess();
  /// Postprocessing the model. Default does nothing. We do not have any
  /// postprocessing for now.
  virtual void postprocess();
  /// Create the root node.
  virtual AlpsTreeNode * createRoot();
  /** This function is called every time the node counts hits
     AlpsParams::intParams::nodeLogInterval. It prints information related to
     search status. In parallel mode only master should log. */
  virtual void nodeLog(AlpsTreeNode * node, bool force);
  /// This is called at the end of the AlpsKnowledgeBroker::rootSearch
  /// Prints solution statistics
  virtual void modelLog();
  //@}

  ///@name Encode and Decode functions
  //@{

  // note(aykut): It is enough to encode the DcoModel fields that are minimal
  // (less network communication).  setupSelf() will be called by Alps to set
  // up the rest of the fields that can be constructed/computed from the set
  // ones.

  // This grabs function "#AlpsEncoded * AlpsKnowledge::encoding() const"
  // inherited from #AlpsKnowledge. It will not get into overload resoulution
  // since we declare "#AlpsEncoded * encode() const" here.
  using AlpsKnowledge::encode;
  /// The method that encodes the this instance of model into the given
  /// #AlpsEncoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// The method that decodes the given #AlpsEncoded object into a new #DcoModel
  /// instance and returns a pointer to it.
  virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
  /// The method that decodes this instance from the given #AlpsEncoded object.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
  //@}

  /// report feasibility of the best solution
  void reportFeasibility();

};

#endif
