#ifndef DcoModel_hpp_
#define DcoModel_hpp_

#include <BcpsModel.h>
#include <OsiConicSolverInterface.hpp>
#include <OsiLorentzCone.hpp>
#include <BcpsBranchStrategy.h>


#include "DcoParams.hpp"
#include "DcoConstraint.hpp"

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
  /** Variable and constraint bounds. */
  //@{
  double * colLB_;
  double * colUB_;
  double * rowLB_;
  double * rowUB_;
  //@}
  /** Number of columns/rows/cones/elements */
  //@{
  int numCols_;
  int numRows_;
  int numCones_;
  int numElems_;
  //@}
  /** Problem matrix (linear constraints only) */
  //@{
  CoinPackedMatrix * matrix_;
  //@}
  /** Objective function. */
  //@{
  double objSense_;
  double * objCoef_;
  //@}
  /** Column types. */
  //@{
  // number of integer columns in the problem
  int numIntegerCols_;
  // indices of integer columns
  // columns are stored in cols_ inherited from BcpsModel
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
  /** Starting var/con bounds for processing each node */
  //@{
  // double * startColLB_;
  // double * startColUB_;
  // double * startRowLB_;
  // double * startRowUB_;
  //@}
  /** Variable selection function. */
  BcpsBranchStrategy * branchStrategy_;
  BcpsBranchStrategy * rampUpBranchStrategy_;
  /** Active node. */
  AlpsTreeNode * activeNode_;
  //------------------------------------------------------
  // PARAMETERS, STATISTICS, and MESSAGE
  //------------------------------------------------------
  /** Dco parameters. */
  DcoParams * dcoPar_;
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
  // todo(aykut) why are we keeping this in the DcoModel class.
  // It seems it is ony used in DcoTreeNode::installSubProblem()
  /** Temporary store old cuts at a node when installing a node. */
  DcoConstraint **oldConstraints_;
  /** The memory size allocated for oldConstraints_. */
  int oldConstraintsSize_;
  /** Number of old constraints. */
  int numOldConstraints_;
  /// Add variables to the model. Helps readInstance function.
  void readAddVariables(CoinMpsIO * reader);
  /// Add linear constraints to the model. Helps readInstance function.
  void readAddLinearConstraints(CoinMpsIO * reader);
  /// Add conic constraints to the model. Helps readInstance function.
  void readAddConicConstraints(CoinMpsIO * reader);
public:
  // PUBLIC DATA FIELDS
    /** Message handler. */
  CoinMessageHandler * dcoMessageHandler_;
  /** Dco messages. */
  CoinMessages * dcoMessages_;
public:
  DcoModel();
  virtual ~DcoModel();
#if defined(__OA__)
  void setSolver(OsiSolverInterface * solver);
#else
  void setSolver(OsiConicSolverInterface * solver);
#endif
#if defined(__OA__)
  OsiSolverInterface * solver() {return solver_;}
#else
  OsiConicSolverInterface * solver() {return solver_;}
#endif
  void approximateCones();
  int getNumCoreVariables() const {return numCols_;}
  int getNumCoreLinearConstraints() const {return numRows_;}
  //@{
  /** Get number of old constraints. */
  int getNumOldConstraints() const { return numOldConstraints_; }
  /** Set number of old constraints. */
  void setNumOldConstraints(int num) { numOldConstraints_ = num; }
  /** Get max number of old constraints. */
  int getOldConstraintsSize() const { return oldConstraintsSize_; }
  /** Set max number of old constraints. */
  void setOldConstraintsSize(int num) { oldConstraintsSize_ = num; }
  /** Access old constraints. */
  DcoConstraint **oldConstraints() { return oldConstraints_; }
  /** set old constraints. */
  void setOldConstraints(DcoConstraint **old) { oldConstraints_ = old; }
  /** Set max number of old constraints. */
  void delOldConstraints() {
    delete [] oldConstraints_;
    oldConstraints_ = NULL;
  }
  //@}
  double * colLB() {return colLB_;}
  double * colUB() {return colUB_;}
  double * rowLB() {return rowLB_;}
  double * rowUB() {return rowUB_;}


  // ALPS VIRTUAL FUNCTIONS
  /** Read in the problem instance */
  virtual void readInstance(char const * dataFile);
  /** Read in Alps parameters. */
  // virtual void readParameters(const int argnum,  char const * const * arglist);
  // /** Write out parameters. */
  // void writeParameters(std::ostream& outstream) const;
  // /** Do necessary work to make model ready for use, such as classify
  //     variable and constraint types.*/
  // virtual bool setupSelf();
  /** Preprocessing the model. Default does nothing. We do not have any
      preprocessing for now. */
  virtual void preprocess();
  // /** Postprocessing the model. Default does nothing. We do not have any
  //     postprocessing for now. */
  // //virtual void postprocess();

  /** Create the root node. */
  virtual AlpsTreeNode * createRoot();

  // /** Problem specific log. */
  // virtual void modelLog() {}

  // /** Node log. */
  // virtual void nodeLog(AlpsTreeNode *node, bool force);

  // /** Return true if all nodes on this process can be fathomed.*/
  // virtual bool fathomAllNodes() { return false; }
  // /** Decode model data from the encoded form and fill member data.*/
  // virtual void decodeToSelf(AlpsEncoded& encoded) {}

  // /** Register knowledge class. */
  // virtual void registerKnowledge() { /* Default does nothing */ }

  // /** Send generated knowledge */
  // virtual void sendGeneratedKnowledge() { /* Default does nothing */ }

  // /** Receive generated knowledge */
  // virtual void receiveGeneratedKnowledge() { /* Default does nothing */ }

  // /** Pack knowledge to be shared with others into an encoded object.
  //     Return NULL means that no knowledge can be shared. */
  // virtual AlpsEncoded* packSharedKnowlege() {
  //   /* Default does nothing */
  //   AlpsEncoded* encoded = NULL;
  //   return encoded;
  // }

  // /** Unpack and store shared knowledge from an encoded object. */
  // virtual void unpackSharedKnowledge(AlpsEncoded&)
  // { /* Default does nothing */ }
  // // end of ALPS VIRTUAL FUNCTIONS








  // BCPS VIRTUAL FUNCTIONS




};

#endif
