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
  // Why do we need this. It seems unnecessary.
  // CoinPackedMatrix * matrix_;
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
  int numElems_;
  //@}
  /** Objective function. */
  //@{
  double objSense_;
  double * objCoef_;
  //@}
  /** Cone data. */
  //@{
  OsiLorentzCone ** cone_;
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
  /** variable type, 0-continuous 1-binary and 2-general integer */
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
  DcoParams * dcoPar_;
  /** Message handler. */
  CoinMessageHandler * dcoMessageHandler_;
  /** Dco messages. */
  CoinMessages * dcoMessages_;
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
  // ALPS VIRTUAL FUNCTIONS
  /** Read in the problem instance */
  virtual void readInstance(char const * dataFile);
  void approximateCones();
  void createObjects();

  /** Read in Alps parameters. */
  // virtual void readParameters(const int argnum,  char const * const * arglist);
  // /** Write out parameters. */
  // void writeParameters(std::ostream& outstream) const;
  // /** Do necessary work to make model ready for use, such as classify
  //     variable and constraint types.*/
  // virtual bool setupSelf();
  // /** Preprocessing the model. Default does nothing. We do not have any
  //     preprocessing for now. */
  // //virtual void preprocess();
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
