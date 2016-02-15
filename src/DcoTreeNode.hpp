#ifndef DcoTreeNode_hpp_
#define DcoTreeNode_hpp_

// todo(aykut): why do we need the following two includes?
// this should be fixed in Alps level.
#include <Alps.h>
#include <AlpsEncoded.h>
#include <AlpsNodeDesc.h>
// #include <BcpsNodeDesc.h>

#include <BcpsTreeNode.h>
#include "DcoNodeDesc.hpp"
#include "DcoModel.hpp"

/*!
   DcoTreeNode inherits BcpsTreeNode. BcpsTreeNode inherits AlpsTreeNode.
   DcoTreeNode -> BcpsTreeNode -> AlpsTreeNode.

   <ul>
   <li> AlpsTreeNode<br>
	AlpTreeNode has many fields. These fields are active_, index_, depth_,
	solEstimate_, quality_, parent_, parentIndex_, numChildren_, children_,
	explicit_, desc_, status_, knowledgeBroker_, sentMark_, diving_. It
	represents for a node in any kind of three search algorithm. It keeps
	node information common to all tree search algorithm.

	Here desc_ is an important field. It keeps a pointer to the node
	description, it is an instance of AlpsNodeDesc class. Projects developed
	on top of Alps can create its own node description by inheriting
	AlpsNodeDesc. Then an AlpsTreeNode will store a pointer to the user defined
	node description.

	Another important field is status_. It is of type AlpsNodeStatus (enum
	type). Values for status are
	<ol>
	<li> AlpsNodeStatusCandidate<br>
	     Indicates node is a candidate to process. These
	     nodes are freshly created nodes and utouched yet.
	<li> AlpsNodeStatusEvaluated<br>
	     Node is processed, but it needs further processing. This happens when
	     a node is selected, cuts are generated for it and then it is inserted
	     to the node pool back. The node with this status is not branched yet.
	<li> AlpsNodeStatusPregnant<br>
	     Node processing is done and it should be branched.
	<li> AlpsNodeStatusBranched<br>
	     Node is processed and branched.
	<li> AlpsNodeStatusFathomed<br>
	     Node is fathomed. It is inferior and will not lead a solution.
	<li> AlpsNodeStatusDiscarded<br>
	     Node is discarded. It is inferior and will not lead a solution.
	</ol>

   <li> BcpsTreeNode<br>
	BcpsTreeNode has one field only and ot is branchObject_. It defines
	many virtual functions. These functions are to be implemented in the
	projects implemented on top of Bcps. These functions are listed as follows,
	<ol>
	<li> generateConstraints()<br>
	     To generate constraints (i.e., cuts) for the corresponding subproblem.
	<li> generateVariables()<br>
	     To generate constraints (lifting) for the corresponding subproblem.
	<li> chooseBranchingObject()<br>
	<li> installSubProblem()<br>
	     Extract node information (bounds, constraints, variables) from
	     this node and load the information into the relaxation solver,
	     such as linear programming solver.
	<li> handleBoundingStatus()<br>
	     This is completely ignored in Blis. I am not sure whether I will use it in DisCO.
	     Explanation in Bcps documentation is as follows,

	     Handle bounding status: (1) relaxed feasible but not integer feasible,
	     (2) integer feasible, (3) infeasible, (4) unbounded,
	     (5) fathomed (for instance, reaching objective limit for LP).

	     Set node status accordingly.
	       \param staus    Input   The solution status of bounding.
	       \param keepOn   Output  Whether to keep on bounding.
	       \param fathomed Output  Whether this node is fathomed.

	<li> process()<br>
	     This methods performs the processing of the node. For branch and cut,
	     this would mean generating cuts, performing bounding operation and
	     branching. This function might update the node status kept in status_
	     (inherited from AlpsTreeNode). Alps search will decide what to do
	     with this node by checking its status.

	     todo(aykut): check the return type of process method. Does it make sense?
	     Does Alps use this return status in any way?

	<li> bound()<br>
	     Solve the corresponding subproblem (continuous relaxation in this
	     case).
	     In Blis this function also updates the pseudocosts for the corresponding
	     integer variables.
	<li> branch()<br>
	     This method must be invoked on a \c pregnant (AlpsNodeStatusPregnant)
	     node (which has all the information needed to create the children)
	     and should create the children's decriptions.
	</ol>
   </ul>

 */

// todo(aykut): how about adding a getBasis() and setBasis() to get basis from
// desc_ and set basis in desc_? to DcoTreeNode class.
// todo(aykut): Bcps has a generic process method? Should we use it?
//

class DcoTreeNode: public BcpsTreeNode {
public:
  ///@name Constructors and Destructors
  //@{
  /// Default constructor.
  DcoTreeNode();
  /// Construct for a given description.
  DcoTreeNode(AlpsNodeDesc * & desc);
  /// Destructor.
  virtual ~DcoTreeNode();
  //@}
  ///@name Virtual functions inherited from AlpsTreeNode
  //@{
  /// Create new tree nodes from the given description.
  virtual AlpsTreeNode * createNewTreeNode(AlpsNodeDesc *& desc) const;
  /// Convert node description to explicit.
  virtual void convertToExplicit();
  /// Convert node description to relative.
  virtual void convertToRelative();
  /// Process node. Alps calls this function when it is picked from node pool.
  virtual int process(bool isRoot=false, bool rampUp=false);
  /// Branch this node. Alps calls this function to create children from this
  /// node. Alps calls this function when the node status is
  /// AlpsNodeStatusPregnant.
  virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
  branch();
  //@}
  ///@name Virtual functions inherited from BcpsTreeNode
  //@{
  /// Generate constraints (cuts) and store them in the given constraint pool.
  virtual int generateConstraints(BcpsModel * bcps_model,
				  BcpsConstraintPool * conPool);
  /// Generate variables (lift the problem) and store them in the given
  /// variable pool.
  virtual int generateVariables(BcpsModel * model,
				BcpsVariablePool * varPool);
  /// Choose a branching object.
  virtual int chooseBranchingObject(BcpsModel * bcps_model);
  /// Install subproblem to the solver.
  virtual int installSubProblem(BcpsModel * bcps_model);
  // todo(aykut) Should we use this?
  //virtual int handleBoundingStatus(int status, bool & keepOn, bool & fathomed);
  /// Bounds the problem by solving the subproblem corresponds to this node.
  virtual int bound(BcpsModel * model);
  //@}
  ///@name Other functions
  //@{
  /// Get node description. Overwrites the one inherited from AlpsTreeNode.
  DcoNodeDesc * getDesc() const;
  /// Get model this node belongs. Overwrites the one inherited from AlpsTreeNode.
  DcoModel * getModel() const;
  //@}
private:
  // these are disabled in AlpsTreeNode
  DcoTreeNode(DcoTreeNode const & other);
  DcoTreeNode & operator=(DcoTreeNode const & rhs);
  /// Sets node status to pregnant and carries necessary operations.
  void processSetPregnant();
  /// This function is called after bound method is called. It checks solver
  /// status.
  void afterBound(DcoSubproblemStatus subproblem_status);
  int boundingLoop(bool isRoot, bool rampUp);
  /// decide what to do, keepBounding?, generate constraints?, generate
  /// variables?. Fromerly known handleBoundStatus
  void decide(DcoSubproblemStatus subproblem_status,
	 bool & generateConstraints,
	 bool & generateVariables);
  /// find number of infeasible integer variables.
  void checkRelaxedCols(int & numInf);
};

#endif
