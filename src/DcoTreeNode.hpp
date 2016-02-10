#ifndef DcoTreeNode_hpp_
#define DcoTreeNode_hpp_

// todo(aykut): why do we need the following two includes?
// this should be fixed in Alps level.
#include <Alps.h>
#include <AlpsEncoded.h>
#include <AlpsNodeDesc.h>
// #include <BcpsNodeDesc.h>

#include <BcpsTreeNode.h>


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


class DcoTreeNode: public BcpsTreeNode {
public:
  DcoTreeNode();
  DcoTreeNode(AlpsNodeDesc * & desc);
  virtual ~DcoTreeNode();
  virtual AlpsTreeNode * createNewTreeNode(AlpsNodeDesc *& desc) const;
  virtual void convertToExplicit();
  virtual void convertToRelative();
  virtual int generateConstraints(BcpsModel * bcps_model,
				  BcpsConstraintPool * conPool);
  virtual int generateVariables(BcpsModel * model,
				BcpsVariablePool * varPool);
  virtual int chooseBranchingObject(BcpsModel * bcps_model);
  virtual int installSubProblem(BcpsModel * bcps_model);
  virtual int handleBoundingStatus(int status, bool & keepOn, bool & fathomed);
  virtual int process(bool isRoot=false, bool rampUp=false);
  /** Bounding procedure to estimate quality of this node. */
  virtual int bound(BcpsModel * model);
  /** This method must be invoked on a \c pregnant node (which has all the
      information needed to create the children) and should create the
      children's decriptions. The stati of the children
      can be any of the ones \c process() can return. */
  virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
  branch();
private:
  // these are disabled in AlpsTreeNode
  DcoTreeNode(DcoTreeNode const & other);
  DcoTreeNode & operator=(DcoTreeNode const & rhs);
};

#endif
