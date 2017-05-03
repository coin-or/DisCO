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
  This class represents a node of the branch and bound tree. This class is a
  fundamental and of highest impartance to understand how DisCO works.

  # Inheritance
   DcoTreeNode inherits BcpsTreeNode. BcpsTreeNode inherits AlpsTreeNode.
   DcoTreeNode -> BcpsTreeNode -> AlpsTreeNode.

  ## Explanation of inherited classes.
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
        BcpsTreeNode has one field only and it is branchObject_. It defines
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

   # Remarks on some important functions

   <ul>
        <li> process()<br>

        This function is called by Alps when it decides to process this
        node. Processing might solve the subproblem represented by the node,
        generating constraints/variables for the discrete problem master
        problem).

        This function should be implemented Bcps level using abstract functions
        (branchConstraintOrPrice(), generateConstraints(),
        generateVariables()). Bcps is already designed that way, but we are not
        using that, but we will be there (date: 05/03/2016, author: aykut).

        <li> boundingLoop()<br>

        This function is called from process(). After solving the subproblem we
        might generate constraints/variables for the problem and solve it
        again. This improves the bound of the subproblem and hopefully bound
        for the master problem. This loops carries this solve, geenrate cuts,
        generate columns, commands many times.

        <li> branchConstrainOrPrice()<br>

        This function is called from boundingLoop() function. During bounding
        loop, this function decides what to do next, i.e., keep bounding,
        branch, generate constraints, generate columns.

        This function is a rename/replacement for handleBoundStatus defined in
        Alps. It is not used that way yet, but has to be (date:
        05/03/2016, author: aykut).

        <li> generateConstraints()<br>

        This function is called if branchConstrainOrPrice() decides generating
        constraints. It generates constraints for the subproblem being solved
        (locally valid cuts), or master problem (globally valid cuts).

        <li> generateVariables()<br>

        This function is called if branchConstrainOrPrice() decides generating
        variables/columns (aka lifting the problem). It generates columns for
        the subproblem being solved (locally valid), or master problem
        (globally valid).

        <li> chooseBranchingObject()<br>

        <li> installSubProblem()<br>

        Installs the subproblem being represented by the node to the solver
        interface.

        <li> bound()<br>

        Solves the subproblem.
   </ul>

   # Generating and applying constraints

   generateConstraints() populates a BcpsConstraintPool
   object. applyConstraints() function takes this object as input and adds cuts
   to the solver.


   # Bcps ideas

   Why generateConstraints in Bcps level returns an int?
   Disco::generateConstraints returns 0 for now.

   applyConstraints() should be in the Bcps level.

   In generateConstraints() Cgl Cut generators create OsiCut instances. We
   convert them to DcoConstraint. In applyConstraints() we convert DcoConstraint
   instances to OsiCut instances to add them to the solver. Why don't we keep
   them in OsiCut instances all along.

   # Alps ideas/questions.

   What happens to fathomed nodes?

   When we create nodes we set quality to the parents objective value. Does it
   ever get updated once node is solved? After cuts are generated?
   Ans: it is updated in the DcoTreeNode::bound() function.

 */



// todo(aykut): how about adding a getBasis() and setBasis() to get basis from
// desc_ and set basis in desc_? to DcoTreeNode class.
// todo(aykut): Bcps has a generic process method? Should we use it?
//

class DcoTreeNode: public BcpsTreeNode {
  /// Keeps statistics to help deciding branch/cut/price.
  struct BcpStats {
    // number of bound iterations
    int numBoundIter_;
    // total number of cuts added.
    int numTotalCuts_;
    // number of cuts in the last time
    int numLastCuts_;
    // total number of OA cuts added
    int numOaCuts_;
    // total number of milp cuts added
    int numMilpCuts_;
    // total objective improvement
    double totalImp_;
    // objective improvement since last bounding
    double lastImp_;
    // objective improvement from last bounding
    double lastObjVal_;
    // objective value at start
    double startObjVal_;
    // number of bunding iteraton for milp cuts
    int numMilpIter_;
    // how many times the cut in the current solver was inactive?
    // indices of cuts are model->numLinearRows(), ... ,
    // model->solver()->getNumRows()-1. inactive_[i] is the number of times cut
    // i is inactive. cut i is the current i+model->numLinearRows() th row of
    // solver.
    std::list<int> inactive_;
    // generator of cuts. generatorIndex_[i] returns the index of the cut generator
    // in model()->conGenerators_, i.e.,
    // model()->conGenerators()[generatorIndex_[i]] is the generator of the cut
    // sitting at index i.
    std::list<int> generatorIndex_;
  };
  BcpStats bcpStats_;
  /// Decide whether the given cut generator should be used, based on the cut
  /// strategy.
  void decide_using_cg(bool & do_use, DcoConGenerator * cg,
                       int numRowsInf, int numColsInf) const;
  /// Copies node description of this node to given child node.
  /// New node is explicitly stored in the memory (no differencing).
  void copyFullNode(DcoNodeDesc * child_node) const;
  /// Sets node status to pregnant and carries necessary operations.
  void processSetPregnant();
  /// This function is called after bound method is called. It checks solver
  /// status.
  void afterBound(DcoSubproblemStatus subproblem_status);
  int boundingLoop(bool isRoot, bool rampUp);
  /// find number of infeasible integer variables.
  void checkRelaxedCols(int & numInf);
  /// update cut stats and clean in necessary
  void checkCuts();
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
  // todo(aykut): We do not use the implementation given by BcpsTreeNode yet, but that is
  // the future.
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
  virtual int generateConstraints(BcpsConstraintPool * conPool);
  /// Generate variables (lift the problem) and store them in the given
  /// variable pool.
  virtual int generateVariables(BcpsVariablePool * varPool);
  /// Choose a branching object.
  virtual int chooseBranchingObject();
  /// Install subproblem corresponding to this node to the solver.
  virtual int installSubProblem();
  /// decide what to do, keepBounding?, generate constraints?, generate
  /// variables?. Fromerly known handleBoundStatus
  virtual void branchConstrainOrPrice(BcpsSubproblemStatus subproblem_status,
                                      bool & keepBounding,
                                      bool & branch,
                                      bool & generateConstraints,
                                      bool & generateVariables);
  /// Bounds the problem by solving the subproblem corresponds to this node.
  virtual BcpsSubproblemStatus bound();
  /// Call heuristics to search for solutions.
  virtual void callHeuristics();
  /// Apply given constraints
  virtual void applyConstraints(BcpsConstraintPool const * conPool);
  //@}

  ///@name Other functions
  //@{
  /// Get node description. Overwrites the one inherited from AlpsTreeNode.
  DcoNodeDesc * getDesc() const;
  //@}

  ///@name Encode and Decode functions for parallel execution
  //@{
  /// Get encode from #AlpsKnowledge
  using AlpsKnowledge::encode;
  /// Pack this into an encoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// Unpack into this from an encoded object.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
  /// Unpack into a new DcoTreeNode object and return a
  /// pointer to it.
  virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
  //@}

private:
  // These are disabled in AlpsTreeNode
  DcoTreeNode(DcoTreeNode const & other);
  DcoTreeNode & operator=(DcoTreeNode const & rhs);
};

#endif
