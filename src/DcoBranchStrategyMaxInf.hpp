#ifndef DcoBranchStrategyMaxInf_hpp_
#define DcoBranchStrategyMaxInf_hpp_

#include <BcpsBranchStrategy.h>

class DcoModel;

/**
   # How does branching work?

   When a node is picked from the pool and it is pregnant Alps calls
   AlpsSearchStrategy::createNewNodes(). AlpsSearchStrategy::createNewNodes()
   calls AlpsNode::branch() and Alps process output of this function further.

   This function
   Call DcoBranchStrategy[Name]::createCandBranchObjects() for the
   corresponding strategy. This will populate branchObjects_ inherited from
   BcpsBranchStrategy (where does bestBranchObject_ get into play?).

   Then we create child nodes using the bestBranchObject_ in
   DcoTreeNode::branch(). return them, Alps will add them into the node pool.

   Where does DcoVariable::createBranchObject get into play?

   This class represents maximum infeasibility branching strategy. Inherits
   BcpsBranchStrategy.

   BcpsBranchStrategy:
   BcpsBranchStrategy is an abstract base class for defining a branch stragey.
   It provides an interface for creating a set of branching candidates and
   comparing them.

   DcoBranchStrategyMaxInf:
   Implements maximum infeasibility branching in interface inherited from
   BcpsBranchStrategy.

   What is a branch object?
   Keeps necessary information for a branching to be perfomed.

   BcpsBranchObject contains the member data required when choosing
   branching entities and excuting actual branching. It also has
   the member funtions to do branching by adjusting bounds, etc.
   in solver. Branching objects can be simple integer variables or more
   complicated objects like SOS.

   BcpsBranchObject fields are type_, model_, objectIndex_, upScore_,
   downScore_, direction_, value_, numBranchesLeft_.

   virtual functions are: clone(), numBranches(), numBranchesLeft(), branch(),
   print(), boundBranch(),

 */
class DcoBranchStrategyMaxInf: public BcpsBranchStrategy {
public:
  DcoBranchStrategyMaxInf();
  DcoBranchStrategyMaxInf(DcoBranchStrategyMaxInf const & other);
  DcoBranchStrategyMaxInf(DcoModel * model);
  virtual ~DcoBranchStrategyMaxInf();
  virtual BcpsBranchStrategy * clone() const;
  virtual int createCandBranchObjects(int numPassesLeft, double ub);
  /** Compare branching object thisOne to bestSoFar. If thisOne is better than
      bestObject, return branching direction(1 or -1), otherwise return 0.  If
      bestSorFar is NULL, then always return branching direction(1 or -1).
  */
  virtual int betterBranchObject(BcpsBranchObject * b,
                                 BcpsBranchObject * bestSoFar);
private:
  DcoBranchStrategyMaxInf & operator=(DcoBranchStrategyMaxInf const & rhs);
};

#endif
