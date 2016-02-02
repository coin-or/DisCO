#ifndef DcoTreeNode_hpp_
#define DcoTreeNode_hpp_

#include <BcpsTreeNode.h>

class DcoTreeNode: public BcpsTreeNode {
public:
  DcoTreeNode();
  virtual ~DcoTreeNode();
  virtual AlpsTreeNode * createNewTreeNode(AlpsNodeDesc*& desc) const;
  virtual void convertToExplicit();
  virtual void convertToRelative();
  virtual int generateConstraints(BcpsModel * model,
				  BcpsConstraintPool * conPool);
  virtual int generateVariables(BcpsModel * model,
				BcpsVariablePool * varPool);
  virtual int chooseBranchingObject(BcpsModel * model);
  virtual int installSubProblem(BcpsModel * model);
  virtual int handleBoundingStatus(int status, bool & keepOn, bool & fathomed);
  virtual int process(bool isRoot = false, bool rampUp = false);
  /** Bounding procedure to estimate quality of this node. */
  virtual int bound(BcpsModel * model);
  /** This method must be invoked on a \c pregnant node (which has all the
      information needed to create the children) and should create the
      children's decriptions. The stati of the children
      can be any of the ones \c process() can return. */
  virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
  branch() = 0;
private:
  // these are disabled in AlpsTreeNode
  DcoTreeNode(DcoTreeNode const & other);
  DcoTreeNode & operator=(DcoTreeNode const & rhs);
};

#endif
