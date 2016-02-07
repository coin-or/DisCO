#include "DcoTreeNode.hpp"
#include "DcoNodeDesc.hpp"

DcoTreeNode::DcoTreeNode() {
}

DcoTreeNode::~DcoTreeNode() {
}

// create tree node from given description
AlpsTreeNode * DcoTreeNode::createNewTreeNode(AlpsNodeDesc *& desc) const {
  //DcoNodeBranchDir branchDir=dynamic_cast<DcoNodeDesc *>(desc)->getBranchedDir();
  int branchInd=dynamic_cast<DcoNodeDesc *>(desc)->branchedInd();
  double lpX = dynamic_cast<DcoNodeDesc *>(desc)->branchedVal();
  DcoModel * model = dynamic_cast<DcoModel*>(desc_->getModel());
  int objInd = model->getIntObjIndices()[branchInd];
  DcoObjectInt * obj = dynamic_cast<DcoObjectInt *>(model->getobjects(objInd));
  double frac = lpX - floor(lpX);
  if (frac > 0.0) {
    model->messageHandler->message(DISCO_NODE_BRANCHONINT,
				   *model->messages())
      << obj->columnIndex() << CoinMessageEol;
  }
  // Set solution estimate for this nodes.
  // double estimate = solEstimate_;
  // if (branchDir == -1) {
  //   estimate -= (1.0-f) * obj->pseudocost().getUpCost();
  // }
  // else {
  //   estimate -= f * obj->pseudocost().getDownCost();
  // }

  // Create a new tree node
  DcoTreeNode * node = new DcoTreeNode(desc);
  desc = NULL;
  return node;
}

void DcoTreeNode::convertToExplicit() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

void DcoTreeNode::convertToRelative() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

int DcoTreeNode::generateConstraints(BcpsModel * model,
				     BcpsConstraintPool * conPool) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::generateVariables(BcpsModel * model,
				   BcpsVariablePool * varPool) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::chooseBranchingObject(BcpsModel * model) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::installSubProblem(BcpsModel * model) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::handleBoundingStatus(int status, bool & keepOn, bool & fathomed) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

int DcoTreeNode::process(bool isRoote, bool rampUp) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

/** Bounding procedure to estimate quality of this node. */
int DcoTreeNode::bound(BcpsModel * model) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return 0;
}

/** This method must be invoked on a \c pregnant node (which has all the
    information needed to create the children) and should create the
    children's decriptions. The stati of the children
    can be any of the ones \c process() can return. */
std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
DcoTreeNode::branch() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > res;
  return res;
}
