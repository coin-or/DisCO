#include "DcoBranchStrategyRel.hpp"
#include "DcoModel.hpp"

DcoBranchStrategyRel::DcoBranchStrategyRel() {
}

DcoBranchStrategyRel::DcoBranchStrategyRel(DcoBranchStrategyRel const & other) {
}

DcoBranchStrategyRel::DcoBranchStrategyRel(DcoModel * model):
  BcpsBranchStrategy(model) {
}

DcoBranchStrategyRel::~DcoBranchStrategyRel() {
}

BcpsBranchStrategy * DcoBranchStrategyRel::clone() const {
  return new DcoBranchStrategyRel(*this);
}

int DcoBranchStrategyRel::createCandBranchObjects(int numPassesLeft,
						  double ub) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

int DcoBranchStrategyRel::betterBranchObject(BcpsBranchObject * b,
					     BcpsBranchObject * bestSoFar) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}
