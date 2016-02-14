#include "DcoBranchStrategyMaxInf.hpp"
#include "DcoModel.hpp"

DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf() {
}

DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf(DcoBranchStrategyMaxInf const & other) {
}

DcoBranchStrategyMaxInf::DcoBranchStrategyMaxInf(DcoModel * model):
  BcpsBranchStrategy(model) {
}

DcoBranchStrategyMaxInf::~DcoBranchStrategyMaxInf() {
}

BcpsBranchStrategy * DcoBranchStrategyMaxInf::clone() const {
  return new DcoBranchStrategyMaxInf(*this);
}

int DcoBranchStrategyMaxInf::createCandBranchObjects(int numPassesLeft,
						  double ub) {
  // get model
  DcoModel * model = dynamic_cast<DcoModel*>(model_);
  // get number of relaxed objects
  int num_relaxed = model->numRelaxedCols();
  // get indices of relaxed object
  int const * relaxed = model->relaxedCols();
  // store branch objects in bobjects
  std::vector<BcpsBranchObject*> bobjects;
  // iterate over relaxed columns and populate bobjects
  for (int i=0; i<num_relaxed; ++i) {
    int preferredDir;
    BcpsObject * curr_object = model->getVariables()[relaxed[i]];
    double infeasibility = curr_object->infeasibility(model, preferredDir);
    // check the amount of infeasibility
    if (infeasibility) {
      BcpsBranchObject * cb =
	curr_object->createBranchObject(model, preferredDir);
      // create a branch object for this
      bobjects.push_back(cb);
    }
  }
  // add this branch object to branchObjects_
  // compare branch objects and keep the best one at bestBranchObject_
  return 0;
}

int DcoBranchStrategyMaxInf::betterBranchObject(BcpsBranchObject * b,
					     BcpsBranchObject * bestSoFar) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}
