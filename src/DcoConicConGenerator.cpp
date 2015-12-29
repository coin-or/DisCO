#include "DcoConicConGenerator.hpp"
#include "OsiCuts.hpp"
#include "DcoModel.hpp"
#include "DcoConstraint.hpp"
#include "DcoHelp.hpp"

#include "OsiSolverInterface.hpp"
//#include "OsiConicSolverInterface.hpp"

/** Copy constructor. */
DcoConicConGenerator::DcoConicConGenerator (
					    const DcoConicConGenerator & other): DcoConGeneratorBase(other) {
  setGenerator(other.generator());
}

/** Assignment operator. */
DcoConicConGenerator &
  DcoConicConGenerator::operator=(
    const DcoConicConGenerator& rhs) {
  if (this != &rhs) {
    DcoConGeneratorBase::operator=(rhs);
    setGenerator(rhs.generator());
  }
  return *this;
}

// returns true if reoptimization is needed.
bool DcoConicConGenerator::generateConstraints(
  BcpsConstraintPool & conPool) {
  OsiCuts * cuts = new OsiCuts();
  OsiSolverInterface * solver = model_->solver();
  // get conic constraint information
  // todo(aykut) things may break in case of conic cuts in the root node.
  // since we use the core cones here.
  int num_cones = model_->getNumCoreCones();
  int * const * const members = model_->getConeMembers();
  OsiLorentzConeType const * types = model_->getConeTypes();
  int const * sizes = model_->getConeSizes();
  // convert cone type to osi types
  generator_->generateCuts(*solver, *cuts, num_cones, types,
			   sizes, members, 3);
  // add cuts to the constraint pool
  int num_cuts = cuts->sizeRowCuts();
  for (int i=0; i<num_cuts; ++i) {
    OsiRowCut & rcut = cuts->rowCut(i);
    int size = rcut.row().getNumElements();
    if (size>0) {
      DcoConstraint * blis_con = DcoOsiCutToConstraint(&rcut);
      conPool.addConstraint(blis_con);
    }
    assert(size>=0);
    // else if (size==0) {
    //   std::cerr << "Empty cut" << std::endl;
    // }
    // else {
    //   std::cerr << "Negative size cut" << std::endl;
    // }
  }
  if (num_cuts) {
    return true;
  }
  else {
    return false;
  }
}
