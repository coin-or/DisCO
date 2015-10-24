#include "DcoConicConGenerator.hpp"
#include "OsiCuts.hpp"
#include "DcoModel.hpp"
#include "DcoConstraint.hpp"
#include "DcoHelp.hpp"

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
  OsiConicSolverInterface * solver = model_->solver();
  generator_->generateCuts(*solver, *cuts);
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
