#include "DcoConicConGenerator.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoConicConstraint.hpp"

#include <CglConicCutGenerator.hpp>

/// Useful constructor.
DcoConicConGenerator::DcoConicConGenerator(DcoModel * model,
                        CglConicCutGenerator * generator,
                        char const * name,
                        DcoCutStrategy strategy,
                        int frequency):
  DcoConGenerator(model, name, strategy, frequency) {
  generator_ = generator;
}

/// Copy constructor.
DcoConicConGenerator::DcoConicConGenerator(DcoConicConGenerator const & other)
  : DcoConGenerator(other) {
  generator_ = other.generator()->clone();
}

/// Destructor.
DcoConicConGenerator::~DcoConicConGenerator() {
  delete generator_;
}

/// Copy assignment operator.
DcoConicConGenerator &
DcoConicConGenerator::operator=(DcoConicConGenerator const & rhs) {
  generator_ = rhs.generator()->clone();
  return *this;
}

/// Generate constraints and add them to the pool.
bool DcoConicConGenerator::generateConstraints(BcpsConstraintPool & conPool) {
  DcoModel * model = DcoConGenerator::model();
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  CoinMessages * messages = model->dcoMessages_;
  // generated cuts will be stored in this
  OsiCuts * cuts = new OsiCuts();
  // cut generator needs solver interface, get it.
  OsiSolverInterface const * solver = model->solver();
  // get conic constraint information
  std::vector<BcpsConstraint*> & rows = model->getConstraints();
  int num_cones = model->numRelaxedRows();

  // cone members, sizes and types
  int ** members = new int*[num_cones];
  int * sizes = new int[num_cones];
  OsiLorentzConeType * types = new OsiLorentzConeType[num_cones];

  // iterate over conic constraints and collect cone information
  for (int i=0; i<num_cones; ++i) {
    DcoConicConstraint * curr = dynamic_cast<DcoConicConstraint*>
      (rows[model->relaxedRows()[i]]);
    sizes[i] = curr->getSize();
    members[i] = new int[sizes[i]];
    std::copy(curr->getMembers(), curr->getMembers()+sizes[i], members[i]);
    DcoLorentzConeType tt = curr->getType();
    if (tt==DcoLorentzCone) {
      types[i] = OSI_QUAD;
    }
    else if (tt==DcoRotatedLorentzCone) {
      types[i] = OSI_RQUAD;
    }
    else {
      message_handler->message(DISCO_UNKNOWN_CONETYPE, *messages)
        << __FILE__ << __LINE__ << CoinMessageEol;
    }
  }
  // call cut generator
  generator_->generateCuts(*solver, *cuts, num_cones, types,
                           sizes, members, 1);

  // debug message
  std::stringstream debug_msg;
  debug_msg << "Conic cut generator generated ";
  debug_msg << cuts->sizeRowCuts();
  debug_msg << " cuts.";
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                              'G', DISCO_DLOG_CUT)
    << CoinMessageEol;
  // end of debug

  // add cuts to the constraint pool
  int num_cuts = cuts->sizeRowCuts();
  for (int i=0; i<num_cuts; ++i) {
    OsiRowCut & rcut = cuts->rowCut(i);
    int num_elem = rcut.row().getNumElements();
    int const * ind = rcut.row().getIndices();
    double const * val = rcut.row().getElements();
    DcoConstraint * con =
      new DcoLinearConstraint(num_elem, ind, val, rcut.lb(), rcut.ub());
      conPool.addConstraint(con);
  }
  delete cuts;
  for (int i=0; i<num_cones; ++i) {
    delete[] members[i];
  }
  delete[] members;
  delete[] sizes;
  delete[] types;
  if (num_cuts) {
    return true;
  }
  else {
    return false;
  }
}
