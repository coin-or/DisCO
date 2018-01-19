/*===========================================================================*
 * This file is part of the Discrete Conic Optimization (DisCO) Solver.      *
 *                                                                           *
 * DisCO is distributed under the Eclipse Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *          Aykut Bulut, Lehigh University                                   *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and      *
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#include "DcoConicConGenerator.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoConicConstraint.hpp"

#include <CglConicCutGenerator.hpp>


extern std::vector<char const *> const dcoConstraintTypeName;

/// Useful constructor.
DcoConicConGenerator::DcoConicConGenerator(DcoModel * model,
                        CglConicCutGenerator * generator,
                        DcoConstraintType type,
                        char const * name,
                        DcoCutStrategy strategy,
                        int frequency):
  DcoConGenerator(model, type, name, strategy, frequency) {
  generator_ = generator;
}

/// Destructor.
DcoConicConGenerator::~DcoConicConGenerator() {
  delete generator_;
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
    sizes[i] = curr->coneSize();
    members[i] = new int[sizes[i]];
    std::copy(curr->coneMembers(), curr->coneMembers()+sizes[i], members[i]);
    DcoLorentzConeType tt = curr->coneType();
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
  message_handler->message(DISCO_CUT_GENERATED, *messages)
    << model->broker()->getProcRank()
    // todo(aykut) fix name
    << dcoConstraintTypeName[type()]
    << cuts->sizeRowCuts()
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
    con->setConstraintType(type());
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
