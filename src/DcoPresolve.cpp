#include "DcoPresolve.hpp"
#include "DcoModel.hpp"
#include "DcoConicConstraint.hpp"


DcoPresolve::DcoPresolve():
  OsiPresolve() {
}

DcoPresolve::~DcoPresolve() {
}

#if defined(__OA__)
OsiSolverInterface *
  DcoPresolve::presolvedModel(OsiSolverInterface & origModel,
                              double feasibilityTolerance,
                              bool keepIntegers,
                              int numberPasses,
                              const char * prohibited,
                              bool doStatus,
                              const char * rowProhibited) {
  return OsiPresolve::presolvedModel(origModel,
                                     feasibilityTolerance,
                                     keepIntegers,
                                     numberPasses,
                                     prohibited,
                                     doStatus,
                                     rowProhibited);
}

  /*! \brief Return a pointer to the presolved model. */
OsiSolverInterface * DcoPresolve::model() const {
  return OsiPresolve::model();
}

/// Return a pointer to the original model
OsiSolverInterface * DcoPresolve::originalModel() const {
  return OsiPresolve::originalModel();
}

void DcoPresolve::postsolve(bool updateStatus) {
  OsiPresolve::postsolve(updateStatus);
}

#else
OsiConicSolverInterface *
  DcoPresolve::presolvedModel(OsiConicSolverInterface & origModel,
                              double feasibilityTolerance,
                              bool keepIntegers,
                              int numberPasses,
                              const char * prohibited,
                              bool doStatus,
                              const char * rowProhibited) {
  std::cerr << "function: " << __FUNCTION__
            << "file: " << __FILE__
            << " line: " << __LINENO__
            << " Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}

/*! \brief Return a pointer to the presolved model. */
OsiConicSolverInterface * DcoPresolve::model() const {
  std::cerr << "function: " << __FUNCTION__
            << "file: " << __FILE__
            << " line: " << __LINENO__
            << " Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}

/// Return a pointer to the original model
OsiConicSolverInterface * DcoPresolve::originalModel() const {
  std::cerr << "function: " << __FUNCTION__
            << "file: " << __FILE__
            << " line: " << __LINENO__
            << " Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}

void DcoPresolve::postsolve(bool updateStatus) {
  std::cerr << "function: " << __FUNCTION__
            << "file: " << __FILE__
            << " line: " << __LINENO__
            << " Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}

#endif

bool DcoPresolve::improve_bounds(DcoModel * model) {
  bool updated = false;
  // iterate over cones and improve bounds
  int num_linear = model->getNumCoreLinearConstraints();
  int num_cones = model->getNumCoreConicConstraints();

  // get bounds, these will get updated through the process.
  double * collb = model->colLB();
  double * colub = model->colUB();

  for (int i=num_linear; i<num_linear+num_cones; ++i) {
    // get constraint
    DcoConicConstraint * curr =
      dynamic_cast<DcoConicConstraint*>(model->getConstraints()[i]);
    // get type
    DcoLorentzConeType type = curr->coneType();
    if (type==DcoLorentzCone) {
      int lead_index = curr->coneMembers()[0];
      double lead_bound = colub[lead_index];
      // set lower and upper bounds of the other cone members
      for (int j=1; j<curr->coneSize(); j++) {
        // if lower bound is less than -lead_bound, update it
        if (collb[j]<-lead_bound) {
          // debug stuff
          std::stringstream debug_msg;
          debug_msg << "Lower bound of col ";
          debug_msg << j;
          debug_msg << " is updated from ";
          debug_msg << collb[j];
          debug_msg << " to ";
          debug_msg << -lead_bound;
          model->dcoMessageHandler_->message(0, "Dco", debug_msg.str().c_str(),
                                             'G', DISCO_DLOG_PRESOLVE)
            << CoinMessageEol;
          // end of debug stuff
          // do the update
          collb[j] = -lead_bound;
          updated = true;
        }
        // if upper bound is larger than lead_bound, update it
        if (colub[j]>lead_bound) {
          // debug stuff
          std::stringstream debug_msg;
          debug_msg << "Upper bound of col ";
          debug_msg << j;
          debug_msg << " is updated from ";
          debug_msg << colub[j];
          debug_msg << " to ";
          debug_msg << lead_bound;
          model->dcoMessageHandler_->message(0, "Dco", debug_msg.str().c_str(),
                                             'G', DISCO_DLOG_PRESOLVE)
            << CoinMessageEol;
          // end of debug stuff
          // do the update
          colub[j] = lead_bound;
          updated = true;
        }
      }
    }
    else if (type==DcoRotatedLorentzCone) {
      model->dcoMessageHandler_->message(0, "Dco",
                                         "Presolve is not implemented for "
                                         "rotated cones yet, skipping...",
                                         'G', DISCO_DLOG_PRESOLVE)
        << CoinMessageEol;
    }
    else {
      model->dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE,
                                         *model->dcoMessages_)
        << type << CoinMessageEol;
    }
  }
  return updated;
}
