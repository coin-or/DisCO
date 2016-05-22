#include "DcoPresolve.hpp"
#include "DcoModel.hpp"
#include "DcoConicConstraint.hpp"

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
    DcoLorentzConeType type = curr->getType();
    if (type==DcoLorentzCone) {
      int lead_index = curr->getMembers()[0];
      double lead_bound = colub[lead_index];
      // set lower and upper bounds of the other cone members
      for (int j=1; j<curr->getSize(); j++) {
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
      model->dcoMessageHandler_->message(DISCO_NOT_IMPLEMENTED,
                                         *model->dcoMessages_)
        << __FILE__ << __LINE__ << CoinMessageEol;
    }
    else {
      model->dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE,
                                         *model->dcoMessages_)
        << type << CoinMessageEol;
    }
  }
  return updated;
}
