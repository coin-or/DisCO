#include "DcoLinearConGenerator.hpp"
#include "DcoModel.hpp"
#include "DcoMessage.hpp"

#include <CglCutGenerator.hpp>

/// Useful constructor.
DcoLinearConGenerator::DcoLinearConGenerator(DcoModel * model,
                        CglCutGenerator * generator,
                        char const * name ,
                        DcoCutStrategy strategy,
                        int frequency):
  DcoConGenerator(model, name, strategy, frequency) {
  generator_ = generator;
}

/// Copy constructor.
DcoLinearConGenerator::DcoLinearConGenerator(DcoLinearConGenerator const & other)
  : DcoConGenerator(other) {
  generator_ = other.generator()->clone();
}

/// Destructor.
DcoLinearConGenerator::~DcoLinearConGenerator() {
  delete generator_;
}

/// Generate constraints and add them to the pool.
bool DcoLinearConGenerator::generateConstraints(BcpsConstraintPool & conPool) {
  CoinMessageHandler * message_handler = model()->dcoMessageHandler_;
  CoinMessages * messages = model()->dcoMessages_;
  // message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
  //   << __FILE__ << __LINE__ << CoinMessageEol;
  message_handler->message(0, "Dco", "Linear cut generator called"
                           ", we do nothing for now.", 'G', DISCO_DLOG_BRANCH)
    << CoinMessageEol;
  return false;
}

/// Copy assignment operator.
DcoLinearConGenerator &
DcoLinearConGenerator::operator=(DcoLinearConGenerator const & rhs) {
  generator_ = rhs.generator()->clone();
  return *this;
}
