#include "DcoConGenerator.hpp"

DcoConGenerator::DcoConGenerator(DcoModel * model,
                                 char const * name,
                                 DcoCutStrategy strategy,
                                 int frequency):
  name_(name), model_(model), strategy_(strategy), frequency_(frequency) {
}

DcoConGenerator::DcoConGenerator(DcoConGenerator const & other):
  name_(other.name()), model_(other.model()),
  strategy_(other.strategy()), frequency_(other.frequency()) {
}

DcoConGenerator & DcoConGenerator::operator=(DcoConGenerator const & rhs) {
  name_ = rhs.name();
  model_ = rhs.model();
  strategy_ = rhs.strategy();
  frequency_ = rhs.frequency();
  numConsGenerated_ = rhs.numConsGenerated();
  numConsUsed_ = rhs.numConsUsed();
  time_ = rhs.time();
  calls_ = rhs.calls();
  noConsCalls_ = rhs.noConsCalls();
  return *this;
}

DcoConGenerator::~DcoConGenerator() {
  model_ = NULL;
}
