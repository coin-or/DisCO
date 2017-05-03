#include "DcoConGenerator.hpp"

void DcoConGeneratorStats::reset() {
  numConsGenerated_ = 0;
  numConsUsed_ = 0;
  time_ = 0.0;
  numCalls_ = 0;
  numNoConsCalls_ = 0;
}

DcoConGenerator::DcoConGenerator(DcoModel * model,
                                 DcoConstraintType type,
                                 char const * name,
                                 DcoCutStrategy strategy,
                                 int frequency):
  model_(model), type_(type), name_(name), strategy_(strategy),
  frequency_(frequency) {
  stats_.reset();
}

DcoConGenerator::~DcoConGenerator() {
  model_ = NULL;
}
