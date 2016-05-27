#include "DcoConGenerator.hpp"

void DcoConGeneratorStats::reset() {
  numConsGenerated_ = 0;
  numConsUsed_ = 0;
  time_ = 0.0;
  numCalls_ = 0;
  numNoConsCalls_ = 0;
}

DcoConGenerator::DcoConGenerator(DcoModel * model,
                                 char const * name,
                                 DcoCutStrategy strategy,
                                 int frequency):
  name_(name), model_(model), strategy_(strategy), frequency_(frequency) {
  stats_.reset();
}

DcoConGenerator::~DcoConGenerator() {
  model_ = NULL;
}
