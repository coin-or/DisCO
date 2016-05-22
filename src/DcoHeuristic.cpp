#include "DcoHeuristic.hpp"

void DcoHeurStats::reset() {
  numCalls_ = 0;
  numNoSolCalls_ = 0;
  time_ = 0.0;
  numSolutions_ = 0;
}

DcoHeuristic::DcoHeuristic(DcoModel * model, char const * name,
                           DcoHeurStrategy strategy,
                           int frequency) {
  model_ = model;
  name_ = name;
  strategy_ = strategy;
  frequency_ = frequency;
  type_ = DcoHeurTypeNotSet;
  stats_.reset();
}
