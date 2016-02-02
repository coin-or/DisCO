#include "DcoSolution.hpp"

DcoSolution::DcoSolution() {
}

DcoSolution::DcoSolution(int size, double const * values, double q) {
}

DcoSolution::~DcoSolution() {
}

BcpsSolution * DcoSolution::selectNonzeros(const double etol) const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}

BcpsSolution * DcoSolution::selectFractional(const double etol) const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}
