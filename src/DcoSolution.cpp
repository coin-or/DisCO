#include "DcoSolution.hpp"

DcoSolution::DcoSolution() {
}

DcoSolution::DcoSolution(int size, double const * values, double quality):
  BcpsSolution(size, values, quality) {
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

/// Encodes the solution into AlpsEncoded object and return pointer to it.
AlpsReturnStatus DcoSolution::encode(AlpsEncoded * encoded) const {
  AlpsReturnStatus status;
  status = AlpsSolution::encode(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected encode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  // encode Bcps fields
  status = BcpsSolution::encode(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected encode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  // Nothing to do for DisCO part.
  return status;
}

/// Decodes the given object into a new solution and returns the pointer to
/// it.
AlpsKnowledge * DcoSolution::decode(AlpsEncoded & encoded) const {
  AlpsReturnStatus status;
  DcoSolution * sol = new DcoSolution();
  status = sol->AlpsSolution::decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected decode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  status = sol->BcpsSolution::decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected decode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  return sol;
}

/// Decode the given AlpsEncoded object into this.
AlpsReturnStatus DcoSolution::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status;
  status = AlpsSolution::decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected decode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  status = BcpsSolution::decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected decode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  return status;
}
