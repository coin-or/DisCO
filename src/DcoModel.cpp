#include "DcoModel.hpp"

DcoModel::DcoModel() {
}

DcoModel::~DcoModel() {
}

#if defined(__OA__)
void DcoModel::setSolver(OsiSolverInterface * solver) {
#else
void DcoModel::setSolver(OsiConicSolverInterface * solver) {
#endif
  solver_ = solver;
}
