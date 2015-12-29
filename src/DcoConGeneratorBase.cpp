#include "DcoConGeneratorBase.hpp"

/** Useful constructor. */
DcoConGeneratorBase::DcoConGeneratorBase(DcoModel * model,
		      const char * name,
		      DcoCutStrategy strategy ,
		      int cutGenerationFrequency,
		      bool normal,
		      bool atSolution,
					 bool infeasiblex) {
  model_ = model;
  name_ = name;
  strategy_ = strategy;
  cutGenerationFrequency_ = cutGenerationFrequency;
  normal_ = normal;
  atSolution_ = atSolution;
  whenInfeasible_ = infeasiblex;
  numConsGenerated_ = 0;
  numConsUsed_ = 0;
  time_ = 0.0;
  calls_ = 0;
  noConsCalls_ = 0;
}

/** Copy constructor. */
DcoConGeneratorBase::DcoConGeneratorBase (const DcoConGeneratorBase & other) {
  model_ = other.model_;
  setName(other.name().c_str());
  setStrategy(other.strategy());
  setCutGenerationFreq(other.cutGenerationFreq());
  setNormal(other.normal());
  setAtSolution(other.atSolution());
  setWhenInfeasible(other.whenInfeasible());
  numConsGenerated_ = other.numConsGenerated();
  numConsUsed_ = other.numConsUsed();
  time_ = other.time();
  calls_ = other.calls();
  noConsCalls_ = other.noConsCalls();
}

/** Assignment operator. */
DcoConGeneratorBase & DcoConGeneratorBase::operator=(
		      const DcoConGeneratorBase & rhs) {
  if (this!=&rhs) {
    model_ = rhs.model_;
    setName(rhs.name().c_str());
    setStrategy(rhs.strategy());
    setCutGenerationFreq(rhs.cutGenerationFreq());
    setNormal(rhs.normal());
    setAtSolution(rhs.atSolution());
    setWhenInfeasible(rhs.whenInfeasible());
    numConsGenerated_ = rhs.numConsGenerated();
    numConsUsed_ = rhs.numConsUsed();
    time_ = rhs.time();
    calls_ = rhs.calls();
    noConsCalls_ = rhs.noConsCalls();
  }
  return *this;
}
