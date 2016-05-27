/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           *
 *                                                                           *
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


//#############################################################################
// This file is modified from COIN/Cbc/CbcHeuristic.hpp
//#############################################################################


#ifndef DcoHeuristic_hpp_
#define DcoHeuristic_hpp_

// STL headers
#include <string>

// Disco headers
#include "Dco.hpp"

class DcoModel;
class DcoSolution;

//#############################################################################

/*!
  This is a class for keeping heuristic statistics.
*/
class DcoHeurStats {
  /// number of calls to the heuristic
  int numCalls_;
  /// number of calls where no solution found.
  int numNoSolCalls_;
  // total CPU time spent.
  // todo(aykut) C++11 has chrono for this, but we are not using C++11 :(
  // but we should soon.
  double time_;
  // number of solutions found by the heuristic
  int numSolutions_;
public:
  ///@name Constructors and Destructors
  //@{
  /// Default constructor
  DcoHeurStats() { reset(); }
  // note(aykut) I do not define copy constructor and copy assignment operator
  // default ones will be fine since the class is simple.
  //@}

  ///@name Update statistics
  //@{
  void addCalls(int c=1) { numCalls_ += c; }
  void addNoSolCalls(int n=1) { numNoSolCalls_ += n; }
  void addTime(double t) { time_ += t; }
  void addNumSolutions(int n=1) { numSolutions_ += n; }
  /// Reset statistics to 0.
  void reset();
  //@}

  ///@name Get statistics
  //@{
  int numCalls() const { return numCalls_; }
  int numNoSolCalls() const { return numNoSolCalls_; }
  double time() const { return time_; }
  int numSolutions() const { return numSolutions_; }
  //@}
};

/*!
  This is an abstract base class for heuristics.
*/

/** Heuristic base class */
class DcoHeuristic {
  /// Pointer to the model.
  DcoModel * model_;
  /// heuristic type.
  DcoHeurType type_;
  /// Name of heuristic.
  std::string name_;
  /** When to call findSolution() routine.
      DcoHeurStrategyNone:     disable
      DcoHeurStrategyRoot:     just root
      DcoHeurStrategyAuto:     automatically decided by DISCO
      DcoHeurStrategyPeriodic: every 't' nodes
      DcoHeurStrategyBeforeRoot: before solving first LP
  */
  DcoHeurStrategy strategy_;
  /// The frequency with which to call the heuristic */
  // todo(aykut) isn't this a part of strategy?
  int frequency_;
  /// Statistics.
  DcoHeurStats stats_;
public:
  ///@name Constructors and Destructor.
  //@{
  /// Useful constructor.
  DcoHeuristic(DcoModel * model, char const * name,
               DcoHeurStrategy strategy, int frequency);
  /// Destructor.
  virtual ~DcoHeuristic() { }
  //@}

  ///@name Get methods.
  //@{
  DcoModel * model() const { return model_; }
  DcoHeurType type() const { return type_; }
  std::string & name() { return name_; }
  DcoHeurStrategy strategy() const { return strategy_; }
  int frequency() const { return frequency_; }
  DcoHeurStats & stats() { return stats_; }
  //@}

  ///@name Set methods
  //@{
  void setType(DcoHeurType type) { type_ = type; }
  //@}
  ///@name Finding solutions.
  //@{
  /// returns a solution if found, NULL otherwise.
  virtual DcoSolution * searchSolution() = 0;
  //@}

private:
  /// Disable default constructor.
  DcoHeuristic();
  /// Disable copy constructor.
  DcoHeuristic(const DcoHeuristic & other);
  /// Disable copy assignment operator
  DcoHeuristic & operator=(const DcoHeuristic & rhs);
};

#endif
