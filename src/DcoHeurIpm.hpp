#ifndef DcoHeurIpm_hpp_
#define DcoHeurIpm_hpp_

// Disco headers
#include "DcoHeuristic.hpp"

/*!
  Calls IPM solver when current relaxation solution is integer feasible but
  conic infeasible.

  Fixes integer variables and solves SOCO to find solutions.
*/

/** Heuristic base class */
class DcoHeurIpm: virtual public DcoHeuristic {
public:
  ///@name Constructors and Destructor.
  //@{
  /// Useful constructor.
  DcoHeurIpm(DcoModel * model, char const * name,
                  DcoHeurStrategy strategy, int frequency);
  /// Destructor.
  virtual ~DcoHeurIpm() { }
  //@}

  ///@name Finding solutions.
  //@{
  /// returns a solution if found, NULL otherwise.
  virtual DcoSolution * searchSolution();
  //@}

private:
  /// Disable default constructor.
  DcoHeurIpm();
  /// Disable copy constructor.
  DcoHeurIpm(const DcoHeurIpm & other);
  /// Disable copy assignment operator
  DcoHeurIpm & operator=(const DcoHeurIpm & rhs);
};

#endif
