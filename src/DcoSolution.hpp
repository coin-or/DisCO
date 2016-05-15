#ifndef DcoSolution_hpp_
#define DcoSolution_hpp_

#include <BcpsSolution.h>

/*!

  DcoSolution represents a solution for the MISOCO problem (master problem)
  represented by the DcoModel. Inherits BcpsSolution.

  # Design concerns:

  BcpsSolution can also keep objects (variables) associated with the solution.
  We are not using that for now, but it can be usefull.

 */

class DcoSolution: public BcpsSolution {
public:
  DcoSolution();
  DcoSolution(int size, double const * values, double quality);
  virtual ~DcoSolution();
  virtual BcpsSolution * selectNonzeros(const double etol=1e-5) const;
  virtual BcpsSolution * selectFractional(const double etol=1e-5) const;
private:
  DcoSolution(DcoSolution const & other);
  DcoSolution & operator=(DcoSolution const & rhs);
};

#endif
