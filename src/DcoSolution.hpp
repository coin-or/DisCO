#ifndef DcoSolution_hpp_
#define DcoSolution_hpp_

#include <BcpsSolution.h>

class DcoSolution: public BcpsSolution {
public:
  DcoSolution();
  DcoSolution(int size, double const * values, double q);
  virtual ~DcoSolution();
  virtual BcpsSolution * selectNonzeros(const double etol=1e-5) const;
  virtual BcpsSolution * selectFractional(const double etol=1e-5) const;
private:
  DcoSolution(DcoSolution const & other);
  DcoSolution & operator=(DcoSolution const & rhs);
};

#endif
