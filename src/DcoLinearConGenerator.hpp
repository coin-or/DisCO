#ifndef DcoLinearConGenerator_hpp_
#define DcoLinearConGenerator_hpp_

#include "DcoConGenerator.hpp"

class CglCutGenerator;

/*!
   DcoLinearConGenerator implements constraint generator interface for linear
   cut generating procedures.
*/

class DcoLinearConGenerator: virtual public DcoConGenerator {
  /// The CglCutGenerator object.
  CglCutGenerator * generator_;
public:
  ///@name Constructors and Destructor
  //@{
  /// Useful constructor.
  DcoLinearConGenerator(DcoModel * model,
                        CglCutGenerator * generator,
                        char const * name = NULL,
                        DcoCutStrategy strategy = DcoCutStrategyAuto,
                        int frequency = 1);
  /// Destructor.
  virtual ~DcoLinearConGenerator();
  //@}

  ///@name Constraint generator functions
  //@{
  /// Generate constraints and add them to the pool.  return true if resolve is
  /// needed (because the state of the solver interface has been modified).
  virtual bool generateConstraints(BcpsConstraintPool & conPool);
  //@}

  // Get cut generator.
  CglCutGenerator const * generator() const { return generator_; }
private:
  /// Default constructor is not allowed.
  DcoLinearConGenerator();
  /// Disable copy constructor.
  DcoLinearConGenerator(DcoLinearConGenerator const & other);
  /// Disable copy assignment operator.
  DcoLinearConGenerator & operator=(DcoLinearConGenerator const & rhs);
};

#endif
