#ifndef DcoConicConGenerator_hpp_
#define DcoConicConGenerator_hpp_

#include "DcoConGenerator.hpp"

class CglConicCutGenerator;

class DcoConicConGenerator: virtual public DcoConGenerator {
  /// The CglCutGenerator object.
  CglConicCutGenerator * generator_;
public:
  ///@name Constructors and Destructor
  //@{
  /// Useful constructor.
  DcoConicConGenerator(DcoModel * model,
                       CglConicCutGenerator * generator,
                       char const * name = NULL,
                       DcoCutStrategy strategy = DcoCutStrategyAuto,
                       int frequency = 1);
  /// Copy constructor.
  DcoConicConGenerator(DcoConicConGenerator const & other);
  /// Destructor.
  virtual ~DcoConicConGenerator();
  //@}

  /// Copy assignment operator.
  DcoConicConGenerator & operator=(DcoConicConGenerator const & rhs);

  ///@name Constraint generator functions
  //@{
  /// Generate constraints and add them to the pool.
  bool generateConstraints(BcpsConstraintPool & conPool);
  //@}

  // Get cut generator.
  CglConicCutGenerator const * generator() const { return generator_; }
private:
  DcoConicConGenerator();
};

#endif
