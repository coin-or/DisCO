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
                       DcoConstraintType type,
                       char const * name = NULL,
                       DcoCutStrategy strategy = DcoCutStrategyAuto,
                       int frequency = 1);
  /// Destructor.
  virtual ~DcoConicConGenerator();
  //@}


  ///@name Constraint generator functions
  //@{
  /// Generate constraints and add them to the pool.
  bool generateConstraints(BcpsConstraintPool & conPool);
  //@}

  // Get cut generator.
  CglConicCutGenerator const * generator() const { return generator_; }
private:
  /// Disable default constructor
  DcoConicConGenerator();
  /// Disable copy constructor.
  DcoConicConGenerator(DcoConicConGenerator const & other);
  /// Disable copy assignment operator.
  DcoConicConGenerator & operator=(DcoConicConGenerator const & rhs);
};

#endif
