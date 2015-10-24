//#############################################################################
// This file is modified from COIN/Cbc/CbcCutGenerator.hpp
//#############################################################################


#ifndef DcoConicConGenerator_h_
#define DcoConicConGenerator_h_

#include "OsiConicSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "DcoConGeneratorBase.hpp"
#include "CglConicCutGenerator.hpp"
#include "BcpsObjectPool.h"

class DcoModel;

class OsiRowCut;
class OsiRowCutDebugger;
class CglCutGenerator;


//#############################################################################

class DcoConicConGenerator: virtual public DcoConGeneratorBase  {
protected:
  /** The CglCutGenerator object. */
  CglConicCutGenerator * generator_;
public:
  /**@name Constructors and destructors */
  //@{
  /** Default constructor. */
  DcoConicConGenerator() : DcoConGeneratorBase(),
                      generator_(NULL) {}
  /** Useful constructor. */
  DcoConicConGenerator(DcoModel * model,
		  CglConicCutGenerator * generator,
		  const char * name = NULL,
		  DcoCutStrategy strategy = DcoCutStrategyAuto,
		  int cutGenerationFrequency = 1,
		  bool normal = true,
		  bool atSolution = false,
		  bool infeasible = false)
    : DcoConGeneratorBase(model, name, strategy, cutGenerationFrequency,
                          normal, atSolution, infeasible) {
    generator_ = generator;
  }

  /** Copy constructor. */
  DcoConicConGenerator (const DcoConicConGenerator &);

  /** Assignment operator. */
  DcoConicConGenerator & operator=(const DcoConicConGenerator& rhs);

  /** Destructor. */
  virtual ~DcoConicConGenerator()
  {
    if (generator_) {
      delete generator_;
      generator_ = NULL;
    }
  }
  //@}
  //@{
  /** \name Get methods */
  /** Get conic cut generator
  */
  CglConicCutGenerator * generator() const { return generator_; }
  //@}

  //@{
  /** \name Set methods */
  /** Set conic cut generator
  */
  void * setGenerator(CglConicCutGenerator * gen) { generator_ = gen; }
  //@}

  /** \name Generate Constraints */
  //@{
  /** Generate cons for the client model.

      Evaluate the state of the client model and decide whether to
      generate cons. The generated cons are inserted into and returned
      in the collection of cons \p cs.

      The routine returns true if reoptimisation is needed (because the
      state of the solver interface has been modified).
  */
  virtual bool generateConstraints(BcpsConstraintPool &conPool);
  //@}

};

#endif
