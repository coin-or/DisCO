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
// This file is modified from COIN/Cbc/CbcCutGenerator.hpp
//#############################################################################


#ifndef DcoConGenerator_h_
#define DcoConGenerator_h_

#include "OsiConicSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "DcoConGeneratorBase.hpp"

class DcoModel;

class OsiRowCut;
class OsiRowCutDebugger;
class CglCutGenerator;


//#############################################################################

/** Interface between Dco and Cut Generation Library.
    \c DcoConGenerator is intended to provide an intelligent
    interface between Dco and the cutting plane algorithms in the CGL.
    A \c DcoConGenerator is bound to a \c CglCutGenerator and
    to an \c DcoModel. It contains parameters which control when and
    how the \c generateCuts method of the \c CglCutGenerator will be called.

    The builtin decision criteria available to use when deciding whether to
    generate cons are: at root, autmatic, every <i>X</i> nodes, when a
    solution is found, and when a subproblem is found to be infeasible.
*/

class DcoConGenerator: virtual public DcoConGeneratorBase {
protected:
  /** The CglCutGenerator object. */
  CglCutGenerator * generator_;
public:
  /**@name Constructors and destructors */
  //@{
  /** Default constructor. */
  DcoConGenerator() : DcoConGeneratorBase(),
                      generator_(NULL) {}
  /** Useful constructor. */
  DcoConGenerator(DcoModel * model,
		  CglCutGenerator * generator,
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
  DcoConGenerator (const DcoConGenerator &);

  /** Assignment operator. */
  DcoConGenerator & operator=(const DcoConGenerator& rhs);

  /** Destructor. */
  virtual ~DcoConGenerator()
  {
    if (generator_) {
      delete generator_;
      generator_ = NULL;
    }
  }
  //@}
  /** \name Get methods */
  //@{
  //** Get cut generator
  CglCutGenerator * generator() const { return generator_; }
  //@}
  /** \name Set methods */
  //@{
  //** Set cut generator
  void setGenerator(CglCutGenerator * gen) { generator_ = gen; }
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
