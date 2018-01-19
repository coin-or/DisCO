/*===========================================================================*
 * This file is part of the Discrete Conic Optimization (DisCO) Solver.      *
 *                                                                           *
 * DisCO is distributed under the Eclipse Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *          Aykut Bulut, Lehigh University                                   *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and      *
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#ifndef DcoPresolve_hpp_
#define DcoPresolve_hpp_

#include "OsiPresolve.hpp"

#if defined(__OA__)
  #include <OsiClpSolverInterface.hpp>
#else
  #include <OsiConicSolverInterface.hpp>
#endif


class DcoModel;

/*!
  This class is for pre-processing of the conic problems. It is built
  on top of OsiPresolve class which has pre-procesiing methods for linear
  problems.

  This class accomodates the two cases when the underlying solver interface is
  a linear (OsiSolverInterface) or conic one (OsiConicSolverInterface).

 */

class DcoPresolve: virtual public OsiPresolve {
#if defined(__OA__)
  OsiSolverInterface * origModel_;
  OsiSolverInterface * presolvedModel_;
#else
  OsiConicSolverInterface * origModel_;
  OsiConicSolverInterface * presolvedModel_;
#endif
public:
  ///@name Constructors and Destructor
  //@{
#if defined(__OA__)
  DcoPresolve(OsiSolverInterface * origModel);
#else
  DcoPresolve(OsiConicSolverInterface * origModel);
#endif
  virtual ~DcoPresolve();
  //@}

  ///@name Presolve Functions
  //@{
  /// Preprocess the given model and store pointer to the pre-processed
  /// copy in presolvedModel_;
  void presolve();
#if defined(__OA__)
  /// Get preprocessed problem.
  OsiSolverInterface * presolvedModel();
#else
  /// Get preprocessed problem.
  OsiConicSolverInterface * presolvedModel();
#endif
  //@}

  bool improve_bounds(DcoModel * model);

  ///@name Postsolve Function
  //@{
  virtual void postsolve(bool updateStatus=true);
  //@}
private:
  /// Disable copy constructor.
  DcoPresolve(DcoPresolve const & other);
  /// Disable copy assignment operator.
  DcoPresolve & operator=(DcoPresolve const & rhs);
};

#endif
