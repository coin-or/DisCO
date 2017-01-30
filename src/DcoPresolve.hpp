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

  Usage is similar to OsiPresolve.

  DcoPresolve pinfo;
  DcoSolverInterface * presolvedModel;
  // Return an OsiSolverInterface loaded with the presolved problem.
  presolvedModel = pinfo.presolvedModel(*origModel,1.0e-8,false,numberPasses) ;
  // solve problem using presolvedModel, use it in branch and bound etc.
  presolvedModel->initialSolve() ;
  // Restate the solution and load it back into origModel.
  pinfo.postsolve(true) ;
  delete presolvedModel ;

  DcoPresolve migth eliminate or assing variables a different
  index. For a presolved model, originalModel()->originalColumns()[i] gives the
  original index of ith varaible of presolved model.

  In DcoModel::preprocess function, following fields are needed to get updated
  after presolve, (we might move this block to DcoPresolve)

  <ul>
    <li> numCols_, colLB_, colUB_,
    <li> objSense_, objCoef_
    <li> numIntegerCols_, integerCols_, isInteger_
    <li> numConicRows_, coneStart_, coneType_, coneMembers_
    <li> numLinearRows_, numRows_, rowLB_, rowUB_
    <li> matrix_
  </ul>

  In serial, ALPS calls readInstance(), then preproces(), then
  setupSelf(). In parallel, ALPS master calls readInstance(), then
  preprocess(), then setupSelf() and then broadcast model to other
  processes. Other processes receive problem data, call decodeToSelf() and then
  setupSelf().

  # Questions

  What happens to eliminated columns? Are tehy fixed to some value? zero?

  What happens to eliminated rows? Are they just redundant rows, or they
  are fixed to some value?

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
  DcoPresolve();
  virtual ~DcoPresolve();
  //@}

#if defined(__OA__)
  virtual OsiSolverInterface * presolvedModel(OsiSolverInterface & origModel,
                                              double feasibilityTolerance=0.0,
                                              bool keepIntegers=true,
                                              int numberPasses=5,
                                              const char * prohibited=NULL,
                                              bool doStatus=true,
                                              const char * rowProhibited=NULL);
  /*! \brief Return a pointer to the presolved model. */
  OsiSolverInterface * model() const;
  /// Return a pointer to the original model
  OsiSolverInterface * originalModel() const;
#else
  virtual OsiConicSolverInterface *
    presolvedModel(OsiConicSolverInterface & origModel,
                   double feasibilityTolerance=0.0,
                   bool keepIntegers=true,
                   int numberPasses=5,
                   const char * prohibited=NULL,
                   bool doStatus=true,
                   const char * rowProhibited=NULL);
  /*! \brief Return a pointer to the presolved model. */
  OsiConicSolverInterface * model() const;
  /// Return a pointer to the original model
  OsiConicSolverInterface * originalModel() const;
#endif
  virtual void postsolve(bool updateStatus=true);
  bool improve_bounds(DcoModel * model);

private:
  /// Disable copy constructor.
  DcoPresolve(DcoPresolve const & other);
  /// Disable copy assignment operator.
  DcoPresolve & operator=(DcoPresolve const & rhs);
};

#endif
