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


// CoinUtils
#include <CoinMpsIO.hpp>

// Disco headers
#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoTreeNode.hpp"
#include "DcoNodeDesc.hpp"
#include "DcoVariable.hpp"
#include "DcoLinearConstraint.hpp"
#include "DcoConicConstraint.hpp"
#include "DcoBranchStrategyMaxInf.hpp"
#include "DcoBranchStrategyPseudo.hpp"
#include "DcoBranchStrategyStrong.hpp"
#include "DcoConGenerator.hpp"
#include "DcoLinearConGenerator.hpp"
#include "DcoConicConGenerator.hpp"
#include "DcoSolution.hpp"
#include "DcoPresolve.hpp"
#include "DcoHeuristic.hpp"
#include "DcoHeurRounding.hpp"
#include "DcoCbfIO.hpp"

// MILP cuts
#include <CglCutGenerator.hpp>
#include <CglProbing.hpp>
#include <CglClique.hpp>
#include <CglOddHole.hpp>
#include <CglFlowCover.hpp>
#include <CglKnapsackCover.hpp>
#include <CglMixedIntegerRounding2.hpp>
#include <CglGomory.hpp>
#include <CglGMI.hpp>
#include <CglTwomir.hpp>

// Conic cuts
#include <CglConicCutGenerator.hpp>
#include <CglConicIPM.hpp>
#include <CglConicIPMint.hpp>
#include <CglConicOA.hpp>

// STL headers
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <cmath>
#include <iomanip>

// ordering of conNames should match the ordering of DcoConstraintType enum
// type.
char const * conNames[] = {
  "NotSet",
  "Core",
  "Clique",
  "FCover",
  "Gomory",
  //"GMI",
  "Knap",
  // Linear MIR
  "MIR",
  "OddHole",
  "Probe",
  "TwoMIR",
  // Conic cuts
  "IPM",
  "IPMint",
  "OA",
  // Conic MIR
  "CMIR",
  "GD1"
};

std::vector<char const *> const
  dcoConstraintTypeName (conNames, conNames + DcoConstraintTypeEnd);

DcoModel::DcoModel() {
  problemName_ = "";
  solver_ = NULL;
  colLB_ = NULL;
  colUB_ = NULL;
  rowLB_ = NULL;
  rowUB_ = NULL;
  numCols_ = 0;
  numRows_ = 0;
  numLinearRows_ = 0;
  numConicRows_ = 0;
  objSense_ = 0.0;
  objCoef_ = NULL;
  numIntegerCols_ = 0;
  integerCols_ = NULL;
  isInteger_ = NULL;
  matrix_ = NULL;
  coneStart_ = NULL;
  coneMembers_ = NULL;
  coneType_ = NULL;

  dcoPar_ = new DcoParams();
  numRelaxedCols_ = 0;
  relaxedCols_ = NULL;
  numRelaxedRows_ = 0;
  relaxedRows_ = NULL;
  dcoMessageHandler_ = new CoinMessageHandler();
  dcoMessages_ = new DcoMessage();
  // set branch strategy
  branchStrategy_ = NULL;
  rampUpBranchStrategy_ = NULL;
  // cut and heuristics objects will be set in setupSelf.

  initOAcuts_ = 0;

  dcoMessageHandler_->setPrefix(0);
  dcoMessageHandler_->message(DISCO_WELCOME, *dcoMessages_)
    << DISCO_VERSION
    << __DATE__
    << CoinMessageEol;
  dcoMessageHandler_->setPrefix(1);
}

DcoModel::~DcoModel() {
  // solver_ is freed in main function.
  if (colLB_) {
    delete[] colLB_;
    colLB_=NULL;
  }
  if (colUB_) {
    delete[] colUB_;
    colUB_=NULL;
  }
  if (rowLB_) {
    delete[] rowLB_;
    rowLB_=NULL;
  }
  if (rowUB_) {
    delete[] rowUB_;
    rowUB_=NULL;
  }
  if (objCoef_) {
    delete[] objCoef_;
    objCoef_=NULL;
  }
  if (integerCols_) {
    delete[] integerCols_;
    integerCols_=NULL;
  }
  if (isInteger_) {
    delete[] isInteger_;
    isInteger_=NULL;
  }
  if (matrix_) {
    delete matrix_;
    matrix_=NULL;
  }
  if (coneStart_) {
    delete[] coneStart_;
    coneStart_=NULL;
  }
  if (coneMembers_) {
    delete[] coneMembers_;
    coneMembers_=NULL;
  }
  if (coneType_) {
    delete[] coneType_;
    coneType_=NULL;
  }
  if (branchStrategy_) {
    delete branchStrategy_;
    branchStrategy_=NULL;
  }
  if (rampUpBranchStrategy_) {
    delete rampUpBranchStrategy_;
    rampUpBranchStrategy_=NULL;
  }
  if (dcoPar_) {
    delete dcoPar_;
    dcoPar_=NULL;
  }
  if (dcoMessageHandler_) {
    delete dcoMessageHandler_;
    dcoMessageHandler_=NULL;
  }
  if (dcoMessages_) {
    delete dcoMessages_;
    dcoMessages_=NULL;
  }
  if (relaxedCols_) {
    delete[] relaxedCols_;
    relaxedCols_=NULL;
  }
  if (relaxedRows_) {
    delete[] relaxedRows_;
    relaxedRows_=NULL;
  }
  std::map<DcoConstraintType, DcoConGenerator*>::iterator it;
  for (it=conGenerators_.begin();
       it!=conGenerators_.end(); ++it) {
    delete it->second;
  }
  conGenerators_.clear();
  for (std::vector<DcoHeuristic*>::iterator it=heuristics_.begin();
       it!=heuristics_.end(); ++it) {
    delete *it;
  }
  heuristics_.clear();
}

#if defined(__OA__)
void DcoModel::setSolver(OsiSolverInterface * solver) {
  solver_ = solver;
}
#else
void DcoModel::setSolver(OsiConicSolverInterface * solver) {
  solver_ = solver;
}
#endif

// reads problem from the given file and sets the fields required by setupself
// only.
// setupSelf needs dcoPar, objSense_, variables_ and constraints_
// dcoPar_ is already set up by the AlpsKnowledgeBroker::initializeSearch().
void DcoModel::readInstance(char const * dataFile) {
  // get input file name
  std::string input_file(dataFile);
  std::string base_name = input_file.substr(0, input_file.rfind('.'));
  std::string extension = input_file.substr(input_file.rfind('.')+1);
  if (!extension.compare("mps")) {
    readInstanceMps(dataFile);
  }
  else if (!extension.compare("cbf")) {
    problemName_ = base_name;
    readInstanceCbf(dataFile);
  }
  else {
    dcoMessageHandler_->message(DISCO_READ_MPSCBFFILEONLY,
                                *dcoMessages_) << CoinMessageEol;
  }

  // == log cone information messages
  if (numConicRows_) {
    dcoMessageHandler_->message(DISCO_READ_CONESTATS1,
                                *dcoMessages_) << numConicRows_
                                               << CoinMessageEol;
    for (int i=0; i<numConicRows_; ++i) {
      dcoMessageHandler_->message(DISCO_READ_CONESTATS2,
                                  *dcoMessages_)
        << i
        << coneStart_[i+1] - coneStart_[i]
        << coneType_[i]
        << CoinMessageEol;
    }
  }
  else {
    dcoMessageHandler_->message(DISCO_READ_NOCONES,
                                *dcoMessages_);
  }

  // log problem information
  std::string sense = (dcoPar_->entry(DcoParams::objSense)==1.0) ? std::string("min") : std::string("min");
  dcoMessageHandler_->message(DISCO_PROBLEM_INFO,
                              *dcoMessages_)
    << problemName_
    << sense.c_str()
    << numCols_
    << numLinearRows_
    << matrix_->getNumElements()
    << numConicRows_
    << numIntegerCols_
    << CoinMessageEol;
}

// this should go into OsiConicSolverInterface or CoinUtils?
void DcoModel::readInstanceCbf(char const * dataFile) {
  // mps file reader
  DcoCbfIO * reader = new DcoCbfIO();
  reader->readCbf(dataFile);
  // set objective sense
  objSense_ = reader->objSense();
  // set dcoPar_
  dcoPar_->setEntry(DcoParams::objSense, objSense_);

  reader->getProblem(colLB_, colUB_, rowLB_, rowUB_, matrix_,
                     numConicRows_, coneStart_, coneMembers_, coneType_);
  numCols_ = matrix_->getNumCols();
  numLinearRows_ = matrix_->getNumRows();
  numRows_ = numLinearRows_ + numConicRows_;

  // add conic row bounds to rowLB and rowUB
  double * tempRowLB = new double[numRows_];
  std::copy(rowLB_, rowLB_+numLinearRows_, tempRowLB);
  std::fill_n(tempRowLB+numLinearRows_, numConicRows_, 0.0);
  delete[] rowLB_;
  rowLB_ = tempRowLB;
  tempRowLB = NULL;
  double * tempRowUB = new double[numRows_];
  std::copy(rowUB_, rowUB_+numLinearRows_, tempRowUB);
  std::fill_n(tempRowUB+numLinearRows_, numConicRows_, reader->getInfinity());
  delete[] rowUB_;
  rowUB_ = tempRowUB;
  tempRowUB = NULL;

  objCoef_ = new double [numCols_]();

  std::copy(reader->objCoef(), reader->objCoef()+reader->getNumCols(),
            objCoef_);

  // get integrality info
  numIntegerCols_ = reader->getNumInteger();
  integerCols_ = new int[numIntegerCols_];
  std::copy(reader->integerCols(), reader->integerCols()+numIntegerCols_,
            integerCols_);
  isInteger_ = new int[numCols_]();
  for (int i=0; i<numIntegerCols_; ++i) {
    isInteger_[integerCols_[i]] = 1;
  }

  delete reader;
}

void DcoModel::readInstanceMps(char const * dataFile) {
  // mps file reader
  CoinMpsIO * reader = new CoinMpsIO;
  // set reader log level
  //reader->messageHandler()->setLogLevel(dcoPar_->entry(DcoParams::logLevel));
  reader->messageHandler()->setLogLevel(0);
  reader->readMps(dataFile, "");
  numCols_ = reader->getNumCols();

  // allocate variable bounds
  colLB_ = new double [numCols_];
  colUB_ = new double [numCols_];
  // set col bounds
  std::copy(reader->getColLower(), reader->getColLower()+numCols_, colLB_);
  std::copy(reader->getColUpper(), reader->getColUpper()+numCols_, colUB_);

  // set objective sense
  // todo(aykut) we should ask reader about the objective sense
  objSense_ = dcoPar_->entry(DcoParams::objSense);

  // set objective coefficients
  objCoef_ = new double [numCols_];
  double const * reader_obj = reader->getObjCoefficients();
  std::copy(reader_obj, reader_obj+numCols_, objCoef_);

  // set integer columns
  // get variable integrality constraints
  numIntegerCols_ = 0;
  integerCols_ = new int[numCols_];
  isInteger_ = new int[numCols_];
  // may return NULL
  char const * is_integer = reader->integerColumns();
  if (is_integer!=NULL and strcmp(is_integer, "")) {
    for (int i=0; i<numCols_; ++i) {
      if (is_integer[i]) {
        integerCols_[numIntegerCols_++] = i;
        isInteger_[i] = 1;
      }
      else {
        isInteger_[i] = 0;
      }
    }
  }
  else {
    dcoMessageHandler_->message(3000, "Dco",
                                "CoinMpsIO::integerColumns() did return "
                                "NULL pointer for "
                                " column types! This looks like a CoinMpsIO "
                                "bug to me. Please report! I will use "
                                " CoinMpsIO::isContinuous(int index) function "
                                "instead.",
                                'W', 0)
      << CoinMessageEol;
    for (int i=0; i<numCols_; ++i) {
      if (!reader->isContinuous(i)) {
        integerCols_[numIntegerCols_++] = i;
        isInteger_[i] = 1;
      }
      else {
        isInteger_[i] = 0;
      }
    }
  }
  // resize integerCols_
  int * temp = new int[numIntegerCols_];
  std::copy(integerCols_, integerCols_+numIntegerCols_, temp);
  delete[] integerCols_;
  integerCols_ = temp;

  // read conic part
  int reader_return = reader->readConicMps(NULL, coneStart_, coneMembers_,
                                           coneType_, numConicRows_);
  // when there is no conic section, status is -3.
  if (reader_return==-3) {
    // no cones in problem
  }
  else if (reader_return!=0) {
    dcoMessageHandler_->message(DISCO_READ_MPSERROR,
                                *dcoMessages_) << reader_return
                                               << CoinMessageEol;
  }
  // store number of constraints in the problem
  numLinearRows_ = reader->getNumRows();
  numRows_ = numLinearRows_ + numConicRows_;

  // allocate row bounds
  rowLB_ = new double [numRows_];
  rowUB_ = new double [numRows_];
  // set row bounds for linear rows
  std::copy(reader->getRowLower(), reader->getRowLower()+numLinearRows_, rowLB_);
  std::copy(reader->getRowUpper(), reader->getRowUpper()+numLinearRows_, rowUB_);
  // set conic row bounds
  std::fill_n(rowLB_+numLinearRows_, numConicRows_, 0.0);
  std::fill_n(rowUB_+numLinearRows_, numConicRows_, reader->getInfinity());

  matrix_ = new CoinPackedMatrix(*reader->getMatrixByRow());
  problemName_ = reader->getProblemName();

  // free Coin MPS reader
  delete reader;
}


void DcoModel::readParameters(const int argnum,
                              const char * const * arglist) {
  AlpsPar()->readFromArglist(argnum, arglist);
  dcoPar_->readFromArglist(argnum, arglist);
  setMessageLevel();
}

// Add variables to *this.
void DcoModel::setupAddVariables() {
  // add variables
  BcpsVariable ** variables = new BcpsVariable*[numCols_];
  for (int i=0; i<numCols_; ++i) {
    variables[i] = new DcoVariable(i, colLB_[i], colUB_[i],
                                   colLB_[i], colUB_[i]);
    if (isInteger_[i]) {
      variables[i]->setIntType('I');
    }
    else {
      variables[i]->setIntType('C');
    }
    variables[i]->setBroker(broker_);
  }
  setVariables(variables, numCols_);
  // variables[i] are now owned by BcpsModel, do not free them.
  delete[] variables;
}

void DcoModel::setupAddLinearConstraints() {
  // == add constraints to *this
  BcpsConstraint ** constraints = new BcpsConstraint*[numLinearRows_];
  int const * indices = matrix_->getIndices();
  double const * values = matrix_->getElements();
  int const * lengths = matrix_->getVectorLengths();
  int const * starts = matrix_->getVectorStarts();
  for (int i=0; i<numLinearRows_; ++i) {
    constraints[i] = new DcoLinearConstraint(lengths[i], indices+starts[i],
                                             values+starts[i], rowLB_[i],
                                             rowUB_[i]);
    constraints[i]->setBroker(broker_);
  }
  setConstraints(constraints, numLinearRows_);
  // constraints[i] are owned by BcpsModel. Do not free them here.
  delete[] constraints;
}

void DcoModel::setupAddConicConstraints() {
  // iterate over cones and add them to the model
  for (int i=0; i<numConicRows_; ++i) {
    if (coneType_[i]!=1 and coneType_[i]!=2) {
      dcoMessageHandler_->message(DISCO_READ_CONEERROR,
                                  *dcoMessages_) << CoinMessageEol;
    }
    int num_members = coneStart_[i+1]-coneStart_[i];
    if (coneType_[i]==2 and num_members<3) {
      dcoMessageHandler_->message(DISCO_READ_ROTATEDCONESIZE,
                                  *dcoMessages_) << CoinMessageEol;
    }
    DcoLorentzConeType type = DcoLorentzCone;
    if (coneType_[i]==2) {
      type = DcoRotatedLorentzCone;
    }
    else if (coneType_[i]!=1) {
      dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE,
                                  *dcoMessages_)
        << __FILE__ << __LINE__ << CoinMessageEol;
    }
    DcoConicConstraint * cc =
      new DcoConicConstraint(type, num_members,
                             coneMembers_+coneStart_[i]);
    cc->setBroker(broker_);
    addConstraint(cc);
  }
}

/// Write out parameters.
void DcoModel::writeParameters(std::ostream& outstream) const {
  outstream << "\n================================================"
            <<std::endl;
  outstream << "ALPS Parameters: " << std::endl;
  AlpsPar_->writeToStream(outstream);
  outstream << "\n================================================"
            <<std::endl;
  outstream << "DISCO Parameters: " << std::endl;
  dcoPar_->writeToStream(outstream);
}

void DcoModel::preprocess() {
  // // set message levels
  // setMessageLevel();


  // // some bounds are improved if updated is true

  // // notes(aykut) this can only improve the bounds if the leading variable is
  // // bounded.  most of the time it is not. Things can get better if we use some
  // // other bound improvement first and it does improve the upper bound of
  // // leading variables.
  // //bool updated = DcoPresolve::improve_bounds(this);

  // write parameters used
  //writeParameters(std::cout);

  // approximation of cones will update numLinearRows_, numRows_, rowLB_,
  // rowUB_, matrix_.
  if (numConicRows_ > 0) {
    approximateCones();
  }
}

void DcoModel::approximateCones() {
#ifdef __OA__
  // need to load problem to the solver.

  // load problem to the solver
  solver_->loadProblem(*matrix_, colLB_, colUB_, objCoef_,
                       rowLB_, rowUB_);
  bool dual_infeasible = false;
  int iter = 0;
  int ipm_iter;
  int oa_iter;
  int num_ipm_cuts = 0;
  int num_oa_cuts = 0;
  // solve problem
  solver_->resolve();
  // get cone data in the required form
  // todo(aykut) think about updating cut library for the input format
  OsiLorentzConeType * coneTypes = new OsiLorentzConeType[numConicRows_];
  int * coneSizes = new int[numConicRows_];
  int const ** coneMembers = new int const *[numConicRows_];
  for (int i=0; i<numConicRows_; ++i) {
    if (coneType_[i]==1) {
      coneTypes[i] = OSI_QUAD;
    }
    else if (coneType_[i]==2) {
      coneTypes[i] = OSI_RQUAD;
    }
    else {
      dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE, *dcoMessages_)
        << __FILE__ << __LINE__ << CoinMessageEol;
    }
    coneSizes[i] = coneStart_[i+1]-coneStart_[i];
    coneMembers[i] = coneMembers_ + coneStart_[i];
  }
  // used to decide on number of iterations in outer approximation
  int largest_cone_size = *std::max_element(coneSizes,
                                            coneSizes+numConicRows_);
  do {
    // generate cuts
    OsiCuts * ipm_cuts = new OsiCuts();
    OsiCuts * oa_cuts = new OsiCuts();
    CglConicCutGenerator * cg_ipm = new CglConicIPM();
    CglConicCutGenerator * cg_oa =
      new CglConicOA(dcoPar_->entry(DcoParams::coneTol));
    // get cone info
    cg_ipm->generateCuts(*solver_, *ipm_cuts, numConicRows_, coneTypes,
                         coneSizes, coneMembers, largest_cone_size);
    // cg_oa->generateCuts(*solver_, *oa_cuts, numCoreCones_, coneTypes_,
    //                   coneSizes_, coneMembers_, largest_cone_size);
    // if we do not get any cuts break the loop
    if (ipm_cuts->sizeRowCuts()==0 && oa_cuts->sizeRowCuts()==0) {
      break;
    }
    // if problem is unbounded do nothing, add cuts to the problem
    // this will make lp relaxation infeasible
    solver_->applyCuts(*ipm_cuts);
    solver_->applyCuts(*oa_cuts);
    solver_->resolve();
    num_ipm_cuts += ipm_cuts->sizeRowCuts();
    delete ipm_cuts;
    delete oa_cuts;
    delete cg_ipm;
    delete cg_oa;
    dual_infeasible = solver_->isProvenDualInfeasible();
    iter++;
  } while(dual_infeasible);
  // add outer approximating cuts for DcoParams::approxNumPass (default is 50)
  // many rounds
  ipm_iter = iter;
  iter = 0;
  int oa_iter_limit = dcoPar_->entry(DcoParams::approxNumPass);
  while(iter<oa_iter_limit) {
    OsiCuts * oa_cuts = new OsiCuts();
    CglConicCutGenerator * cg_oa =
      new CglConicOA(dcoPar_->entry(DcoParams::coneTol));
    cg_oa->generateCuts(*solver_, *oa_cuts, numConicRows_, coneTypes,
                        coneSizes, coneMembers, 1);
    int num_cuts = oa_cuts->sizeRowCuts();
    num_oa_cuts += num_cuts;
    if (num_cuts==0) {
      // ifno cuts are produced break early
      delete oa_cuts;
      delete cg_oa;
      break;
    }
    solver_->applyCuts(*oa_cuts);
    solver_->resolve();
    delete oa_cuts;
    delete cg_oa;
    iter++;
  }
  oa_iter = iter;
  std::cout << "===== Preprocessing Summary =====" << std::endl;
  std::cout << "IPM iterations " << ipm_iter << std::endl;
  std::cout << "IPM cuts " << num_ipm_cuts << std::endl;
  std::cout << "OA iterations " << oa_iter << std::endl;
  std::cout << "OA cuts " << num_oa_cuts << std::endl;
  std::cout << "Linear relaxation objective value "
            << solver_->getObjValue() << std::endl;
  std::cout << "=================================" << std::endl;
  delete[] coneTypes;
  delete[] coneSizes;
  delete[] coneMembers;

  double cutOaSlack = dcoPar_->entry(DcoParams::cutOaSlack1);
  // print cut activity
  // {
  //   int numCuts = solver_->getNumRows() - numLinearRows_;
  //   double const * activity = solver_->getRowActivity();
  //   double const * lb = solver_->getRowLower();
  //   double const * ub = solver_->getRowUpper();
  //   CoinWarmStartBasis const * ws =
  //     dynamic_cast<CoinWarmStartBasis*> (solver_->getWarmStart());
    //for (int i=numLinearRows_; i<solver_->getNumRows(); ++i) {
      // double norm = 0.0;
      // {
      //   // compute norm of cut
      //   CoinPackedMatrix const * mat = solver_->getMatrixByRow();
      //   int start = mat->getVectorStarts()[i];
      //   int size = mat->getVectorLengths()[i];
      //   double const * value = mat->getElements() + start;
      //   for (int j=0; j<size; ++j) {
      //     norm += fabs(value[j]);
      //   }
      // }

      // std::cout << "row: " << std::setw(4) << i
      //           << " lb: " << std::setw(12) << lb[i]
      //           << " activity: " << std::setw(12) << activity[i]
      //           << " ub: " << std::setw(12) << ub[i]
      //           << " status: " << (ws->getArtifStatus(i)==CoinWarmStartBasis::basic)
      //           << " remove: " << ((ub[i]-activity[i])>cutOaSlack)
      //           << " norm: " << norm
      //           << std::endl;
    //}
  //}


  // remove inactive cuts
  {
    int numCuts = solver_->getNumRows() - numLinearRows_;
    // get warm start basis
    CoinWarmStartBasis * ws =
      dynamic_cast<CoinWarmStartBasis*> (solver_->getWarmStart());
    if (ws==NULL) {
      // nothing to do if there is no warm start information
      std::cerr << "Disco warning: No warm start object exists in solver. "
                << "Unable to clean cuts." << std::endl;
    }
    else if (numCuts==0) {
      // nothing to do if there no any cuts
    }
    else {
      int numDel = 0;
      int * delInd = new int[numCuts];
      double const * activity = solver_->getRowActivity();
      //double const * lb = solver_->getRowLower();
      double const * ub = solver_->getRowUpper();

      // iterate over cuts and colect inactive ones
      for (int i=0; i<numCuts; ++i) {
        if (ws->getArtifStatus(numLinearRows_+i) == CoinWarmStartBasis::basic) {
          // cut is inactive
          // do we really want to remove? check how much it is inactive
          if (ub[numLinearRows_+i] - activity[numLinearRows_+i] > cutOaSlack) {
            // remove since it is inactive too much
            delInd[numDel++] = i+numLinearRows_;
          }
        }
      }
      // delete cuts from solver
      if (numDel) {
        std::cout << "Approx cones: "
                  << " removed: " << numDel
                  << " remain: " << numCuts-numDel
                  << std::endl;
        solver_->deleteRows(numDel, delInd);
        // resolve to correct status
        solver_->resolve();
      }
      delete[] delInd;
    }
    delete ws;
  }
  initOAcuts_ = solver_->getNumRows() - numLinearRows_;

  // get updated data from solver
  // delete matrix_;
  // matrix_ = new CoinPackedMatrix(*solver_->getMatrixByRow());
  // // update number of rows
  // numLinearRows_ = solver_->getNumRows();
  // numRows_ = numLinearRows_+numConicRows_;
  // // update row bounds
  // delete[] rowLB_;
  // rowLB_ = new double[numRows_];
  // std::copy(solver_->getRowLower(),
  //           solver_->getRowLower()+numLinearRows_, rowLB_);
  // std::fill_n(rowLB_+numLinearRows_, numConicRows_, 0.0);
  // delete[] rowUB_;
  // rowUB_ = new double[numRows_];
  // std::copy(solver_->getRowUpper(),
  //           solver_->getRowUpper()+numLinearRows_, rowUB_);
  // std::fill_n(rowUB_+numLinearRows_, numConicRows_, solver_->getInfinity());
#endif
}

//todo(aykut) why does this return to bool?
// should be fixed in Alps level.

// relase redundant memory once the fields are set.
bool DcoModel::setupSelf() {
  // set message level again. We call this in readInstance() too.
  // this is so due to parallelization.
  setMessageLevel();
  // set relaxed array for integer columns
  numRelaxedCols_ = numIntegerCols_;
  relaxedCols_ = new int[numRelaxedCols_];
  std::copy(integerCols_, integerCols_+numIntegerCols_,
            relaxedCols_);
  // set iteration count to 0
  numRelaxIterations_ = 0;
#ifdef __OA__
  solver_->reset();
  solver_->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo, NULL);
  // clp specific options for getting unboundedness directions
  dynamic_cast<OsiClpSolverInterface*>(solver_)->getModelPtr()->setMoreSpecialOptions(0);
  dynamic_cast<OsiClpSolverInterface*>(solver_)->getModelPtr()->setLogLevel(0);
#else
  solver_->setHintParam(OsiDoReducePrint, true, OsiHintTry);
#endif
  // load problem to the solver
  solver_->loadProblem(*matrix_, colLB_, colUB_, objCoef_,
                       rowLB_, rowUB_);
  // set integer variables
  solver_->setInteger(integerCols_, numIntegerCols_);

#if defined(__OA__)
  // we relax conic constraints when OA is used.
  // set relaxed array for conic constraints
  numRelaxedRows_ = numConicRows_;
  relaxedRows_ = new int[numRelaxedRows_];
  // notes(aykut) we assume conic rows start after linear rows.
  for (int i=0; i<numRelaxedRows_; ++i) {
    relaxedRows_[i] = numLinearRows_+i;
  }

  // set leading variable lower bounds to 0
  for (int i=0; i<numConicRows_; ++i) {
    if (coneType_[i]==1) {
      colLB_[coneMembers_[coneStart_[i]]] = 0.0;
      solver_->setColLower(coneMembers_[coneStart_[i]], 0.0);
    }
    else if (coneType_[i]==2) {
      colLB_[coneMembers_[coneStart_[i]]] = 0.0;
      colLB_[coneMembers_[coneStart_[i]+1]] = 0.0;
      solver_->setColLower(coneMembers_[coneStart_[i]], 0.0);
      solver_->setColLower(coneMembers_[coneStart_[i]+1], 0.0);
    }
    else {
      dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE,
                                  *dcoMessages_)
        << coneType_[i] << CoinMessageEol;
    }
  }
#else
  // add conic constraints to the solver
  for (int i=0; i<numConicRows_; ++i) {
    // do not relax conic constraints, add them to the conic solver
    OsiLorentzConeType osi_type;
    if (coneType_[i]==1) {
      osi_type = OSI_QUAD;
    }
    else if (coneType_[i]==2) {
      osi_type = OSI_RQUAD;
    }
    solver_->addConicConstraint(osi_type, coneStart_[i+1]-coneStart_[i],
                              coneMembers_+coneStart_[i]);
  }
#endif

  // create disco variables
  setupAddVariables();
  // create disco constraints, linear
  setupAddLinearConstraints();
  // create disco constraints, conic
  setupAddConicConstraints();

  // set branch strategy
  setBranchingStrategy();

  // add constraint generators
#ifdef __OA__
  addConstraintGenerators();
#endif

  // add heuristics
  addHeuristics();

  // set cutoff for solvers if user provides one.
  double cutoff = dcoPar_->entry(DcoParams::cutoff);
  if (cutoff!=ALPS_INC_MAX) {
    solver_->setDblParam(OsiDualObjectiveLimit, objSense_*cutoff);
  }
  return true;
}

// set message level
void DcoModel::setMessageLevel() {
  // get Alps log level
  int alps_log_level = AlpsPar()->entry(AlpsParams::msgLevel);
  int dco_log_level = dcoPar_->entry(DcoParams::logLevel);

#if defined(DISCO_DEBUG) || defined(DISCO_DEBUG_BRANCH) ||defined(DISCO_DEBUG_CUT) || defined(DISCO_DEBUG_PROCESS)
  // reset Disco log level, we will set it depending on the debug macros are
  // defined.
  dco_log_level = 0;
#endif

#ifdef DISCO_DEBUG_BRANCH
  dco_log_level += DISCO_DLOG_BRANCH;
#endif

#ifdef DISCO_DEBUG_CUT
  dco_log_level += DISCO_DLOG_CUT;
#endif

#ifdef DISCO_DEBUG_PROCESS
  dco_log_level += DISCO_DLOG_PROCESS;
#endif

#ifdef DISCO_DEBUG
  // debug branch, cut, process
  dco_log_level = DISCO_DLOG_BRANCH;
  dco_log_level += DISCO_DLOG_CUT;
  dco_log_level += DISCO_DLOG_PROCESS;
#endif

  // todo(aykut) create different parameters for Alps and Bcps.
  broker()->messageHandler()->setLogLevel(alps_log_level);
  bcpsMessageHandler()->setLogLevel(alps_log_level);
  dcoMessageHandler_->setLogLevel(dco_log_level);
}

void DcoModel::addConstraintGenerators() {
  // get cut strategy
  cutStrategy_ = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutStrategy));
  // get cut generation frequency
  cutGenerationFrequency_ =
    dcoPar_->entry(DcoParams::cutGenerationFrequency);
  if (cutGenerationFrequency_ < 1) {
    // invalid cut fraquency given, change it to 1.
    dcoMessageHandler_->message(DISCO_INVALID_CUT_FREQUENCY,
                                *dcoMessages_)
      << cutGenerationFrequency_
      << 1
      << CoinMessageEol;
    cutGenerationFrequency_ = 1;
  }

  // get cut strategies from parameters
  DcoCutStrategy cliqueStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutCliqueStrategy));
  DcoCutStrategy fCoverStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutFlowCoverStrategy));
  DcoCutStrategy gomoryStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutGomoryStrategy));
  DcoCutStrategy knapStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutKnapsackStrategy));
  DcoCutStrategy mirStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutMirStrategy));
  DcoCutStrategy oddHoleStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutOddHoleStrategy));
  DcoCutStrategy probeStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutProbingStrategy));
  DcoCutStrategy twoMirStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutTwoMirStrategy));
  DcoCutStrategy ipmStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutIpmStrategy));
  DcoCutStrategy ipmintStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutIpmIntStrategy));
  DcoCutStrategy oaStrategy = static_cast<DcoCutStrategy>
    (dcoPar_->entry(DcoParams::cutOaStrategy));

  // get cut frequencies from parameters
  int cliqueFreq = dcoPar_->entry(DcoParams::cutCliqueFreq);
  int fCoverFreq = dcoPar_->entry(DcoParams::cutFlowCoverFreq);
  int gomoryFreq = dcoPar_->entry(DcoParams::cutGomoryFreq);
  int knapFreq = dcoPar_->entry(DcoParams::cutKnapsackFreq);
  int mirFreq = dcoPar_->entry(DcoParams::cutMirFreq);
  int oddHoleFreq = dcoPar_->entry(DcoParams::cutOddHoleFreq);
  int probeFreq = dcoPar_->entry(DcoParams::cutProbingFreq);
  int twoMirFreq = dcoPar_->entry(DcoParams::cutTwoMirFreq);
  int ipmFreq = dcoPar_->entry(DcoParams::cutIpmFreq);
  int ipmintFreq = dcoPar_->entry(DcoParams::cutIpmIntFreq);
  int oaFreq = dcoPar_->entry(DcoParams::cutOaFreq);

  //----------------------------------
  // Add cut generators.
  //----------------------------------
  // Add probe cut generator
  if (probeStrategy == DcoCutStrategyNotSet) {
    // Disable by default
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      probeStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      probeStrategy = cutStrategy_;
      probeFreq = cutGenerationFrequency_;
    }
    else {
      probeStrategy = cutStrategy_;
    }
  }
  if (probeStrategy != DcoCutStrategyNone) {
    CglProbing *probing = new CglProbing;
    probing->setUsingObjective(true);
    probing->setMaxPass(1);
    probing->setMaxPassRoot(5);
    // Number of unsatisfied variables to look at
    probing->setMaxProbe(10);
    probing->setMaxProbeRoot(1000);
    // How far to follow the consequences
    probing->setMaxLook(50);
    probing->setMaxLookRoot(500);
    // Only look at rows with fewer than this number of elements
    probing->setMaxElements(200);
    probing->setRowCuts(3);
    addConGenerator(probing, DcoConstraintTypeProbe, probeStrategy, probeFreq);
  }

  // Add clique cut generator.
  if (cliqueStrategy == DcoCutStrategyNotSet) {
    // Disable by default
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      cliqueStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      cliqueFreq = cutGenerationFrequency_;
      cliqueStrategy = DcoCutStrategyPeriodic;
    }
    else { // Root or Auto
      cliqueStrategy = cutStrategy_;
    }
  }
  if (cliqueStrategy != DcoCutStrategyNone) {
    CglClique *cliqueCut = new CglClique ;
    cliqueCut->setStarCliqueReport(false);
    cliqueCut->setRowCliqueReport(false);
    addConGenerator(cliqueCut, DcoConstraintTypeClique, cliqueStrategy,
                    cliqueFreq);
  }

  // Add odd hole cut generator.
  if (oddHoleStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // Disable by default
      oddHoleStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      oddHoleStrategy = DcoCutStrategyPeriodic;
      oddHoleFreq = cutGenerationFrequency_;
    }
    else {
      oddHoleStrategy = cutStrategy_;
    }
  }
  if (oddHoleStrategy != DcoCutStrategyNone) {
    CglOddHole *oldHoleCut = new CglOddHole;
    oldHoleCut->setMinimumViolation(0.005);
    oldHoleCut->setMinimumViolationPer(0.00002);
    // try larger limit
    oldHoleCut->setMaximumEntries(200);
    addConGenerator(oldHoleCut, DcoConstraintTypeOddHole, oddHoleStrategy,
                    oddHoleFreq);
  }

  // Add flow cover cut generator.
  if (fCoverStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // Disable by default
      fCoverStrategy = DcoCutStrategyNone;
      fCoverFreq = cutGenerationFrequency_;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      fCoverStrategy = cutStrategy_;
      fCoverFreq = cutGenerationFrequency_;
    }
    else {
      fCoverStrategy = cutStrategy_;
    }
  }
  if (fCoverStrategy != DcoCutStrategyNone) {
    CglFlowCover *flowGen = new CglFlowCover;
    addConGenerator(flowGen, DcoConstraintTypeFCover, fCoverStrategy,
                    fCoverFreq);
  }

  // Add knapsack cut generator.
  if (knapStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // Disable by default
      knapStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      knapStrategy = cutStrategy_;
      knapFreq = cutGenerationFrequency_;
    }
    else {
      knapStrategy = cutStrategy_;
    }
  }
  if (knapStrategy != DcoCutStrategyNone) {
    CglKnapsackCover *knapCut = new CglKnapsackCover;
    addConGenerator(knapCut, DcoConstraintTypeKnap, knapStrategy, knapFreq);
  }

  // Add MIR cut generator.
  if (mirStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // Disable by default
      mirStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      mirStrategy = cutStrategy_;
      mirFreq = cutGenerationFrequency_;
    }
    else {
      mirStrategy = cutStrategy_;
    }
  }
  if (mirStrategy != DcoCutStrategyNone) {
    CglMixedIntegerRounding2 *mixedGen = new CglMixedIntegerRounding2;
    addConGenerator(mixedGen, DcoConstraintTypeMIR, mirStrategy, mirFreq);
  }

  // Add Gomory cut generator.
  if (gomoryStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // Disable by default
      gomoryStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      gomoryStrategy = cutStrategy_;
      gomoryFreq = cutGenerationFrequency_;
    }
    else {
      gomoryStrategy = cutStrategy_;
    }
  }
  if (gomoryStrategy != DcoCutStrategyNone) {
    CglGomory *gomoryCut = new CglGomory;
    //CglGMI *gomoryCut = new CglGMI;
    // try larger limit
    gomoryCut->setLimit(40);
    addConGenerator(gomoryCut, DcoConstraintTypeGomory, gomoryStrategy,
                    gomoryFreq);
  }

  // Add Tow MIR cut generator.
  // Disable forever, not useful.
  twoMirStrategy = DcoCutStrategyNone;
  if (twoMirStrategy != DcoCutStrategyNone) {
    CglTwomir *twoMirCut =  new CglTwomir;
    addConGenerator(twoMirCut, DcoConstraintTypeTwoMIR, twoMirStrategy,
                    twoMirFreq);
  }

  // Add IPM cut generator
  if (ipmStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // Disable by default
      ipmStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      ipmStrategy = cutStrategy_;
      ipmFreq = cutGenerationFrequency_;
    }
    else {
      ipmStrategy = cutStrategy_;
    }
  }
  if (ipmStrategy != DcoCutStrategyNone) {
    CglConicCutGenerator * ipm_gen = new CglConicIPM();
    addConGenerator(ipm_gen, DcoConstraintTypeIPM, ipmStrategy, ipmFreq);
  }

  // Add IPM integer cut generator
  if (ipmintStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // Disable by default
      ipmintStrategy = DcoCutStrategyNone;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      ipmintStrategy = cutStrategy_;
      ipmintFreq = cutGenerationFrequency_;
    }
    else {
      ipmintStrategy = cutStrategy_;
    }
  }
  if (ipmintStrategy != DcoCutStrategyNone) {
    CglConicCutGenerator * ipm_int_gen = new CglConicIPMint();
    addConGenerator(ipm_int_gen, DcoConstraintTypeIPMint, ipmintStrategy,
                    ipmintFreq);
  }

  // Add Outer approximation cut generator
  if (oaStrategy == DcoCutStrategyNotSet) {
    if (cutStrategy_ == DcoCutStrategyNotSet) {
      // periodic by default
      oaStrategy = DcoCutStrategyPeriodic;
    }
    else if (cutStrategy_ == DcoCutStrategyPeriodic) {
      oaStrategy = cutStrategy_;
      oaFreq = cutGenerationFrequency_;
    }
    else {
      oaStrategy = cutStrategy_;
    }
  }
  if (oaStrategy != DcoCutStrategyNone && numConicRows_) {
    CglConicCutGenerator * oa_gen =
      new CglConicOA(dcoPar_->entry(DcoParams::coneTol));
    addConGenerator(oa_gen, DcoConstraintTypeOA, oaStrategy, oaFreq);
  }

  // Adjust cutStrategy_ according to the strategies of each cut generators.
  // set it to the most allowing one.
  // if there is at least one periodic strategy, set it to periodic.
  // if no periodic and there is at least one root, set it to root
  // set it to None otherwise.
  cutStrategy_ = DcoCutStrategyNone;
  cutGenerationFrequency_ = 100;
  bool periodic_exists = false;
  bool root_exists = false;
  std::map<DcoConstraintType, DcoConGenerator*>::iterator it;
  for (it=conGenerators_.begin(); it!=conGenerators_.end(); ++it) {
    DcoCutStrategy curr = it->second->strategy();
    if (curr==DcoCutStrategyPeriodic) {
      periodic_exists = true;
      break;
    }
    else if (curr==DcoCutStrategyRoot) {
      root_exists = true;
    }
  }
  if (periodic_exists) {
    // greatest divisor of all frequencies
    cutStrategy_ = DcoCutStrategyPeriodic;
    cutGenerationFrequency_ = 1;
  }
  else if (root_exists) {
    cutStrategy_ = DcoCutStrategyRoot;
    // this is not relevant, since we will generate only in root.
    cutGenerationFrequency_ = 100;
  }
}


void DcoModel::addConGenerator(CglCutGenerator * cgl_gen,
                               DcoConstraintType type,
                               DcoCutStrategy dco_strategy,
                               int frequency) {
  char const * name = dcoConstraintTypeName[type];
  DcoConGenerator * con_gen = new DcoLinearConGenerator(this, cgl_gen, type,
                                                        name,
                                                        dco_strategy,
                                                        frequency);
  conGenerators_[type] = con_gen;
}

/// Add constraint generator.
void DcoModel::addConGenerator(CglConicCutGenerator * cgl_gen,
                               DcoConstraintType type,
                               DcoCutStrategy dco_strategy,
                               int frequency) {
  char const * name = dcoConstraintTypeName[type];
  DcoConGenerator * con_gen = new DcoConicConGenerator(this, cgl_gen, type,
                                                       name,
                                                       dco_strategy,
                                                       frequency);
  conGenerators_[type] = con_gen;
}

void DcoModel::addHeuristics() {
  // todo(aykut) this function (heuristic adding process) can be improved
  // since global parameters are ignored with this design.

  heuristics_.clear();
  // get global heuristic strategy
  heurStrategy_ = static_cast<DcoHeurStrategy>
    (dcoPar_->entry(DcoParams::heurStrategy));
  // get global heuristic call frequency
  heurFrequency_ =
    dcoPar_->entry(DcoParams::heurCallFrequency);
  if (heurFrequency_ < 1) {
    // invalid heur fraquency given, change it to 1.
    dcoMessageHandler_->message(DISCO_INVALID_HEUR_FREQUENCY,
                                *dcoMessages_)
      << heurFrequency_
      << 1
      << CoinMessageEol;
    heurFrequency_ = 1;
  }

  // get rounding heuristics strategy from parameters
  DcoHeurStrategy roundingStrategy = static_cast<DcoHeurStrategy>
    (dcoPar_->entry(DcoParams::heurRoundStrategy));
  // get rounding heuristic frequency from parameters
  int roundingFreq = dcoPar_->entry(DcoParams::heurRoundFreq);

  // add heuristics
  // == add rounding heuristics
  if (roundingStrategy != DcoHeurStrategyNone) {
    DcoHeuristic * round = new DcoHeurRounding(this, "rounding",
                                               roundingStrategy, roundingFreq);
    heuristics_.push_back(round);
  }


  // Adjust heurStrategy_ according to the strategies/frequencies of each
  // heuristic. Set it to the most allowing one.
  // if there is at least one periodic strategy, set it to periodic.
  // if no periodic and there is at least one root, set it to root
  // set it to None otherwise.
  heurStrategy_ = DcoHeurStrategyNone;
  heurFrequency_ = -1;
  bool periodic_exists = false;
  bool root_exists = false;
  std::vector<DcoHeuristic*>::iterator it;
  for (it=heuristics_.begin(); it!=heuristics_.end(); ++it) {
    DcoHeurStrategy curr = (*it)->strategy();
    if (curr==DcoHeurStrategyPeriodic) {
      periodic_exists = true;
      break;
    }
    else if (curr==DcoHeurStrategyRoot) {
      root_exists = true;
    }
  }
  if (periodic_exists) {
    heurStrategy_ = DcoHeurStrategyPeriodic;
    heurFrequency_ = 1;
  }
  else if (root_exists) {
    heurStrategy_ = DcoHeurStrategyRoot;
    // this is not relevant, since we will generate only in root.
    heurFrequency_ = -1;
  }
}

void DcoModel::setBranchingStrategy() {
    // set branching startegy
  int brStrategy = static_cast<DcoBranchingStrategy>
    (dcoPar_->entry(DcoParams::branchStrategy));
  switch(brStrategy) {
  case DcoBranchingStrategyMaxInfeasibility:
    branchStrategy_ = new DcoBranchStrategyMaxInf(this);
    break;
  case DcoBranchingStrategyPseudoCost:
    branchStrategy_ = new DcoBranchStrategyPseudo(this);
    break;
  // case DcoBranchingStrategyReliability:
  //   branchStrategy_ = new DcoBranchStrategyRel(this, reliability);
  //   break;
  case DcoBranchingStrategyStrong:
     branchStrategy_ = new DcoBranchStrategyStrong(this);
     break;
  // case DcoBranchingStrategyBilevel:
  //   branchStrategy_ = new DcoBranchStrategyBilevel(this);
  //   break;
  default:
    dcoMessageHandler_->message(DISCO_UNKNOWN_BRANCHSTRATEGY,
                                *dcoMessages_) << brStrategy << CoinMessageEol;
    throw CoinError("Unknown branch strategy.", "setupSelf","DcoModel");
  }

  // set ramp up branch strategy
  brStrategy = static_cast<DcoBranchingStrategy>
    (dcoPar_->entry(DcoParams::branchStrategyRampUp));
  switch(brStrategy) {
  case DcoBranchingStrategyMaxInfeasibility:
    rampUpBranchStrategy_ = new DcoBranchStrategyMaxInf(this);
    break;
  case DcoBranchingStrategyPseudoCost:
    rampUpBranchStrategy_ = new DcoBranchStrategyPseudo(this);
    break;
  // case DcoBranchingStrategyReliability:
  //   rampUpBranchStrategy_ = new DcoBranchStrategyRel(this, reliability);
  //   break;
  case DcoBranchingStrategyStrong:
     rampUpBranchStrategy_ = new DcoBranchStrategyStrong(this);
     break;
  // case DcoBranchingStrategyBilevel:
  //   rampUpBranchStrategy_ = new DcoBranchStrategyBilevel(this);
  //   break;
  default:
    dcoMessageHandler_->message(DISCO_UNKNOWN_BRANCHSTRATEGY,
                                *dcoMessages_) << brStrategy << CoinMessageEol;
    throw std::exception();
    throw CoinError("Unknown branch strategy.", "setupSelf","DcoModel");
  }
}

void DcoModel::postprocess() {
}


AlpsTreeNode * DcoModel::createRoot() {
  DcoTreeNode * root = new DcoTreeNode();
  DcoNodeDesc * desc = new DcoNodeDesc(this);
  root->setDesc(desc);
  std::vector<BcpsVariable *> & cols = getVariables();
  std::vector<BcpsConstraint *> & rows = getConstraints();
  int numCols = getNumCoreVariables();
  int numRows = getNumCoreConstraints();
  int * varIndices1 = new int [numCols];
  int * varIndices2 = new int [numCols];
  int * varIndices3 = NULL;
  int * varIndices4 = NULL;
  double * vlhe = new double [numCols];
  double * vuhe = new double [numCols];
  double * vlse = NULL;
  double * vuse = NULL;
  int * conIndices1 = new int [numRows];
  int * conIndices2 = new int [numRows];
  int * conIndices3 = NULL;
  int * conIndices4 = NULL;
  double * clhe = new double [numRows];
  double * cuhe = new double [numRows];
  double * clse = NULL;
  double * cuse = NULL;
  //-------------------------------------------------------------
  // Get var bounds and indices.
  //-------------------------------------------------------------
  for (int i=0; i<numCols; ++i) {
    vlhe[i] = cols[i]->getLbHard();
    vuhe[i] = cols[i]->getUbHard();
    varIndices1[i] = i;
    varIndices2[i] = i;
  }
  //-------------------------------------------------------------
  // Get con bounds and indices.
  //-------------------------------------------------------------
  for (int i=0; i<numRows; ++i) {
    clhe[i] = rows[i]->getLbHard();
    cuhe[i] = rows[i]->getUbHard();
    conIndices1[i] = i;
    conIndices2[i] = i;
  }
  int * tempInd = NULL;
  BcpsObject ** tempObj = NULL;
  desc->assignVars(0, tempInd,
                   0, tempObj,
                   false, numCols, varIndices1, vlhe, /*Var hard lb*/
                   false, numCols, varIndices2, vuhe, /*Var hard ub*/
                   false, 0, varIndices3, vlse,       /*Var soft lb*/
                   false, 0, varIndices4, vuse);      /*Var soft ub*/
  desc->assignCons(0, tempInd,
                   0, tempObj,
                   false, numRows, conIndices1, clhe, /*Con hard lb*/
                   false, numRows, conIndices2, cuhe, /*Con hard ub*/
                   false, 0,conIndices3,clse,         /*Con soft lb*/
                   false, 0,conIndices4,cuse);        /*Con soft ub*/
  //-------------------------------------------------------------
  // Mark it as an explicit node.
  //-------------------------------------------------------------
  root->setExplicit(1);
  return root;
}


DcoSolution * DcoModel::feasibleSolution(int & numInfColumns,
                                         double  & colInf,
                                         int & numInfRows,
                                         double & rowInf) {
  // set stats to 0
  numInfColumns = 0;
  numInfRows = 0;

  // for debug purposes we will keep the largest column and row infeasibility
  colInf = 0.0;
  rowInf = 0.0;

  // check feasibility of relxed columns, ie. integrality constraints
  // get vector of variables
  std::vector<BcpsVariable*> & cols = getVariables();
  for (int i=0; i<numRelaxedCols_; ++i) {
    // get column relaxedCol_[i]
    DcoVariable * curr = dynamic_cast<DcoVariable*> (cols[relaxedCols_[i]]);
    // check feasibility
    int preferredDir;
    double infeas = curr->infeasibility(this, preferredDir);
    if (infeas>0) {
      numInfColumns++;
      if (colInf<infeas) {
        colInf = infeas;
      }
    }
  }

  // check feasibility of relaxed rows
  // get vector of constraints
  std::vector<BcpsConstraint*> & rows = getConstraints();
  for (int i=0; i<numRelaxedRows_; ++i) {
    // get row relaxedRows_[i]
    DcoConstraint * curr = dynamic_cast<DcoConstraint*> (rows[relaxedRows_[i]]);
    // check feasibility
    int preferredDir;
    double infeas = curr->infeasibility(this, preferredDir);
    if (infeas>0) {
      numInfRows++;
      if (rowInf<infeas) {
        rowInf = infeas;
      }
    }
  }
  // report largest column and row infeasibilities
  dcoMessageHandler_->message(DISCO_INFEAS_REPORT, *dcoMessages_)
    << broker()->getProcRank()
    << colInf
    << rowInf
    << CoinMessageEol;

  // create DcoSolution instance if feasbile
  DcoSolution * dco_sol = 0;
  if (numInfColumns==0 && numInfRows==0) {
    double const * sol = solver()->getColSolution();
    double quality = solver()->getObjValue();
    dco_sol = new DcoSolution(numCols_, sol, quality);
    dco_sol->setBroker(broker_);

    // log debug information
    dcoMessageHandler_->message(DISCO_SOL_FOUND, *dcoMessages_)
      << broker()->getProcRank()
      << quality
      << CoinMessageEol;
  }
  return dco_sol;
}

//todo(aykut) When all node bounds are worse than incumbent solution
// this function reports negative gap.
// this happens since Alps takes nodes that will be fathomed into account,
// in (getBestNode function).
void DcoModel::nodeLog(AlpsTreeNode * node, bool force) {
  if ((broker_->getProcType() != AlpsProcessTypeMaster) &&
      (broker_->getProcType() != AlpsProcessTypeSerial)) {
    // do nothing if not serial code nor master processor
    return;
  }
  // todo(aykut) somehow AlpsKnowledgeBrokerMPI does not call this
  // function properly. It calls this once and only header is printed.
  // Check this. Disable in parallel mode.
  if (broker_->getProcType() == AlpsProcessTypeMaster) {
    return ;
  }
  // number of processed nodes
  int num_processed = broker()->getNumNodesProcessed();
  // number of nodes left
  int num_left  = broker()->updateNumNodesLeft();
  // log interval
  int interval =
    broker()->getModel()->AlpsPar()->entry(AlpsParams::nodeLogInterval);
  // need to print header if this is the first call of this function.
  dcoMessageHandler_->setPrefix(0);
  if (broker()->getNumNodeLog()==0 or
      broker()->getNumNodeLog()%50==0) {
    //todo(aykut) fix this function's documentation in Alps. It does not
    // count how many times this function is called, but how many times
    // this function is logged.
    broker()->setNumNodeLog(broker()->getNumNodeLog()+1);
    dcoMessageHandler_->message(DISCO_NODE_LOG_HEADER,
                                *dcoMessages_)
      << CoinMessageEol;
  }
  else if (force || (num_processed % interval == 0)) {
    double lb = ALPS_INFINITY;
    AlpsTreeNode * best_node = broker()->getBestNode();
    if (best_node) {
      lb = best_node->getQuality();
    }
    else {
      // no nodes stored in the broker.
    }
    if (broker()->hasKnowledge(AlpsKnowledgeTypeSolution)) {
      double ub = broker()->getIncumbentValue();
      // check whether lb is larger than ub, this might happen if there are
      // nodes that are not fathomed yet.
      if (lb>ub) {
        lb = ub;
      }
      std::stringstream lb_ss;
      lb_ss << std::setw(14) << std::left << std::scientific << lb;
      broker()->setNumNodeLog(broker()->getNumNodeLog()+1);
      double gap = 100.0*((ub-lb)/fabs(ub));
      // we need to convert lb, ub, gap into strings, CoinMessageHandler is
      // not capable for this yet.
      std::stringstream ub_ss;
      ub_ss << std::setw(14) << std::left << std::scientific << ub;
      std::stringstream gap_ss;
      gap_ss.precision(2);
      gap_ss << std::setw(6) << std::fixed << std::left << gap;
      dcoMessageHandler_->message(DISCO_NODE_LOG, *dcoMessages_)
        << num_processed
        << num_left
        << lb_ss.str().c_str()
        << ub_ss.str().c_str()
        << gap_ss.str().c_str()
        << static_cast<int>(broker()->timer().getCpuTime())
        << CoinMessageEol;
    }
    else {
      std::stringstream lb_ss;
      lb_ss << std::setw(14) << std::left << std::scientific << lb;
      broker()->setNumNodeLog(broker()->getNumNodeLog()+1);
      dcoMessageHandler_->message(DISCO_NODE_LOG_NO_SOL, *dcoMessages_)
        << num_processed
        << num_left
        << lb_ss.str().c_str()
        << static_cast<int>(broker()->timer().getCpuTime())
        << CoinMessageEol;
    }
  }
  dcoMessageHandler_->setPrefix(1);
}

/// This is called at the end of the AlpsKnowledgeBroker::rootSearch
/// Prints solution statistics
void DcoModel::modelLog() {
  if (broker_->getProcType() == AlpsProcessTypeSerial) {
#if defined(__OA__) || defined(__COLA__)
    // report solver iterations

    // notes(aykut) This will overflow if number of iterations is large, but I
    // have no option since CoinMessageHandler does not overload operator <<
    // for long long int.
    int num_iter = static_cast<int> (numRelaxIterations_);
    dcoMessageHandler_->message(DISCO_SOLVER_ITERATIONS, *dcoMessages_)
      << num_iter
      << CoinMessageEol;
#endif
    // report cut generator statistics
    std::map<DcoConstraintType, DcoConGenerator*>::iterator it;
    for (it=conGenerators_.begin(); it != conGenerators_.end(); ++it) {
      DcoConGenerator * curr = it->second;
      if (curr->stats().numCalls() > 0) {
        dcoMessageHandler_->message(DISCO_CUT_STATS_FINAL,
                                        *dcoMessages_)
          << curr->name()
          << curr->stats().numCalls()
          << curr->stats().numConsGenerated()
          << curr->stats().numConsUsed()
          << curr->stats().time()
          << curr->strategy()
          << CoinMessageEol;
      }
    }
    // report heuristic statistics
    for (unsigned int k=0; k<heuristics_.size(); ++k) {
      if (heuristics(k)->stats().numCalls() > 0) {
        dcoMessageHandler_->message(DISCO_HEUR_STATS_FINAL,
                                    *dcoMessages_)
          << heuristics(k)->name()
          << heuristics(k)->stats().numCalls()
          << heuristics(k)->stats().numSolutions()
          << heuristics(k)->stats().time()
          << heuristics(k)->strategy()
          << CoinMessageEol;
      }
    }
  }
  else if (broker_->getProcType()==AlpsProcessTypeMaster) {
    dcoMessageHandler_->message(0, "Dco",
                                "Don't know how to log in parallel mode.",
                                'G', DISCO_DLOG_MPI)
      << CoinMessageEol;
  }
}

void DcoModel::reportFeasibility() {
  // return if there is no solution to report
  if (broker_->getNumKnowledges(AlpsKnowledgeTypeSolution)==0) {
    return;
  }
  if (broker_->getProcType()!=AlpsProcessTypeSerial and
      broker_->getProcType()!=AlpsProcessTypeMaster) {
    return;
  }
  // get solution
  DcoSolution * dcosol =
    dynamic_cast<DcoSolution*>
    (broker_->getBestKnowledge(AlpsKnowledgeTypeSolution).first);
  // conic constraints
  double const * sol = dcosol->getValues();

  // integrality
  dcoMessageHandler_->message(0, "Dco", "Integrality Report",
                              'G', DISCO_DLOG_PROCESS)
    << CoinMessageEol;
  std::stringstream msg;
  for (int i=0; i<numRelaxedCols_; ++i) {
    msg << "x["
        << relaxedCols_[i]
        << "] = "
        << sol[relaxedCols_[i]];
    dcoMessageHandler_->message(0, "Dco", msg.str().c_str(),
                                'G', DISCO_DLOG_PROCESS)
      << CoinMessageEol;
    msg.str(std::string());
  }

  // conic constraints
  dcoMessageHandler_->message(0, "Dco", "Conic Feasibility",
                              'G', DISCO_DLOG_PROCESS)
    << CoinMessageEol;
  for (int i=numLinearRows_; i<numLinearRows_+numConicRows_; ++i) {
    DcoConicConstraint * con =
      dynamic_cast<DcoConicConstraint*> (constraints_[i]);
    int const * members = con->coneMembers();
    DcoLorentzConeType type = con->coneType();
    int size = con->coneSize();

    double * values = new double[size];
    for (int j=0; j<size; ++j) {
      values[j] = sol[members[j]];
    }
    double term1 = 0.0;
    double term2 = 0.0;
    if (type==DcoLorentzCone) {
      term1 = values[0];
      term2 = std::inner_product(values+1, values+size, values+1, 0.0);
      term2 = sqrt(term2);
    }
    else if (type==DcoRotatedLorentzCone) {
      term1 = 2.0*sol[members[0]]*sol[members[1]];
      term2 = std::inner_product(values+2, values+size, values+2, 0.0);
    }
    else {
      dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE,
                                  *dcoMessages_)
        << type << CoinMessageEol;
    }
    msg << "Cone "
        << i - numLinearRows_
        << " "
        << term1-term2;
    dcoMessageHandler_->message(0, "Dco", msg.str().c_str(),
                                'G', DISCO_DLOG_PROCESS)
      << CoinMessageEol;
    msg.str(std::string());
    delete[] values;
  }

  // print maximimum violation for integrality
  {
    double max_violation = 0.0;
    for (int i=0; i<numRelaxedCols_; ++i) {
      double value = sol[relaxedCols_[i]];
      double dist_to_floor = value - floor(value);
      double dist_to_ceil = ceil(value) - value;
      double viol = std::min(dist_to_floor, dist_to_ceil);
      if (viol > max_violation) {
        max_violation = viol;
      }
    }
    dcoMessageHandler_->message(DISCO_SOL_INT_FEAS_REPORT,
                                *dcoMessages_)
      << max_violation << CoinMessageEol;
  }

  // print maximimum violation for conic constraints
  {
    double max_violation = 0.0;
    for (int i=numLinearRows_; i<numLinearRows_+numConicRows_; ++i) {
      DcoConicConstraint * con =
        dynamic_cast<DcoConicConstraint*> (constraints_[i]);
      int const * members = con->coneMembers();
      DcoLorentzConeType type = con->coneType();
      int size = con->coneSize();
      double * values = new double[size];
      for (int j=0; j<size; ++j) {
        values[j] = sol[members[j]];
      }
      double term1 = 0.0;
      double term2 = 0.0;
      if (type==DcoLorentzCone) {
        term1 = values[0];
        term2 = std::inner_product(values+1, values+size, values+1, 0.0);
        term2 = sqrt(term2);
      }
      else if (type==DcoRotatedLorentzCone) {
        term1 = 2.0*sol[members[0]]*sol[members[1]];
        term2 = std::inner_product(values+2, values+size, values+2, 0.0);
      }
      else {
        dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE,
                                    *dcoMessages_)
          << type << CoinMessageEol;
      }
      double viol = term2 - term1;
      if (viol > max_violation) {
        max_violation = viol;
      }
      delete[] values;
    }
    dcoMessageHandler_->message(DISCO_SOL_CONE_FEAS_REPORT,
                                *dcoMessages_)
      << max_violation << CoinMessageEol;
  }
}

/// The method that encodes this instance of model into the given
/// #AlpsEncoded object.
AlpsReturnStatus DcoModel::encode(AlpsEncoded * encoded) const {
  AlpsReturnStatus status;
  // encode Alps parts
  status = AlpsModel::encode(encoded);
  if (status!=AlpsReturnStatusOk) {
    dcoMessageHandler_->message(DISCO_UNEXPECTED_ENCODE_STATUS, *dcoMessages_)
      << __FILE__ << __LINE__ << CoinMessageEol;
  }
  // encode problem name
  int probname_size = static_cast<int>(problemName_.size());
  encoded->writeRep(probname_size);
  encoded->writeRep(problemName_.c_str(), probname_size);
  // encode number of constraints
  encoded->writeRep(numCols_);
  encoded->writeRep(colLB_, numCols_);
  encoded->writeRep(colUB_, numCols_);
  encoded->writeRep(numLinearRows_);
  encoded->writeRep(numConicRows_);
  encoded->writeRep(rowLB_, numRows_);
  encoded->writeRep(rowUB_, numRows_);
  encoded->writeRep(objSense_);
  encoded->writeRep(objCoef_, numCols_);
  encoded->writeRep(numIntegerCols_);
  encoded->writeRep(integerCols_, numIntegerCols_);
  encoded->writeRep(isInteger_, numCols_);
  // encode cone info
  if (numConicRows_) {
    encoded->writeRep(coneStart_, numConicRows_+1);
    encoded->writeRep(coneType_, numConicRows_);
    encoded->writeRep(coneMembers_, coneStart_[numConicRows_]);
  }
  // encode matrix
  encoded->writeRep(matrix_->getNumElements());
  encoded->writeRep(matrix_->getVectorStarts(), numLinearRows_);
  encoded->writeRep(matrix_->getVectorLengths(), numLinearRows_);
  encoded->writeRep(matrix_->getIndices(), matrix_->getNumElements());
  encoded->writeRep(matrix_->getElements(), matrix_->getNumElements());
  encoded->writeRep(initOAcuts_);
  // encode parameters
  dcoPar_->pack(*encoded);

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "Proc[" << broker_->getProcRank() << "]"
            << " model " << this << " encoded." << std::endl;
  dcoMessageHandler_->message(0, "Dco", debug_msg.str().c_str(),
                              'G', DISCO_DLOG_MPI)
    << CoinMessageEol;
  // end of debug stuff

  return status;
}

/// The method that decodes the given #AlpsEncoded object into a new #DcoModel
/// instance and returns a pointer to it.
AlpsKnowledge * DcoModel::decode(AlpsEncoded & encoded) const {
  std::cerr << "not implemented yet." << std::endl;
  return NULL;
}

AlpsReturnStatus DcoModel::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status;
  // decode Alps parts
  status = AlpsModel::decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    dcoMessageHandler_->message(DISCO_UNEXPECTED_DECODE_STATUS, *dcoMessages_)
      << __FILE__ << __LINE__ << CoinMessageEol;
  }
  // decode problem name
  int probname_size = 0;
  char * probname = NULL;
  encoded.readRep(probname_size);
  encoded.readRep(probname, probname_size);
  problemName_ = probname;
  delete[] probname;
  // decode rest
  encoded.readRep(numCols_);
  encoded.readRep(colLB_, numCols_);
  encoded.readRep(colUB_, numCols_);
  encoded.readRep(numLinearRows_);
  encoded.readRep(numConicRows_);
  numRows_ = numLinearRows_+numConicRows_;
  encoded.readRep(rowLB_, numRows_);
  encoded.readRep(rowUB_, numRows_);
  encoded.readRep(objSense_);
  encoded.readRep(objCoef_, numCols_);
  encoded.readRep(numIntegerCols_);
  encoded.readRep(integerCols_, numIntegerCols_);
  encoded.readRep(isInteger_, numCols_);
  // decode cone info
  if (numConicRows_) {
    int cone_start_size;
    encoded.readRep(coneStart_, cone_start_size);
    assert(cone_start_size==numConicRows_+1);
    encoded.readRep(coneType_, numConicRows_);
    encoded.readRep(coneMembers_, coneStart_[numConicRows_]);
  }
  // decode matrix
  int num_elem;
  int * starts;
  int * lengths;
  int * indices;
  double * elements;
  encoded.readRep(num_elem);
  encoded.readRep(starts, numLinearRows_);
  encoded.readRep(lengths, numLinearRows_);
  encoded.readRep(indices, num_elem);
  encoded.readRep(elements, num_elem);
  matrix_ = new CoinPackedMatrix(false, numCols_, numLinearRows_, num_elem,
                                 elements, indices, starts, lengths, 0.0, 0.0);
  encoded.readRep(initOAcuts_);
  delete[] starts;
  delete[] lengths;
  delete[] indices;
  delete[] elements;
  dcoPar_->unpack(encoded);

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "Proc[" << broker_->getProcRank() << "]"
            << " model decoded into " << this << "." << std::endl;
  dcoMessageHandler_->message(0, "Dco", debug_msg.str().c_str(),
                              'G', DISCO_DLOG_MPI)
    << CoinMessageEol;
  // end of debug stuff

  return status;
}

void DcoModel::addNumRelaxIterations() {
  numRelaxIterations_ += solver_->getIterationCount();
}
