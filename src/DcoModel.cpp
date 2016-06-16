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
#include "DcoConGenerator.hpp"
#include "DcoLinearConGenerator.hpp"
#include "DcoConicConGenerator.hpp"
#include "DcoSolution.hpp"
#include "DcoPresolve.hpp"
#include "DcoHeuristic.hpp"
#include "DcoHeurRounding.hpp"

// MILP cuts
#include <CglCutGenerator.hpp>
#include <CglProbing.hpp>
#include <CglClique.hpp>
#include <CglOddHole.hpp>
#include <CglFlowCover.hpp>
#include <CglKnapsackCover.hpp>
#include <CglMixedIntegerRounding2.hpp>
#include <CglGomory.hpp>
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

DcoModel::DcoModel() {
  solver_ = NULL;
  colLB_ = NULL;
  colUB_ = NULL;
  rowLB_ = NULL;
  rowUB_ = NULL;
  numCols_ = 0;
  numRows_ = 0;
  numLinearRows_ = 0;
  numConicRows_ = 0;
  matrix_ = NULL;
  objSense_ = 0.0;
  objCoef_ = NULL;
  numIntegerCols_ = 0;
  integerCols_ = NULL;

  currRelGap_ = 1e5;
  currAbsGap_ = 1e5;
  bestQuality_ = COIN_DBL_MAX;
  activeNode_ = NULL;
  dcoPar_ = new DcoParams();
  numNodes_ = 0;
  numIterations_ = 0;;
  aveIterations_ = 0;
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
}

DcoModel::~DcoModel() {
  // solver_ is freed in main function.
  if (matrix_) {
    delete matrix_;
  }
  if (colLB_) {
    delete[] colLB_;
  }
  if (colUB_) {
    delete[] colUB_;
  }
  if (rowLB_) {
    delete[] rowLB_;
  }
  if (rowUB_) {
    delete[] rowUB_;
  }
  if (objCoef_) {
    delete[] objCoef_;
  }
  if (branchStrategy_) {
    delete branchStrategy_;
  }
  if (rampUpBranchStrategy_) {
    delete rampUpBranchStrategy_;
  }
  if (activeNode_) {
    delete activeNode_;
  }
  if (dcoPar_) {
    delete dcoPar_;
  }
  if (dcoMessageHandler_) {
    delete dcoMessageHandler_;
  }
  if (dcoMessages_) {
    delete dcoMessages_;
  }
  if (relaxedCols_) {
    delete[] relaxedCols_;
  }
  if (relaxedRows_) {
    delete[] relaxedRows_;
  }
  for (std::vector<DcoConGenerator*>::iterator it=conGenerators_.begin();
       it!=conGenerators_.end(); ++it) {
    delete *it;
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

void DcoModel::readInstance(char const * dataFile) {
  // get input file name
  std::string input_file(dataFile);
  std::string base_name = input_file.substr(0, input_file.rfind('.'));
  std::string extension = input_file.substr(input_file.rfind('.')+1);
  if (extension.compare("mps")) {
    dcoMessageHandler_->message(DISCO_READ_MPSFILEONLY,
                                *dcoMessages_) << CoinMessageEol;
  }

  // read mps file
  CoinMpsIO * reader = new CoinMpsIO;

  // get log level parameters
  int dcoLogLevel =  dcoPar_->entry(DcoParams::logLevel);
  reader->messageHandler()->setLogLevel(dcoLogLevel);
  reader->readMps(dataFile, "");
  numCols_ = reader->getNumCols();
  numLinearRows_ = reader->getNumRows();
  matrix_ = new CoinPackedMatrix(*(reader->getMatrixByCol()));

  // == allocate variable bounds
  colLB_ = new double [numCols_];
  colUB_ = new double [numCols_];

  // == allocate row bounds
  rowLB_ = new double [numLinearRows_];
  rowUB_ = new double [numLinearRows_];

  // == copy bounds
  std::copy(reader->getColLower(), reader->getColLower()+numCols_, colLB_);
  std::copy(reader->getColUpper(), reader->getColUpper()+numCols_, colUB_);
  std::copy(reader->getRowLower(), reader->getRowLower()+numLinearRows_, rowLB_);
  std::copy(reader->getRowUpper(), reader->getRowUpper()+numLinearRows_, rowUB_);

  // == set objective sense
  // todo(aykut) we should ask reader about the objective sense
  objSense_ = dcoPar_->entry(DcoParams::objSense);

  // == allocate objective coefficients
  objCoef_ = new double [numCols_];
  double const * reader_obj = reader->getObjCoefficients();
  if (objSense_ > 0.0) {
    std::copy(reader_obj, reader_obj+numCols_, objCoef_);
  }
  else {
    for (int i = 0; i<numCols_; ++i) {
      objCoef_[i] = -reader_obj[i];
    }
  }
  // == load data to solver
  solver_->loadProblem(*matrix_, colLB_, colUB_, objCoef_,
                       rowLB_, rowUB_);
  // == add variables to the model
  readAddVariables(reader);
  // == add linear constraints to the model
  readAddLinearConstraints(reader);
  // == add conic constraints to the model
  readAddConicConstraints(reader);
  delete reader;
}


void DcoModel::readParameters(const int argnum,
                              const char * const * arglist) {
  AlpsPar()->readFromArglist(argnum, arglist);
  dcoPar_->readFromArglist(argnum, arglist);
}

// Add variables to *this.
void DcoModel::readAddVariables(CoinMpsIO * reader) {
  // get variable integrality constraints
  numIntegerCols_ = 0;
  DcoIntegralityType * i_type = new DcoIntegralityType[numCols_];

  // may return NULL, especially when there are no integer variables.
  char const * is_integer = reader->integerColumns();
  if (is_integer!=NULL) {
    for (int i=0; i<numCols_; ++i) {
      if (is_integer[i]) {
        i_type[i] = DcoIntegralityTypeInt;
        numIntegerCols_++;
      }
      else {
        i_type[i] = DcoIntegralityTypeCont;
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
      if (reader->isContinuous(i)) {
        i_type[i] = DcoIntegralityTypeCont;
      }
      else {
        i_type[i] = DcoIntegralityTypeInt;
        numIntegerCols_++;
      }
    }
  }
  // add variables
  BcpsVariable ** variables = new BcpsVariable*[numCols_];
  for (int i=0; i<numCols_; ++i) {
    variables[i] = new DcoVariable(i, colLB_[i], colUB_[i], colLB_[i],
                                   colUB_[i], i_type[i]);
  }
  setVariables(variables, numCols_);
  // variables[i] are now owned by BcpsModel, do not free them.
  delete[] variables;
  delete[] i_type;
}

void DcoModel::readAddLinearConstraints(CoinMpsIO * reader) {
  // == add constraints to *this
  BcpsConstraint ** constraints = new BcpsConstraint*[numLinearRows_];
  CoinPackedMatrix const * matrix = reader->getMatrixByRow();
  int const * indices = matrix->getIndices();
  double const * values = matrix->getElements();
  int const * lengths = matrix->getVectorLengths();
  int const * starts = matrix->getVectorStarts();
  for (int i=0; i<numLinearRows_; ++i) {
    constraints[i] = new DcoLinearConstraint(lengths[i], indices+starts[i],
                                             values+starts[i], rowLB_[i],
                                             rowUB_[i]);
  }
  setConstraints(constraints, numLinearRows_);
  // constraints[i] are owned by BcpsModel. Do not free them here.
  delete[] constraints;
}

void DcoModel::readAddConicConstraints(CoinMpsIO * reader) {
  int nOfCones = 0;
  int * coneStart = NULL;
  int * coneMembers = NULL;
  int * coneType = NULL;
  int reader_return = reader->readConicMps(NULL, coneStart, coneMembers,
                                       coneType, nOfCones);
  // when there is no conic section status is -3.
  if (reader_return==-3) {
    dcoMessageHandler_->message(DISCO_READ_NOCONES,
                                *dcoMessages_);
  }
  else if (reader_return!=0) {
    dcoMessageHandler_->message(DISCO_READ_MPSERROR,
                                *dcoMessages_) << reader_return
                                              << CoinMessageEol;
  }
  // store number of cones in the problem
  numConicRows_ = nOfCones;
  // log cone information messages
  if (nOfCones) {
    dcoMessageHandler_->message(DISCO_READ_CONESTATS1,
                                *dcoMessages_) << nOfCones
                                               << CoinMessageEol;
    for (int i=0; i<nOfCones; ++i) {
      dcoMessageHandler_->message(DISCO_READ_CONESTATS2,
                                  *dcoMessages_)
        << i
        << coneStart[i+1] - coneStart[i]
        << coneType[i]
        << CoinMessageEol;
    }
  }
  // iterate over cones and add them to the model
  for (int i=0; i<nOfCones; ++i) {
    if (coneType[i]!=1 and coneType[i]!=2) {
      dcoMessageHandler_->message(DISCO_READ_CONEERROR,
                                  *dcoMessages_) << CoinMessageEol;
    }
    int num_members = coneStart[i+1]-coneStart[i];
    if (coneType[i]==2 and num_members<3) {
      dcoMessageHandler_->message(DISCO_READ_ROTATEDCONESIZE,
                                  *dcoMessages_) << CoinMessageEol;
    }
    DcoLorentzConeType type;
    if (coneType[i]==1) {
      type = DcoLorentzCone;
    }
    else if (coneType[i]==2) {
      type = DcoRotatedLorentzCone;
    }
    addConstraint(new DcoConicConstraint(type, num_members,
                                         coneMembers+coneStart[i]));
    // cone_[i] = new OsiLorentzCone(type, num_members,
    //                            coneMembers+coneStart[i]);
#ifndef __OA__
    OsiLorentzConeType osi_type;
    if (coneType[i]==1) {
      osi_type = OSI_QUAD;
    }
    else if (coneType[i]==2) {
      osi_type = OSI_RQUAD;
    }
    solver_->addConicConstraint(osi_type, coneStart[i+1]-coneStart[i],
                                coneMembers+coneStart[i]);
#endif
  }
  delete[] coneStart;
  delete[] coneMembers;
  delete[] coneType;
  // update rowLB_ and rowUB_, add conic row bounds
  double * temp = rowLB_;
  rowLB_ = new double[numLinearRows_+numConicRows_];
  std::copy(temp, temp+numLinearRows_, rowLB_);
  std::fill_n(rowLB_+numLinearRows_, numConicRows_, 0.0);
  delete[] temp;
  temp = rowUB_;
  rowUB_ = new double[numLinearRows_+numConicRows_];
  std::copy(temp, temp+numLinearRows_, rowUB_);
  std::fill_n(rowUB_+numLinearRows_, numConicRows_, DISCO_INFINITY);
  delete[] temp;
}

/** Write out parameters. */
void
DcoModel::writeParameters(std::ostream& outstream) const {
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
  // set message levels
  setMessageLevel();


  // some bounds are improved if updated is true

  // notes(aykut) this can only improve the bounds if the leading variable is
  // bounded.  most of the time it is not. Things can get better if we use some
  // other bound improvement first and it does improve the upper bound of
  // leading variables.
  //bool updated = DcoPresolve::improve_bounds(this);

  // write parameters used
  //writeParameters(std::cout);
  // generate ipm cuts first
  // then generate 50 rounds of oa cuts

  // this process will change the following data members, conLB_, conUB_, numRows_, numElems_
  // what about cone objects and their positions? This should also be updated.
  approximateCones();

  // number of linear rows stored in the solver
  int solver_rows = solver_->getNumRows();

  // update row bounds
  // == lower bound
  delete[] rowLB_;
  rowLB_ = new double[solver_rows+numConicRows_];
  std::copy(solver_->getRowLower(), solver_->getRowLower()+solver_rows, rowLB_);
  std::fill_n(rowLB_+solver_rows, numConicRows_, 0.0);
  // == upper bound
  delete[] rowUB_;
  rowUB_ = new double[solver_rows+numConicRows_];
  std::copy(solver_->getRowUpper(), solver_->getRowUpper()+solver_rows, rowUB_);
  std::fill_n(rowUB_+solver_rows, numConicRows_, DISCO_INFINITY);

  // re-order constraints data member.
  // == copy conic constraint pointers
  std::vector<BcpsConstraint*>
    conic_constraints(constraints_.begin()+numLinearRows_, constraints_.end());
  // pop conic constraints
  for (int i=0; i<numConicRows_; ++i) {
    constraints_.pop_back();
  }
  // == add new linear constraints to *this
  CoinPackedMatrix const * matrix = solver_->getMatrixByRow();
  int const * indices = matrix->getIndices();
  double const * values = matrix->getElements();
  int const * lengths = matrix->getVectorLengths();
  int const * starts = matrix->getVectorStarts();
  for (int i=numLinearRows_; i<solver_rows; ++i) {
    BcpsConstraint * curr = new DcoLinearConstraint(lengths[i],
                                                    indices+starts[i],
                                                    values+starts[i],
                                                    rowLB_[i],
                                                    rowUB_[i]);
    constraints_.push_back(curr);
    curr = NULL;
  }
  // == add conic constraints to the end
  for (int i=0; i<numConicRows_; ++i) {
    constraints_.push_back(conic_constraints[i]);
  }
  conic_constraints.clear();

  // update the rest of the row related data members
  numLinearRows_ = solver_rows;
  numRows_ = numLinearRows_ + numConicRows_;
}

void DcoModel::approximateCones() {
#ifdef __OA__
  bool dual_infeasible = false;
  int iter = 0;
  int ipm_iter;
  int oa_iter;
  int num_ipm_cuts = 0;
  int num_oa_cuts = 0;
  // solve problem
  solver_->resolve();
  // get cone data in the required form
  OsiLorentzConeType * coneTypes = new OsiLorentzConeType[numConicRows_];
  int * coneSizes = new int[numConicRows_];
  int ** coneMembers = new int*[numConicRows_];
  for (int i=0; i<numConicRows_; ++i) {
    DcoConicConstraint * con =
      dynamic_cast<DcoConicConstraint*>(constraints_[numLinearRows_+i]);
    DcoLorentzConeType type = con->getType();
    if (type==DcoLorentzCone) {
      coneTypes[i] = OSI_QUAD;
    }
    else if (type==DcoRotatedLorentzCone) {
      coneTypes[i] = OSI_RQUAD;
    }
    else {
      dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE, *dcoMessages_)
        << __FILE__ << __LINE__ << CoinMessageEol;
    }
    coneSizes[i] = con->getSize();
    coneMembers[i] = new int[con->getSize()];
    std::copy(con->getMembers(), con->getMembers()+con->getSize(), coneMembers[i]);
  }
  do {
    // generate cuts
    OsiCuts * ipm_cuts = new OsiCuts();
    OsiCuts * oa_cuts = new OsiCuts();
    CglConicCutGenerator * cg_ipm = new CglConicIPM();
    CglConicCutGenerator * cg_oa =
      new CglConicOA(dcoPar_->entry(DcoParams::coneTol));
    // get cone info
    int largest_cone_size = *std::max_element(coneSizes, coneSizes+numConicRows_);
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
  // add outer apprixmating cuts for 50 rounds
  ipm_iter = iter;
  iter = 0;
  // todo(aykut): parametrize 50
  while(iter<5) {
    OsiCuts * oa_cuts = new OsiCuts();
    CglConicCutGenerator * cg_oa =
      new CglConicOA(dcoPar_->entry(DcoParams::coneTol));
    cg_oa->generateCuts(*solver_, *oa_cuts, numConicRows_, coneTypes,
                        coneSizes, coneMembers, 1);
    int num_cuts = oa_cuts->sizeRowCuts();
    num_oa_cuts += num_cuts;
    if (num_cuts==0) {
      // ifno cuts are produced break early
      break;
    }
    solver_->applyCuts(*oa_cuts);
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
  std::cout << "Linear relaxation objective value " << solver_->getObjValue() << std::endl;
  std::cout << "=================================" << std::endl;
  delete coneTypes;
  delete coneSizes;
  for (int i=0; i<numConicRows_; ++i) {
    delete[] coneMembers[i];
  }
  delete[] coneMembers;
#endif
}

//todo(aykut) why does this return to bool?
// should be fixed in Alps level.
bool DcoModel::setupSelf() {

  // set integer column indices
  std::vector<BcpsVariable*> & cols = getVariables();
  int numCols = getNumCoreVariables();
  integerCols_ = new int[numIntegerCols_];
  for (int i=0, k=0; i<numCols; ++i) {
    BcpsIntegral_t i_type = cols[i]->getIntType();
    if (i_type=='I' or i_type=='B') {
      integerCols_[k] = i;
      k++;
    }
  }
  // set relaxed array for integer columns
  numRelaxedCols_ = numIntegerCols_;
  relaxedCols_ = new int[numRelaxedCols_];
  std::copy(integerCols_, integerCols_+numIntegerCols_,
            relaxedCols_);
#if defined(__OA__)
  // we relax conic constraint too when OA is used.
  // set relaxed array for conic constraints
  numRelaxedRows_ = numConicRows_;
  relaxedRows_ = new int[numRelaxedRows_];
  // todo(aykut) we assume conic rows start after linear rows.
  // if not iterate over rows and determine their type (DcoConstraint::type())
  for (int i=0; i<numRelaxedRows_; ++i) {
    relaxedRows_[i] = numLinearRows_+i;
  }

  // set leading variable lower bounds to 0
  for (int i=numLinearRows_; i<numLinearRows_+numConicRows_; ++i) {
    DcoConicConstraint * con =
      dynamic_cast<DcoConicConstraint*> (constraints_[i]);
    if (con->getType()==DcoLorentzCone) {
      //todo(aykut) assumes variable indices in reader and variable
      // indices in disco are same
      variables_[con->getMembers()[0]]->setLbHard(0.0);
    }
    else if (con->getType()==DcoRotatedLorentzCone) {
      variables_[con->getMembers()[0]]->setLbHard(0.0);
      variables_[con->getMembers()[1]]->setLbHard(0.0);
    }
    else {
      dcoMessageHandler_->message(DISCO_UNKNOWN_CONETYPE,
                                  *dcoMessages_)
        << con->getType() << CoinMessageEol;
    }
  }
#endif

  // set branch strategy
  setBranchingStrategy();

  // add constraint generators
  addConstraintGenerators();

  // add heuristics
  addHeuristics();

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
  getKnowledgeBroker()->messageHandler()->setLogLevel(alps_log_level);
  bcpsMessageHandler()->setLogLevel(alps_log_level);
  dcoMessageHandler_->setLogLevel(dco_log_level);


  // dynamic_cast<OsiClpSolverInterface*>(solver)->setHintParam(OsiDoReducePrint,false,OsiHintDo, 0);
  solver_->setHintParam(OsiDoReducePrint,false,OsiHintDo, 0);
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
    addConGenerator(probing, "Probing", probeStrategy, probeFreq);
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
    addConGenerator(cliqueCut, "Clique", cliqueStrategy, cliqueFreq);
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
    addConGenerator(oldHoleCut, "OddHole", oddHoleStrategy, oddHoleFreq);
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
    addConGenerator(flowGen, "Flow Cover", fCoverStrategy, fCoverFreq);
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
    addConGenerator(knapCut, "Knapsack", knapStrategy, knapFreq);
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
    addConGenerator(mixedGen, "MIR", mirStrategy, mirFreq);
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
    // try larger limit
    gomoryCut->setLimit(300);
    addConGenerator(gomoryCut, "Gomory", gomoryStrategy, gomoryFreq);
  }

  // Add Tow MIR cut generator.
  // Disable forever, not useful.
  twoMirStrategy = DcoCutStrategyNone;
  if (twoMirStrategy != DcoCutStrategyNone) {
    CglTwomir *twoMirCut =  new CglTwomir;
    addConGenerator(twoMirCut, "Two MIR", twoMirStrategy, twoMirFreq);
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
    addConGenerator(ipm_gen, "IPM", ipmStrategy, ipmFreq);
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
    addConGenerator(ipm_int_gen, "IPMint", ipmintStrategy, ipmintFreq);
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
  if (oaStrategy != DcoCutStrategyNone) {
    CglConicCutGenerator * oa_gen =
      new CglConicOA(dcoPar_->entry(DcoParams::coneTol));
    addConGenerator(oa_gen, "OA", oaStrategy, oaFreq);
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
  std::vector<DcoConGenerator*>::iterator it;
  for (it=conGenerators_.begin(); it!=conGenerators_.end(); ++it) {
    DcoCutStrategy curr = (*it)->strategy();
    if (curr==DcoCutStrategyPeriodic) {
      periodic_exists = true;
      break;
    }
    else if (curr==DcoCutStrategyRoot) {
      root_exists = true;
    }
  }
  if (periodic_exists) {
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
                               char const * name,
                               DcoCutStrategy dco_strategy,
                               int frequency) {
  DcoConGenerator * con_gen = new DcoLinearConGenerator(this, cgl_gen, name,
                                                        dco_strategy,
                                                        frequency);
  conGenerators_.push_back(con_gen);
}

/// Add constraint generator.
void DcoModel::addConGenerator(CglConicCutGenerator * cgl_gen,
                               char const * name,
                               DcoCutStrategy dco_strategy,
                               int frequency) {
  DcoConGenerator * con_gen = new DcoConicConGenerator(this, cgl_gen, name,
                                                       dco_strategy,
                                                       frequency);
  conGenerators_.push_back(con_gen);
}

void DcoModel::addHeuristics() {
  // todo(aykut) this function (heuristic adding process) can be improved
  // since global parameters are ignored with this design.

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
  // case DcoBranchingStrategyStrong:
  //   branchStrategy_ = new DcoBranchStrategyStrong(this);
  //   break;
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
  // case DcoBranchingStrategyStrong:
  //   rampUpBranchStrategy_ = new DcoBranchStrategyStrong(this);
  //   break;
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
                                         int & numInfRows) {
  // set stats to 0
  numInfColumns = 0;
  numInfRows = 0;

  // for debug purposes we will keep the largest column and row infeasibility
  double col_inf = 0.0;
  double row_inf = 0.0;

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
      if (col_inf<infeas) {
        col_inf = infeas;
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
      if (row_inf<infeas) {
        row_inf = infeas;
      }
    }
  }

  // debug stuff
  // report largest column and row infeasibilities
  std::stringstream debug_msg;
  debug_msg << "Column infeasibility "
            << col_inf
            << " Row infeasibility "
            << row_inf;
  dcoMessageHandler_->message(0, "Dco", debug_msg.str().c_str(),
                              'G', DISCO_DLOG_PROCESS)
    << CoinMessageEol;
  // end of debug stuff


  // create DcoSolution instance if feasbile
  DcoSolution * dco_sol = 0;
  if (numInfColumns==0 && numInfRows==0) {
    double const * sol = solver()->getColSolution();
    double quality = solver()->getObjValue();
    dco_sol = new DcoSolution(numCols_, sol, quality);

    // debug stuff
    // flush stream
    debug_msg.str(std::string());
    debug_msg << "Solution found. ";
    debug_msg << "Obj value ";
    debug_msg << quality;
    dcoMessageHandler_->message(0, "Dco", debug_msg.str().c_str(),
                                'G', DISCO_DLOG_PROCESS)
      << CoinMessageEol;
    // end of debug stuff

  }
  return dco_sol;
}

// todo(aykut) why does this return int?
// todo(aykut) what if objsense is -1 and problem is maximization?
int DcoModel::storeSolution(DcoSolution * sol) {
  double quality = sol->getQuality();
  // Store in Alps pool, assumes minimization.
  getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
                                     sol,
                                     objSense_ * quality);
  if (quality<bestQuality_) {
    bestQuality_ = quality;
    solver_->setDblParam(OsiDualObjectiveLimit, objSense_*quality);
  }

  // // debug write mps file to disk
  // std::cout << "writing problem to disk..." << std::endl;
  // std::stringstream problem;
  // problem << broker_->getNumNodesProcessed();
  // solver_->writeMps(problem.str().c_str(), "mps", 0.0);
  // // end of problem writing
  return AlpsReturnStatusOk;
}

double DcoModel::bestQuality() {
  return bestQuality_;
  //return getKnowledgeBroker()->getBestQuality();
}

/// This is called at the end of the AlpsKnowledgeBroker::rootSearch
/// Prints solution statistics
void DcoModel::modelLog() {
  if (broker_->getProcType() == AlpsProcessTypeSerial) {
    for (int k=0; k<conGenerators_.size(); ++k) {
      if (conGenerators(k)->stats().numCalls() > 0) {
        dcoMessageHandler_->message(DISCO_CUT_STATS_FINAL,
                                        *dcoMessages_)
          << conGenerators(k)->name()
          << conGenerators(k)->stats().numCalls()
          << conGenerators(k)->stats().numConsGenerated()
          << conGenerators(k)->stats().time()
          << conGenerators(k)->strategy()
          << CoinMessageEol;
      }
    }
    for (int k=0; k<heuristics_.size(); ++k) {
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
  int prefered;
  for (int i=numLinearRows_; i<numLinearRows_+numConicRows_; ++i) {
    DcoConicConstraint * con =
      dynamic_cast<DcoConicConstraint*> (constraints_[i]);
    int const * members = con->getMembers();
    DcoLorentzConeType type = con->getType();
    int size = con->getSize();

    double * values = new double[size];
    for (int i=0; i<size; ++i) {
      values[i] = sol[members[i]];
    }
    double term1;
    double term2;
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
  }
}
