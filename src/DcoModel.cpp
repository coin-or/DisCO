#include <CoinMpsIO.hpp>

#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoTreeNode.hpp"
#include "DcoNodeDesc.hpp"
#include "DcoVariable.hpp"
#include "DcoLinearConstraint.hpp"
#include "DcoConicConstraint.hpp"

DcoModel::DcoModel() {
  // set parameters
  dcoPar_ = new DcoParams();
  dcoMessageHandler_ = new CoinMessageHandler();
#ifdef DISCO_DEBUG
  dcoMessageHandler_->setLogLevel(4);
#endif
  dcoMessages_ = new DcoMessage();
  // initialize other fields
  colLB_ = NULL;
  colUB_ = NULL;
  rowLB_ = NULL;
  rowUB_ = NULL;
  numCols_ = 0;
  numRows_ = 0;
  numElems_ = 0;
  matrix_ = NULL;
  objSense_ = 0.0;
  objCoef_ = NULL;
  numCones_ = 0;
  numIntegerCols_ = 0;
  intColIndices_ = NULL;
  numSolutions_ = 0;
  // todo(aykut) this might create problems.
  incObjValue_ = 0.0;
  incumbentSol_ = NULL;
  // todo(aykut) following two might create problems.
  cutoff_ = 0.0;
  cutoffInc_ = 0.0;
  branchStrategy_ = NULL;
  rampUpBranchStrategy_ = NULL;
  activeNode_ = NULL;
  numNodes_ = 0;
  /** Number of lp(Simplex) iterations. */
  numIterations_ = 0;;
  /** Average number of lp iterations to solve a subproblem. */
  aveIterations_ = 0;
  /** The number of passes during bounding procedure.*/
  boundingPass_ = 0;
  /** Current relative optimal gap. */
  currRelGap_ = 1e5;
  /** Current absolute optimal gap. */
  currAbsGap_ = 1e5;
  // todo(aykut) why are we keeping this in the DcoModel class.
  // It seems it is ony used in DcoTreeNode::installSubProblem()
  /** Temporary store old cuts at a node when installing a node. */
  oldConstraints_ = 0;
  /** The memory size allocated for oldConstraints_. */
  oldConstraintsSize_ = 0;
  /** Number of old constraints. */
  numOldConstraints_ = 0;
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
  if (incumbentSol_) {
    delete[] incumbentSol_;
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
  for (int i=0; i<numOldConstraints_; ++i) {
    delete oldConstraints_[i];
  }
  delete[] oldConstraints_;
  if (dcoMessageHandler_) {
    delete dcoMessageHandler_;
  }
  if (dcoMessages_) {
    delete dcoMessages_;
  }
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

// read conic mps
void DcoModel::readInstance(char const * dataFile) {
  // get log level parameters
  int dcoLogLevel =  dcoPar_->entry(DcoParams::logLevel);
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
  reader->messageHandler()->setLogLevel(dcoLogLevel);
  reader->readMps(dataFile, "");
  numCols_ = reader->getNumCols();
  numRows_ = reader->getNumRows();
  numElems_ = reader->getNumElements();
  matrix_ = new CoinPackedMatrix(*(reader->getMatrixByCol()));
  // == allocate variable bounds
  colLB_ = new double [numCols_];
  colUB_ = new double [numCols_];
  // == allocate row bounds
  rowLB_ = new double [numRows_];
  rowUB_ = new double [numRows_];
  // == copy bounds
  std::copy(reader->getColLower(), reader->getColLower()+numCols_, colLB_);
  std::copy(reader->getColUpper(), reader->getColUpper()+numCols_, colUB_);
  std::copy(reader->getRowLower(), reader->getRowLower()+numRows_, rowLB_);
  std::copy(reader->getRowUpper(), reader->getRowUpper()+numRows_, rowUB_);
  // == set objective sense
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

// Add variables to *this.
void DcoModel::readAddVariables(CoinMpsIO * reader) {
  // get variable integrality constraints
  char const * is_integer = reader->integerColumns();
  DcoIntegralityType * i_type = new DcoIntegralityType[numCols_];
  for (int i=0; i<numCols_; ++i) {
    if (is_integer[i]) {
      i_type[i] = DcoIntegralityTypeInt;
    }
    else {
      i_type[i] = DcoIntegralityTypeCont;
    }
  }
  // add variables
  BcpsVariable ** variables = new BcpsVariable*[numCols_];
  for (int i=0; i<numCols_; ++i) {
    variables[i] = new DcoVariable(colLB_[i], colUB_[i], colLB_[i],
				   colLB_[i], i_type[i]);
  }
  setVariables(variables, numCols_);
  delete[] i_type;
  for (int i=0; i<numCols_; ++i) {
    delete variables[i];
  }
  delete[] variables;
}

void DcoModel::readAddLinearConstraints(CoinMpsIO * reader) {
  // == add constraints to *this
  BcpsConstraint ** constraints = new BcpsConstraint*[numRows_+numCones_];
  CoinPackedMatrix const * matrix = reader->getMatrixByRow();
  int const * indices = matrix->getIndices();
  double const * values = matrix->getElements();
  int const * lengths = matrix->getVectorLengths();
  int const * starts = matrix->getVectorStarts();
  for (int i=0; i<numRows_; ++i) {
    constraints[i] = new DcoLinearConstraint(lengths[i], indices+starts[i],
					     values+starts[i], rowLB_[i],
					     rowUB_[i]);
  }
  setConstraints(constraints, numRows_+numCones_);
  for (int i=0; i<numRows_+numCones_; ++i) {
    delete constraints[i];
  }
  delete[] constraints;
}

// read conic part and add conic constraints to the model.
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
  numCones_ = nOfCones;
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
    //				  coneMembers+coneStart[i]);
#ifndef __OA__
    OsiLorentzConeType osi_type;
    if (coneType[i]==1) {
      type = OSI_QUAD;
    }
    else if (coneType[i]==2) {
      type = OSI_RQUAD;
    }
    solver_->addConicConstraint(type, coneStart[i+1]-coneStart[i],
				coneMembers+coneStart[i]);
#endif
  }
  delete[] coneStart;
  delete[] coneMembers;
  delete[] coneType;
}

void DcoModel::preprocess() {
#ifdef __OA__
  approximateCones();
  // update row information
  numRows_ = solver_->getNumRows();
  // update num elements
  numElems_ = solver_->getNumElements();
  // update row bounds
  delete[] rowLB_;
  delete[] rowUB_;
  rowLB_ = new double [numRows_];
  rowUB_ = new double [numRows_];
  std::copy(solver_->getRowLower(),
	    solver_->getRowLower()+numRows_,
	    rowLB_);
  std::copy(solver_->getRowUpper(),
	    solver_->getRowUpper()+numRows_,
	    rowUB_);
#endif
}

void DcoModel::approximateCones() {
}


AlpsTreeNode * DcoModel::createRoot() {
  DcoTreeNode * root = new DcoTreeNode();
  DcoNodeDesc * desc = new DcoNodeDesc(this);
  root->setDesc(desc);
  std::vector<BcpsVariable *> cols = getVariables();
  std::vector<BcpsConstraint *> rows = getConstraints();
  int numCols = static_cast<int> (cols.size());
  int numRows = static_cast<int> (rows.size());
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
