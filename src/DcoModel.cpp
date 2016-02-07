#include <CoinMpsIO.hpp>

#include "DcoModel.hpp"
#include "DcoMessage.hpp"
#include "DcoTreeNode.hpp"
#include "DcoNodeDesc.hpp"

DcoModel::DcoModel() {
  // set parameters
  dcoPar_ = new DcoParams();
  dcoMessageHandler_ = new CoinMessageHandler();
  dcoMessages_ = new DcoMessage();
}

DcoModel::~DcoModel() {
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
  int reader_return;
  CoinMpsIO * reader = new CoinMpsIO;
  reader->messageHandler()->setLogLevel(dcoLogLevel);
  reader_return = reader->readMps(dataFile, "");
  numCols_ = reader->getNumCols();
  numRows_ = reader->getNumRows();
  numElems_ = reader->getNumElements();
  CoinPackedMatrix * matrix = new CoinPackedMatrix(*(reader->getMatrixByCol()));
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
  // == set variable type
  colType_ = new char[numCols_];
  for(int i = 0; i<numCols_; ++i) {
    if (reader->isContinuous(i)) {
      colType_[i] = 0;
    }
    else {
      if (colLB_[i] == 0 && colUB_[i] == 1.0) {
	colType_[i] = 1;
      }
      else {
	colType_[i] = 2;
      }
    }
  }
  // == load data to solver
  solver_->loadProblem(*matrix, colLB_, colUB_, objCoef_,
		       rowLB_, rowUB_);
  // == read conic part
  int nOfCones = 0;
  int * coneStart = NULL;
  int * coneIdx = NULL;
  int * coneType = NULL;
  reader_return = reader->readConicMps(NULL, coneStart, coneIdx,
				       coneType, nOfCones);
  // == when there is no conic section status is -3.
  if (reader_return==-3) {
    dcoMessageHandler_->message(DISCO_READ_NOCONES,
				*dcoMessages_);
  }
  else if (reader_return!=0) {
    dcoMessageHandler_->message(DISCO_READ_MPSERROR,
				*dcoMessages_) << reader_return
					      << CoinMessageEol;
  }
  // == allocate memory for cone data members of the class
  numCones_ = nOfCones;
  cone_ = new OsiLorentzCone*[nOfCones];
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
    OsiLorentzConeType type;
    if (coneType[i]==1) {
      type = OSI_QUAD;
    }
    else if (coneType[i]==2) {
      type = OSI_RQUAD;
    }
    cone_[i] = new OsiLorentzCone(type, coneStart[i+1]-coneStart[i],
				  coneIdx+coneStart[i]);
#ifndef __OA__
    solver_->addConicConstraint(type, coneStart[i+1]-coneStart[i],
				coneIdx+coneStart[i]);
#endif
  }
  if (nOfCones) {
    dcoMessageHandler_->message(DISCO_READ_CONESTATS1,
				*dcoMessages_) << nOfCones
					      << CoinMessageEol;
    for (int i=0; i<nOfCones; ++i) {
      dcoMessageHandler_->message(DISCO_READ_CONESTATS2,
				  *dcoMessages_)
	<< i
	<< coneStart[i+1]-coneStart[i]
	<< cone_[i]->type()
	<< CoinMessageEol;
    }
  }
  delete [] coneStart;
  delete [] coneIdx;
  delete [] coneType;
  delete reader;
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
  // create columns and rows of the problem
  createObjects();
}

void DcoModel::approximateCones() {
}

void DcoModel::createObjects() {

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
