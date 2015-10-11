/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * ALPS is distributed under the Eclipse Public License as part of the       *
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

#include <cstdio>
#include <vector>

#include "CoinFinite.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"

#include "BcpsObject.h"

#include "DcoBranchObjectInt.hpp"
#include "DcoBranchObjectBilevel.hpp"
#include "DcoBranchStrategyMaxInf.hpp"
#include "DcoBranchStrategyPseudo.hpp"
#include "DcoBranchStrategyRel.hpp"
#include "DcoBranchStrategyStrong.hpp"
#include "DcoBranchStrategyBilevel.hpp"

#include "DcoConstraint.hpp"
#include "DcoConicConstraint.hpp"
#include "DcoHeurRound.hpp"
#include "DcoModel.hpp"
#include "DcoObjectInt.hpp"
#include "DcoSolution.hpp"
#include "DcoTreeNode.hpp"
#include "DcoVariable.hpp"

#include "CglConicOA.hpp"
#include "DcoConicConGenerator.hpp"

#define DISCO_MIN_SHARE_CON 5
#define DISCO_MAX_SHARE_CON 25
#define DISCO_MIN_SHARE_VAR 5
#define DISCO_MAX_SHARE_VAR 25

//#############################################################################

// Returns the greatest common denominator of two
// positive integers, a and b, found using Euclid's algorithm
static int gcd(int a, int b)
{
    int remainder = -1;
    // make sure a<=b (will always remain so)
    if(a > b) {
	// Swap a and b
	int temp = a;
	a = b;
	b = temp;
    }
    // if zero then gcd is nonzero (zero may occur in rhs of packed)
    if (!a) {
	if (b) {
	    return b;
	} else {
	    printf("**** gcd given two zeros!!\n");
	    abort();
	}
    }
    while (remainder) {
	remainder = b % a;
	b = a;
	a = remainder;
    }
    return b;
}

//#############################################################################

void
DcoModel::init()
{
    origLpSolver_ = NULL;
    presolvedLpSolver_ = NULL;
    lpSolver_ = NULL;

    numCols_ = 0;
    numRows_ = 0;
    numElems_ = 0;

    colMatrix_ = 0;

    varLB_ = NULL;
    varUB_ = NULL;
    conLB_ = NULL;
    conUB_ = NULL;

    startVarLB_ = NULL;
    startVarUB_ = NULL;
    startConLB_ = NULL;
    startConUB_ = NULL;

    tempVarLBPos_ = NULL;
    tempVarUBPos_ = NULL;
    tempConLBPos_ = NULL;
    tempConUBPos_ = NULL;

    // cone data
    numCoreCones_ = 0;
    coneMembers_ = NULL;
    coneTypes_ = NULL;
    coneSizes_ = NULL;

    objSense_ = 1.0;
    objCoef_ = NULL;

    colType_ = 0;
    numIntObjects_ = 0;
    intColIndices_ = NULL;
    intObjIndices_ = NULL;

    presolve_ = new DcoPresolve();

    numSolutions_ = 0;
    numHeurSolutions_ = 0;

    incumbent_ = NULL;

    activeNode_ = NULL;
    numStrong_ = 0;
    numBranchResolve_ = 0;

    cutoff_ = COIN_DBL_MAX;
    incObjValue_ = COIN_DBL_MAX;
    blisMessageHandler_ = new CoinMessageHandler();
    blisMessages_ = DcoMessage();

    objects_ = NULL;
    numObjects_ = 0;

    numNodes_ = 0;
    numIterations_ = 0;
    aveIterations_ = 0;

    branchStrategy_ = NULL;
    rampUpBranchStrategy_ = NULL;
    priority_ = NULL;
    nodeWeight_ = 1.0;

    isRoot_ = true;
    boundingPass_ = 0;
    DcoPar_ = new DcoParams;

    // Tolerance
    integerTol_ = 1.0e-5;
    optimalRelGap_ = 1.0e-4;
    optimalAbsGap_ = 1.0e-6;
    currRelGap_ = ALPS_OBJ_MAX;
    currAbsGap_ = ALPS_OBJ_MAX;

    // Heuristic
    heurStrategy_ = DcoHeurStrategyAuto;
    heurCallFrequency_ = 1;
    numHeuristics_ = 0;
    heuristics_ = NULL;

    // Cons related
    cutStrategy_ = DcoCutStrategyRoot;
    cutGenerationFrequency_ = 1;
    numCutGenerators_ = 0;
    generators_ = NULL;
    constraintPool_ = NULL;
    constraintPoolReceive_ = NULL;
    constraintPoolSend_ = NULL;
    oldConstraints_ = NULL;
    oldConstraintsSize_ = 0;
    numOldConstraints_ = 0;
    conRandoms_ = NULL;

    sharedObjectMark_ = NULL;

    // conic cut related
    conicCutStrategy_ = DcoConicCutStrategyRoot;
    conicCutGenerationFrequency_ = 1;

    // AT - Begin - Set flags for presolve
    presolved=false;
    problemSetup=false;
    // AT - End
}

//#############################################################################

// Read from file (currently MPS format. TODO: LP format).

void DcoModel::readInstance(const char* dataFile) {
  int j;
  int msgLevel =  AlpsPar_->entry(AlpsParams::msgLevel);
  char ext[5];
  //------------------------------------------------------
  // Read in data from file.
  //------------------------------------------------------
  int last_dot(0);
  for (j = 0;; j++){
    if (dataFile[j] == '\0')
      break;
    if (dataFile[j] == '.') {
      last_dot = j;
    }
  }
  if(last_dot){
    strcpy(ext, dataFile + last_dot + 1);
  }else{
    strcpy(ext, "");
  }

  int rc = 1;
  if (strcmp(ext, "mps") != 0) {
    std::cerr << "Mps files only!" << std::endl;
    throw std::exception();
  }
  CoinMpsIO *reader = new CoinMpsIO;
  reader->messageHandler()->setLogLevel(msgLevel);
  rc = reader->readMps(dataFile, "");
  //------------------------------------------------------
  // Get problem data.
  //------------------------------------------------------
  numCols_ = reader->getNumCols();
  numRows_ = reader->getNumRows();
  numElems_ = reader->getNumElements();

  colMatrix_ = new CoinPackedMatrix();
  *colMatrix_ = *(reader->getMatrixByCol());

  varLB_ = new double [numCols_];
  varUB_ = new double [numCols_];

  conLB_ = new double [numRows_];
  conUB_ = new double [numRows_];

  memcpy(varLB_, reader->getColLower(), sizeof(double) * numCols_);
  memcpy(varUB_, reader->getColUpper(), sizeof(double) * numCols_);

  memcpy(conLB_, reader->getRowLower(), sizeof(double) * numRows_);
  memcpy(conUB_, reader->getRowUpper(), sizeof(double) * numRows_);

  objSense_ = DcoPar_->entry(DcoParams::objSense);

  objCoef_ = new double [numCols_];
  if (objSense_ > 0.0) {
    memcpy(objCoef_, reader->getObjCoefficients(),
	   sizeof(double) * numCols_);
  }
  else {
    const double *readerObj =  reader->getObjCoefficients();
    for (j = 0; j < numCols_; ++j) {
      objCoef_[j] = - readerObj[j];
    }
  }

  //------------------------------------------------------
  // Set colType_
  //------------------------------------------------------

  colType_ = new char [numCols_];

  for(j = 0; j < numCols_; ++j) {
    if (reader->isContinuous(j)) {
      colType_[j] = 'C';
    }
    else {
      if (varLB_[j] == 0 && varUB_[j] == 1.0) {
	colType_[j] = 'B';
      }
      else {
	colType_[j] = 'I';
      }
    }
  }

  // origLpSolver_->readMps(dataFile, "");
  // solver()->readMps(dataFile, "");

  //+++++++++++++++++++++++++++++++++++++++++++++
  // AT Data read. Now we must
  // 1) check if presolve
  // 2) if presolve, create lp : lpSolver_->loadProblem(*colMatrix_,
  //    varLB_, varUB_, 	   objCoef_,  conLB_, conUB_);
  // 3) call presolve
  // 4) post solve the presolvedLP. (not blis postprocess)
  // 5) update data so that we use use the presolved problem: numcol,
  //    var bounds, obj coef and conLB, conUP;
  //+++++++++++++++++++++++++++++++++++++++++++++

  //AT
  std::cout<<" About to presolve for the whole tree"<<std::endl;
  presolveForTheWholeTree();
  //AT
  std::cout<<" Model presolved for the whole tree"<<std::endl;

  //-------------------------------------------------------------
  // Create variables and constraints.
  //-------------------------------------------------------------

  OsiConicSolverInterface * si = origLpSolver_;
  si->loadProblem(*colMatrix_, varLB_, varUB_, objCoef_,
  			     conLB_, conUB_);
  // read conic part
  int nOfCones = 0;
  int * coneStart = NULL;
  int * coneIdx = NULL;
  int * coneType = NULL;
  rc = reader->readConicMps(NULL, coneStart, coneIdx, coneType, nOfCones);
  // when there is no conic section status is -3.
  if (rc==-3) {
    std::cout << "OsiConic: No conic section is mps file." << std::endl;
  }
  else if (rc!=0) {
    std::cerr << "OsiConic: readConicMps returned code " << rc << std::endl;
    throw std::exception();
  }
  // allocate memory for cone data members of the class
  numCoreCones_ = nOfCones;
  coneTypes_ = new int[nOfCones];
  coneSizes_ = new int[nOfCones];
  coneMembers_ = new int*[nOfCones];
  for (int i=0; i<nOfCones; ++i) {
    if (coneType[i]!=1 and coneType[i]!=2) {
      std::cerr << "Invalid cone type!" << std::endl;
      throw std::exception();
    }
    int num_members = coneStart[i+1]-coneStart[i];
    if (coneType[i]==2 and num_members<3) {
      std::cerr<< "Rotated cones should have at least 3 members!" << std::endl;
      throw std::exception();
    }
    // get members
    coneSizes_[i] = num_members;
    coneMembers_[i] = new int[num_members];
    int k=0;
    for (int j=coneStart[i]; j<coneStart[i+1]; ++j) {
      coneMembers_[i][k] = coneIdx[j];
      k++;
    }
    OsiLorentzConeType type;
    if (coneType[i]==1) {
      type = OSI_QUAD;
      coneTypes_[i] = 0;
    }
    else if (coneType[i]==2) {
      type = OSI_RQUAD;
      coneTypes_[i] = 1;
    }
    origLpSolver_->addConicConstraint(type, num_members, coneMembers_[i]);
  }
  // check log level and print ccordingly
  if (nOfCones) {
    printf("Conic section has %d cones\n",nOfCones);
    for (int iCone=0;iCone<nOfCones;iCone++) {
      printf("Cone %d has %d entries (type %d) ",iCone,coneStart[iCone+1]-coneStart[iCone],
	     coneType[iCone]);
      for (int j=coneStart[iCone];j<coneStart[iCone+1];j++)
	printf("%d ",coneIdx[j]);
      printf("\n");
    }
  }
  delete [] coneStart;
  delete [] coneIdx;
  delete [] coneType;
  delete reader;
  createObjects();

  // Generate and add conic cuts
  // if (isRoot) {
  //   OsiConicSolverInterface * si = model->solver();
  //   CglConicGD1 cg(si);
  //   bool equalObj;
  //   CoinRelFltEq eq(0.0001);
  //   int num_cut = 0;
  //   double obj;
  //   do {
  //     // Get current solution value
  //     obj = si->getObjValue();
  //     // Generate and apply cuts
  //     cg.generateAndAddCuts(*si);
  //     //si->writeMps("after_cut");
  //     si->resolve();
  //     equalObj = eq(si->getObjValue(), obj);
  //   } while (!equalObj);
  // }



}

//############################################################################
void DcoModel::createObjects() {
  //------------------------------------------------------
  // Create variables and constraints.
  //------------------------------------------------------
  const double *elements = colMatrix_->getElements();
  const int *indices = colMatrix_->getIndices();
  const int *lengths = colMatrix_->getVectorLengths();
  const CoinBigIndex *starts = colMatrix_->getVectorStarts();
  int beg = 0;
  for (int j=0; j<numCols_; ++j) {
    beg = starts[j];
    DcoVariable * var = new DcoVariable(varLB_[j],
					varUB_[j],
					varLB_[j],
					varUB_[j],
					objCoef_[j],
					lengths[j],
					indices + beg,
					elements + beg);
    var->setObjectIndex(j);
    var->setRepType(BCPS_CORE);
    var->setStatus(BCPS_NONREMOVALBE);
    var->setIntType(colType_[j]);
    variables_.push_back(var);
    // this is released in Bcps deconstructors.
    var = NULL;
  }
  for (int j=0; j<numRows_; ++j) {
    DcoConstraint *con = new DcoConstraint(conLB_[j],
					   conUB_[j],
					   conLB_[j],
					   conUB_[j]);
    con->setObjectIndex(j);
    con->setRepType(BCPS_CORE);
    con->setStatus(BCPS_NONREMOVALBE);
    constraints_.push_back(con);
    // todo(aykut) it seems memory allocated to con is never released.
    con = NULL;
  }
  // Initialize the cone information
  numCoreCones_ = origLpSolver_->getNumCones();
  coneMembers_ = new int*[numCoreCones_];
  coneTypes_ = new int[numCoreCones_];
  coneSizes_ = new int[numCoreCones_];
  for(int j=0; j<numCoreCones_; j++) {
    OsiLorentzConeType type;
    origLpSolver_->getConicConstraint(j, type, coneSizes_[j], coneMembers_[j]);
    if (type==OSI_QUAD) {
      coneTypes_[j] = 0;
    }
    else {
      coneTypes_[j] = 1;
    }
  }
  // add conic constraints as objects
  for (int j=0; j<numCoreCones_; ++j) {
    DcoConicConstraint *con = new DcoConicConstraint(coneSizes_[j],
						     coneTypes_[j],
						     coneMembers_[j]);
    con->setObjectIndex(j+numRows_);
    con->setRepType(BCPS_CORE);
    con->setStatus(BCPS_NONREMOVALBE);
    constraints_.push_back(con);
    // memory for this will get released in ~BcpsConstraint since
    // ~DcoConicConstraint is virtual.
    con = NULL;
  }

  // Set all objects as core by default.
  numCoreVariables_ = numCols_;
  numCoreConstraints_ = numRows_;
}

//#############################################################################

/** For parallel code, only the master calls this function.
 *  1) Set colMatrix_, varLB_, varUB_, conLB_, conUB
 *     numCols_, numRows_
 *  2) Set objCoef_ and objSense_
 *  3) Set colType_ ('C', 'I', or 'B')
 *  4) Set variables_ and constraints_
 *  NOTE: Dco takes over the memory ownship of vars and cons, which
 *        means users must NOT free vars or cons.
 */
void
DcoModel::importModel(std::vector<DcoVariable *> vars,
                       std::vector<DcoConstraint *> cons)
{

    CoinBigIndex i, j;
    int k, size;

    int* varIndices = NULL;
    double *varValues = NULL;

    numCols_ = static_cast<int> (vars.size());
    numRows_ = static_cast<int> (cons.size());

    varLB_ = new double [numCols_];
    varUB_ = new double [numCols_];

    conLB_ = new double [numRows_];
    conUB_ = new double [numRows_];

    objCoef_ = new double [numCols_];
    colType_ = new char [numCols_];

    // Get numElems_ and colType_
    for (i = 0; i < numCols_; ++i) {
        numElems_ += vars[i]->getSize();
        colType_[i] = vars[i]->getIntType();
    }

#if 0
    std::cout << "numCols_ = " << numCols_
              << "; numRows_ = " << numRows_ << std::endl;
#endif

    //------------------------------------------------------
    // Set matrix, bounds
    //------------------------------------------------------

    // For colMatrix_, need free memory
    CoinBigIndex * start = new CoinBigIndex [numCols_+1];
    int* indices = new int[numElems_];
    double* values = new double[numElems_];
    int* length = new int [numCols_];

    // Get varLB_, varUB_, objCoef_, and matrix from variables
    for (numElems_ = 0, i = 0; i < numCols_; ++i) {
        varLB_[i] = vars[i]->getLbHard();
        varUB_[i] = vars[i]->getUbHard();
        objCoef_[i] = vars[i]->getObjCoef();

        start[i] = numElems_;

        varValues = vars[i]->getValues();
        varIndices = vars[i]->getIndices();
        size = vars[i]->getSize();
        for (j = 0; j < size; ++j, ++numElems_){
            indices[numElems_] = varIndices[j];
            values[numElems_] = varValues[j];
        }
    }

    // Last position.
    start[numCols_] = numElems_;

    for (k = 0; k < numCols_; ++k) {
        length[k] = start[(k+1)] - start[k];
    }

    colMatrix_ =  new CoinPackedMatrix(true, // column-majored
                                       numRows_,
                                       numCols_,
                                       numElems_,
                                       values,
                                       indices,
                                       start,
                                       length);


    // Get conLB_ and conUB_
    for (k = 0; k < numRows_; ++k) {
        conLB_[k] = cons[k]->getLbHard();
        conUB_[k] = cons[k]->getUbHard();
    }

    //------------------------------------------------------
    // Set variables_ and constraints_
    //------------------------------------------------------

    for (k = 0; k < numCols_; ++k) {
        variables_.push_back(vars[k]);
    }

    for (k = 0; k < numRows_; ++k) {
        constraints_.push_back(cons[k]);
    }

    //------------------------------------------------------
    // Free memory
    //------------------------------------------------------

    delete [] start;
    delete [] length;
    delete [] indices;
    delete [] values;
}

//############################################################################

/** Read in parameters. */
void
DcoModel::readParameters(const int argnum, const char * const * arglist)
{
    AlpsPar_->readFromArglist(argnum, arglist);

    int msgLevel = AlpsPar_->entry(AlpsParams::msgLevel);
    if (msgLevel > 0) {
	std::cout << "Reading in ALPS parameters ..." << std::endl;
	std::cout << "Reading in DISCO parameters ..." << std::endl;
    }
    bcpsMessageHandler_->setLogLevel(msgLevel);
    blisMessageHandler_->setLogLevel(msgLevel);
    DcoPar_->readFromArglist(argnum, arglist);
}

//##############################################################################

/** Write out parameters. */
void
DcoModel::writeParameters(std::ostream& outstream) const
{
    outstream << "\n================================================"
              <<std::endl;
    outstream << "ALPS Parameters: " << std::endl;
    AlpsPar_->writeToStream(outstream);
    outstream << "\n================================================"
              <<std::endl;
    outstream << "DISCO Parameters: " << std::endl;
    DcoPar_->writeToStream(outstream);
}

//############################################################################

/** Do necessary work to make model usable and presolve.
    Return success or not.
    This function is called when constructing knowledge broker. */
bool
DcoModel::setupSelf()
{
    int j;

    bcpsMessageHandler_->setLogLevel(broker_->getMsgLevel());
    blisMessageHandler_->setLogLevel(broker_->getMsgLevel());

    if (broker_->getMsgLevel() > 0) {
        if (broker_->getProcRank() == broker_->getMasterRank()) {
	   if (strcmp(BCPS_VERSION, "trunk")){
	      std::cout << "Bcps Version: " << BCPS_VERSION << std::endl;
	   }else{
	      std::cout << "Bcps Version: Trunk (unstable) \n";
	   }
#ifdef BCPS_SVN_REV
	   std::cout << "Bcps Revision Number: " << BCPS_SVN_REV << std::endl;
#endif
	   if (strcmp(DISCO_VERSION, "trunk")){
	      std::cout << "Dco Version: " << DISCO_VERSION << std::endl;
	   }else{
	      std::cout << "Dco Version: Trunk (unstable) \n";
	   }
#ifdef DISCO_SVN_REV
	   std::cout << "Dco Revision Number: " << DISCO_SVN_REV << std::endl;
#endif
        }
    }

    //------------------------------------------------------
    // Set numIntObjects_, intColIndices_, intObjectIndices_
    //------------------------------------------------------

    // TODO: now integer, later other objects
    intObjIndices_ = new int [numCols_];
    memset(intObjIndices_, 0, sizeof(int) * numCols_);

    numIntObjects_ = 0;
    intColIndices_ = new int [numCols_];

    for(j = 0; j < numCols_; ++j) {
	if (colType_[j] == 'I' || colType_[j] == 'B') {
            intColIndices_[numIntObjects_++] = j;
	}
    }
    if (numIntObjects_ == 0) {
        if (broker_->getMsgLevel() > 0) {
            bcpsMessageHandler_->message(DISCO_W_LP, blisMessages())
                << CoinMessageEol;
        }
    }

    //------------------------------------------------------
    // Load data to LP solver.
    //------------------------------------------------------

    if (!lpSolver_) {
        // preprocessing causes this check.
        lpSolver_ = origLpSolver_;
    }

    //AT here problem is set for first time for lp solver -> should be moved to read problem??
    if(!presolved)
        lpSolver_->loadProblem(*colMatrix_,
                               varLB_, varUB_,
                               objCoef_,
                               conLB_, conUB_);

    lpSolver_->setObjSense(1.0);
    lpSolver_->setInteger(intColIndices_, numIntObjects_);

    // AT - Begin - solver ready, in principle we could presolve here
    problemSetup=true;
    // AT - End

    //------------------------------------------------------
    // Create integer objects and analyze objective coef.
    //------------------------------------------------------

    createIntgerObjects(true);

    // Do this after loading LP.
    analyzeObjective();

    //------------------------------------------------------
    // Allocate memory.
    //------------------------------------------------------

    startVarLB_ = new double [numCols_];
    startVarUB_ = new double [numCols_];

    startConLB_ = new double [numRows_];
    startConUB_ = new double [numRows_];

    tempVarLBPos_ = new int [numCols_];
    tempVarUBPos_ = new int [numCols_];

    tempConLBPos_ = new int [numRows_];
    tempConUBPos_ = new int [numRows_];

    //------------------------------------------------------
    // Get parameters.
    //------------------------------------------------------

    integerTol_ = DcoPar_->entry(DcoParams::integerTol);
    optimalRelGap_ = DcoPar_->entry(DcoParams::optimalRelGap);
    optimalAbsGap_ = DcoPar_->entry(DcoParams::optimalAbsGap);

    int relibility = DcoPar_->entry(DcoParams::pseudoRelibility);
    cutoff_ =  DcoPar_->entry(DcoParams::cutoff);

    //------------------------------------------------------
    // Modify parameters.
    //------------------------------------------------------

    // Disable Alps message
    // AlpsPar()->setEntry(AlpsParams::msgLevel, 1);

    AlpsPar()->setEntry(AlpsParams::printSystemStatus, 0);

#ifdef DISCO_DEBUG_MORE
    std::string problemName;
    lpSolver_->getStrParam(OsiProbName, problemName);
    printf("DISCO: setupSelf: Problem name - %s\n", problemName.c_str());
    lpSolver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
#endif

    //------------------------------------------------------
    // Set branch strategy.
    //------------------------------------------------------

    int brStrategy = DcoPar_->entry(DcoParams::branchStrategy);

    if (brStrategy == DcoBranchingStrategyMaxInfeasibility) {
        // Max inf
        branchStrategy_ = new DcoBranchStrategyMaxInf(this);
    }
    else if (brStrategy == DcoBranchingStrategyPseudoCost) {
        // Pseudocost
        branchStrategy_ = new DcoBranchStrategyPseudo(this, 1);
    }
    else if (brStrategy == DcoBranchingStrategyReliability) {
        // Relibility
        branchStrategy_ = new DcoBranchStrategyRel(this, relibility);
    }
    else if (brStrategy == DcoBranchingStrategyStrong) {
        // Strong
        branchStrategy_ = new DcoBranchStrategyStrong(this);
    }
    else if (brStrategy == DcoBranchingStrategyBilevel) {
        // Bilevel
        branchStrategy_ = new DcoBranchStrategyBilevel(this);
    }
    else {
        throw CoinError("Unknown branch strategy.", "setupSelf","DcoModel");
    }

    brStrategy = DcoPar_->entry(DcoParams::branchStrategyRampUp);

    if (brStrategy == DcoBranchingStrategyMaxInfeasibility) {
        // Max inf
      rampUpBranchStrategy_ = new DcoBranchStrategyMaxInf(this);
    }
    else if (brStrategy == DcoBranchingStrategyPseudoCost) {
        // Pseudocost
        rampUpBranchStrategy_ = new DcoBranchStrategyPseudo(this, 1);
    }
    else if (brStrategy == DcoBranchingStrategyReliability) {
        // Relibility
        rampUpBranchStrategy_ = new DcoBranchStrategyRel(this, relibility);
    }
    else if (brStrategy == DcoBranchingStrategyStrong) {
        // Strong
        rampUpBranchStrategy_ = new DcoBranchStrategyStrong(this);
    }
    else if (brStrategy == DcoBranchingStrategyBilevel) {
        // Bilevel
        rampUpBranchStrategy_ = new DcoBranchStrategyBilevel(this);
    }
    else {
        throw CoinError("Unknown branch strategy.", "setupSelf","DcoModel");
    }

    //------------------------------------------------------
    // Add heuristics.
    //------------------------------------------------------

    heurStrategy_ =
       static_cast<DcoHeurStrategy> (DcoPar_->entry(DcoParams::heurStrategy));
    heurCallFrequency_ =
       static_cast<DcoHeurStrategy> (DcoPar_->entry(DcoParams::heurCallFrequency));
    DcoHeurStrategy heurRoundStrategy =
       static_cast<DcoHeurStrategy> (DcoPar_->entry(DcoParams::heurRoundStrategy));
    int callFreq = 1;

    if (heurRoundStrategy == DcoHeurStrategyNotSet) {
        heurRoundStrategy = heurStrategy_;
	callFreq = heurCallFrequency_;
    }
    if (heurRoundStrategy != DcoHeurStrategyNone &&
	heurRoundStrategy != DcoHeurStrategyNotSet) {
        // Add rounding heuristic
        DcoHeurRound *heurRound = new DcoHeurRound(this,
                                                     "Rounding",
                                                     heurRoundStrategy,
						     callFreq);
        addHeuristic(heurRound);
    }

    // Adjust heurStrategy
    for (j = 0; j < numHeuristics_; ++j) {
        if (heuristics_[j]->strategy() != DcoHeurStrategyNone) {
            // Doesn't matter what's the strategy, we just want to
            // call heuristics.
            heurStrategy_ = DcoHeurStrategyAuto;
            break;
        }
    }

    //------------------------------------------------------
    // Cut generators settings.
    //------------------------------------------------------

    // Compute dense cutoff.

    const CoinPackedMatrix * rowMatrix = lpSolver_->getMatrixByRow();
    const int * rowLen = rowMatrix->getVectorLengths();
    double maxLen = 0.0, minLen = ALPS_DBL_MAX, sumLen = 0.0;
    double aveLen, diffLen, stdLen;
    double denseConFactor = DcoPar_->entry(DcoParams::denseConFactor);

    for (j = 0; j < numRows_; ++j) {
	if (rowLen[j] > maxLen) maxLen = rowLen[j];
	if (rowLen[j] < minLen) minLen = rowLen[j];
	sumLen += rowLen[j];
    }
    assert(numRows_ > 0);
    aveLen = sumLen / numRows_;
    sumLen = 0.0;

    for (j = 0; j < numRows_; ++j) {
	diffLen = rowLen[j] - aveLen;
	diffLen *= diffLen;
	sumLen += diffLen;
    }
    stdLen = sqrt(sumLen/numRows_);
    if (denseConFactor > 10e5) {
        denseConCutoff_ = ALPS_INT_MAX;
    }
    else {
        denseConCutoff_ = static_cast<int>(aveLen + denseConFactor*stdLen);
        denseConCutoff_ = ALPS_MIN(numCols_/2, denseConCutoff_);
        denseConCutoff_ = ALPS_MAX(100, denseConCutoff_);
    }

#ifdef DISCO_DEBUG
    std::cout << "aveLen=" << aveLen << ", minLen=" << minLen
	      << ", maxLen=" << maxLen << ", stdLen=" << stdLen
	      << ", denseConCutoff_=" << denseConCutoff_ << std::endl;
#endif

    // NOTE: maxNumCons is valid only for automatic strategy.
    double cutFactor = DcoPar_->entry(DcoParams::cutFactor);

    if (cutFactor > 1.0e5) {
        maxNumCons_ = ALPS_INT_MAX;
	// old constraint size will automatically double if no space.
	oldConstraintsSize_ = 10000;
    }
    else {
        maxNumCons_ = (int)((cutFactor - 1.0) * numRows_);
	oldConstraintsSize_ = maxNumCons_;
    }

    constraintPool_ = new BcpsConstraintPool();
    constraintPoolReceive_ = new BcpsConstraintPool();
    constraintPoolSend_ = new BcpsConstraintPool();
    oldConstraints_ = new DcoConstraint* [oldConstraintsSize_];

    cutStrategy_ = static_cast<DcoCutStrategy>
	(DcoPar_->entry(DcoParams::cutStrategy));
    cutGenerationFrequency_ = static_cast<DcoCutStrategy>
	(DcoPar_->entry(DcoParams::cutGenerationFrequency));

    if (cutGenerationFrequency_ < 1) {
	std::cout << "WARNING: Input cut generation frequency is "
		  << cutGenerationFrequency_
		  << ", which is not allowed. Changed it to 1" << std::endl;
	cutGenerationFrequency_ = 1;
    }
#if 0
    std::cout << "Initially, cutStrategy_ = " << cutStrategy_
	      << "; freq = " << cutGenerationFrequency_
	      << std::endl;
#endif

    DcoCutStrategy cliqueStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutCliqueStrategy));
    DcoCutStrategy fCoverStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutFlowCoverStrategy));
    DcoCutStrategy gomoryStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutGomoryStrategy));
    DcoCutStrategy knapStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutKnapsackStrategy));
    DcoCutStrategy mirStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutMirStrategy));
    DcoCutStrategy oddHoleStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutOddHoleStrategy));
    DcoCutStrategy probeStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutProbingStrategy));
    DcoCutStrategy twoMirStrategy = static_cast<DcoCutStrategy>
       (DcoPar_->entry(DcoParams::cutTwoMirStrategy));

    int cliqueFreq = DcoPar_->entry(DcoParams::cutCliqueFreq);
    int fCoverFreq = DcoPar_->entry(DcoParams::cutFlowCoverFreq);
    int gomoryFreq = DcoPar_->entry(DcoParams::cutGomoryFreq);
    int knapFreq = DcoPar_->entry(DcoParams::cutKnapsackFreq);
    int mirFreq = DcoPar_->entry(DcoParams::cutMirFreq);
    int oddHoleFreq = DcoPar_->entry(DcoParams::cutOddHoleFreq);
    int probeFreq = DcoPar_->entry(DcoParams::cutProbingFreq);
    int twoMirFreq = DcoPar_->entry(DcoParams::cutTwoMirFreq);

    // conic cut parameters
    conicCutStrategy_ = static_cast<DcoConicCutStrategy>
	(DcoPar_->entry(DcoParams::conicCutStrategy));
    conicCutGenerationFrequency_ = DcoPar_->entry(DcoParams::cutGenerationFrequency);

    if (conicCutGenerationFrequency_ < 1) {
	std::cout << "WARNING: Input conic cut generation frequency is "
		  << cutGenerationFrequency_
		  << ", which is not allowed. Changed it to 1" << std::endl;
	conicCutGenerationFrequency_ = 1;
    }


    //------------------------------------------------------
    // Add cut generators.
    //------------------------------------------------------


    //----------------------------------
    // Add probe cut generator.
    //----------------------------------

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
	addCutGenerator(probing, "Probing", probeStrategy, probeFreq);
    }

    //----------------------------------
    // Add clique cut generator.
    //----------------------------------

    if (cliqueStrategy == DcoCutStrategyNotSet) {
        // Only at root by default
	if (cutStrategy_ == DcoCutStrategyNotSet) {
	    cliqueStrategy = DcoCutStrategyRoot;
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
	addCutGenerator(cliqueCut, "Clique", cliqueStrategy, cliqueFreq);
    }

    //----------------------------------
    // Add odd hole cut generator.
    //----------------------------------

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
	addCutGenerator(oldHoleCut, "OddHole",oddHoleStrategy,oddHoleFreq);
    }

    //----------------------------------
    // Add flow cover cut generator.
    //----------------------------------

    if (fCoverStrategy == DcoCutStrategyNotSet) {
	if (cutStrategy_ == DcoCutStrategyNotSet) {
	    fCoverStrategy = DcoCutStrategyAuto;
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
	addCutGenerator(flowGen, "Flow Cover", fCoverStrategy, fCoverFreq);
    }

    //----------------------------------
    // Add knapsack cut generator.
    //----------------------------------

    if (knapStrategy == DcoCutStrategyNotSet) {
	if (cutStrategy_ == DcoCutStrategyNotSet) {
	    // Only at root by default
	    knapStrategy = DcoCutStrategyRoot;
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
	addCutGenerator(knapCut, "Knapsack", knapStrategy, knapFreq);
    }

    //----------------------------------
    // Add MIR cut generator.
    //----------------------------------

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
	addCutGenerator(mixedGen, "MIR", mirStrategy, mirFreq);
    }

    //----------------------------------
    // Add Gomory cut generator.
    //----------------------------------

    if (gomoryStrategy == DcoCutStrategyNotSet) {
	if (cutStrategy_ == DcoCutStrategyNotSet) {
	    // Only at root by default
	    gomoryStrategy = DcoCutStrategyRoot;
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
	addCutGenerator(gomoryCut, "Gomory", gomoryStrategy, gomoryFreq);
    }

    //----------------------------------
    // Add Tow MIR cut generator.
    //----------------------------------

    // Disable forever, not useful.
    twoMirStrategy = DcoCutStrategyNone;
    if (twoMirStrategy != DcoCutStrategyNone) {
        CglTwomir *twoMirCut =  new CglTwomir;
	addCutGenerator(twoMirCut, "Two MIR", twoMirStrategy, twoMirFreq);
    }

#if defined(__COLA__)
    //----------------------------------
    // Add Outer approximation cut generator
    //----------------------------------
    CglConicCutGenerator * oa_gen = new CglConicOA();
    DcoConGeneratorBase * dco_cg = new DcoConicConGenerator(this,
					   oa_gen, "oa_gen");
    addCutGenerator(dco_cg);
    oa_gen = NULL;
#endif

    //--------------------------------------------
    // Adjust cutStrategy_ according to the strategies of
    // each cut generators.
    //--------------------------------------------

    if (numCutGenerators_ > 0) {
	DcoCutStrategy strategy0 = cutGenerators(0)->strategy();
	DcoCutStrategy strategy1;
	for (j = 1; j < numCutGenerators_; ++j) {
	    strategy1 = cutGenerators(j)->strategy();
	    if (strategy1 != strategy0) {
	       // A generator has different strategy.
		break;
	    }
	}
	if (j == numCutGenerators_) {
	    // All cut generators has same strategy.
	    cutStrategy_ = strategy0;
	}
        else {
            // Assume to generate cons at each node since generators
	    // has various strategies.
            cutStrategy_ = DcoCutStrategyPeriodic;
	    cutGenerationFrequency_ = 1;
        }
    }

    // adjust conic cut strategy
    if (conicCutStrategy_==DcoConicCutStrategyNotSet) {
      conicCutStrategy_ = DcoConicCutStrategyNone;
    }

#if 0
    std::cout << "AFTER: cutStrategy_ = " << cutStrategy_ << std::endl;
#endif

    //--------------------------------------------
    // Random vector
    //--------------------------------------------

    // TODO, if generating variable, then need more space.
    conRandoms_ = new double [numCols_];
    double deseed = 12345678.0;
    double DSEED2 = 2147483647.0;

    for (j = 0; j < numCols_; ++j) {
        deseed *= 16807.0;
        int jseed = (int) (deseed / DSEED2);
        deseed -= (double) jseed * DSEED2;
        double random = deseed / DSEED2;
        conRandoms_[j] = random;
        //std::cout << "conRandoms_[" << j << "]="
        //      <<conRandoms_[j]<< std::endl;
    }


    return true;
}

//############################################################################

// AT - Begin
void DcoModel::presolveForTheWholeTree() {
 	 //if (!lpSolver_) {
        // preprocessing causes this check.
     //   lpSolver_ = origLpSolver_;
    //}
    int numPasses = 50;
    double feaTol = 1.0e-3;
    bool keepIntegers = true;
    char *prohibited = 0;
    bool doPresolve = DcoPar_->entry(DcoParams::presolve);
#ifdef DISCO_DEBUG
	std::cout << "Presolve = "<< doPresolve<<" problem setup " <<problemSetup << std::endl;
#endif

	if (doPresolve) {
		// Initialize solver with original data
		//AT
		std::cout<<" About to initialize problem with original data"<<std::endl;
		origLpSolver_->loadProblem(*colMatrix_,varLB_, varUB_,objCoef_,  conLB_,conUB_);
		//AT
		std::cout<<" Problem initialized "<<std::endl;
		std::cout<<" Preprocessing "<<std::endl;
		presolvedLpSolver_ = presolve_->preprocess( *origLpSolver_,
			feaTol,
			keepIntegers,
			numPasses,
			prohibited);
//AT DEBUG
		#ifdef DISCO_DEBUG
		for(int i =0;i<presolvedLpSolver_->getNumCols();i++)
			if(presolvedLpSolver_->isInteger(i)) {
				std::cout<<" Variable "<<i<<" Is still Integer"<<std::endl;
			}
		#endif

		std::cout<<" Preprocessing  done"<<std::endl;
		presolvedLpSolver_->initialSolve();
		std::cout<<" Solved "<<std::endl;

		colMatrix_=presolvedLpSolver_->getMutableMatrixByCol();
		numCols_=presolvedLpSolver_->getNumCols();
		numRows_=presolvedLpSolver_->getNumRows();
		  memcpy(varLB_, presolvedLpSolver_->getColLower(), sizeof(double) * numCols_);
    memcpy(varUB_, presolvedLpSolver_->getColUpper(), sizeof(double) * numCols_);
    memcpy(conLB_, presolvedLpSolver_->getRowLower(), sizeof(double) * numRows_);
    memcpy(conUB_, presolvedLpSolver_->getRowUpper(), sizeof(double) * numRows_);

		 if (objSense_ > 0.0) {
        memcpy(objCoef_, presolvedLpSolver_->getObjCoefficients(), sizeof(double) * numCols_);
    }
    else {
        const double *mpsObj =  presolvedLpSolver_->getObjCoefficients();
        for (int j = 0; j < numCols_; ++j) {
            objCoef_[j] = - mpsObj[j];

		}
    }


		lpSolver_ = presolvedLpSolver_->clone(true);//->clone();
		setSolver(lpSolver_);

		presolved=true;
   }

}
// AT - End

//############################################################################
// AT - Begin -- not used
void DcoModel::preprocess()
{
    int numPasses = 50;
    double feaTol = 1.0e-3;
    bool keepIntegers = true;
    char *prohibited = 0;
	// Do presolve only if problem already setup - In Alp preprocess is called before...
    bool doPresolve = DcoPar_->entry(DcoParams::presolve) && problemSetup;
#ifdef DISCO_DEBUG
	std::cout << "Presolve = "<< doPresolve<<" problem setup " <<problemSetup << std::endl;
#endif
//AT disabled for now
	if(true) return;
	if (doPresolve) {
		presolvedLpSolver_ = presolve_->preprocess( *origLpSolver_,
						   feaTol,
						   keepIntegers,
						   numPasses,
						   prohibited);
		presolvedLpSolver_->initialSolve();
		lpSolver_ = presolvedLpSolver_->clone();
		#ifdef DISCO_DEBUG
			std::cout << "Presolve done  "<< presolvedLpSolver_->isProvenOptimal() << std::endl;
		#endif
		presolved=true;
   }
    else {
	lpSolver_ = origLpSolver_;
    }
}
// AT - End

//############################################################################


// AT - Begin - not used before, do postsolve
void DcoModel::postprocess()
{

//std::cout<<" Lp before post  col "<<lpSolver_->getNumCols()<<" Rows "<<lpSolver_->getNumRows()<<std::endl;
  if(!DcoPar_->entry(DcoParams::presolve))
    return;
  std::cout<<" POST SOLVING "<<std::endl;
  std::cout<<" Original Model  col "<<origLpSolver_->getNumCols()<<" Rows "<<origLpSolver_->getNumRows()<<std::endl;

  numCols_=origLpSolver_->getNumCols();
  DcoSolution * sol =dynamic_cast<DcoSolution*>(getKnowledgeBroker()->getBestKnowledge(AlpsKnowledgeTypeSolution).first);
  presolve_->model()->setColSolution(sol->getValues());
  presolve_->postprocess();
  std::cout<<" Sol size "<<sol->getSize()<< " Quality "<<sol->getQuality() <<std::endl;
  std::cout<< " Original model val: "<<presolve_->originalModel()->getObjValue()<<std::endl;
  const double * values = presolve_->originalModel()->getColSolution();

  DcoSolution *  currentSol = new DcoSolution(presolve_->originalModel()->getNumCols(),values,presolve_->originalModel()->getObjValue());
  sol->setQuality(currentSol->getQuality());
  sol->setSize(currentSol->getSize());
  sol->setValues(currentSol->getValues(),currentSol->getSize());
  /*std::cout<<"Num Knowledge " << getKnowledgeBroker()->getKnowledgePool(AlpsKnowledgeTypeSolution)->getNumKnowledges()<<std::endl;

  getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
		currentSol,
		objSense_ * currentSol->getQuality());
  getKnowledgeBroker()->getKnowledgePool(AlpsKnowledgeTypeSolution)->cl;
  std::cout<<"Num Knowledge " << getKnowledgeBroker()->getKnowledgePool(AlpsKnowledgeTypeSolution)->getNumKnowledges()<<std::endl;
	*/

 // lpSolver_=presolve_->model();//->clone();
//  std::cout<<" Lp after post  col "<<lpSolver_->getNumCols()<<" Rows "<<lpSolver_->getNumRows()<<std::endl;
 // std::cout<<" Presolved after post  col "<<presolvedLpSolver_ ->getNumCols()<<" Rows "<<presolvedLpSolver_ ->getNumRows()<<std::endl;
  //setSolver(lpSolver_);

 // lpSolver_=origLpSolver_->clone();
#ifdef DISCO_DEBUG
  std::cout<<"POST process\n Num col "<<numCols_<<" Presolved "<<lpSolver_->getNumCols()
	   <<"\nNum row "<<numRows_<<" Presolved "<<lpSolver_->getNumRows()
	   <<std::endl;
#endif
}
// AT - End

//#############################################################################

int
DcoModel::storeSolution(DcoSolutionType how, DcoSolution* sol)
{
    double quality = sol->getQuality();

    double cutoff = getCutoff();

    // Update cutoff and lp cutoff.
    setCutoff(quality);

    ++numSolutions_;

    // Store in pool, assume minimization.
    getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
                                       sol,
                                       objSense_ * sol->getQuality());
#if 0
        std::cout << "Found a better solution"
                  <<", old cutoff = " << cutoff
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif

    // Record how the solution was found
    switch (how) {
    case DcoSolutionTypeBounding:
#ifdef DISCO_DEBUG
        std::cout << "Rounding heuristics found a better solution"
                  <<", old cutoff = " << cutoff
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    case DcoSolutionTypeBranching:
#ifdef DISCO_DEBUG
        std::cout << "Branching found a better solution"
                  <<", old cutoff = " << cutoff
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    case DcoSolutionTypeDiving:
        ++numHeurSolutions_;
#ifdef DISCO_DEBUG
        std::cout << "Branching found a better solution"
                  <<", old cutoff = " << cutoff
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    case DcoSolutionTypeHeuristic:
        ++numHeurSolutions_;
        if (broker_->getMsgLevel() > 200) {
            std::cout << "Heuristics found a better solution"
                      <<", old cutoff = " << cutoff
                      << ", new cutoff = " << getCutoff()  << std::endl;
        }
        break;
    case DcoSolutionTypeStrong:
#ifdef DISCO_DEBUG
        std::cout << "Strong branching found a better solution"
                  <<", old cutoff = " << cutoff
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    default:
#ifdef DISCO_DEBUG
        std::cout << "Nowhere found a better solution"
                  <<", old cutoff = " << cutoff
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    }

    return AlpsReturnStatusOk;
}

//############################################################################
// todo(aykut) rename it to createIntegerObjects not intger but integer.
void
DcoModel::createIntgerObjects(bool startAgain)
{
    assert(lpSolver_);

    if (numIntObjects_ && !startAgain && objects_) return;

    int iCol;
    int numCols = getNumCols();

    const double *colLB = lpSolver_->getColLower();
    const double *colUB = lpSolver_->getColUpper();
    DcoObjectInt *intObject = NULL;

    if (intColIndices_) {
        delete [] intColIndices_;
    }
    numIntObjects_ = 0;

    for (iCol = 0; iCol < numCols; ++iCol) {
	if (lpSolver_->isInteger(iCol)) ++numIntObjects_;
    }

#if 0
    std::cout << "==== findInteger: numIntObjects_ = "
              << numIntObjects_ << std::endl;
#endif

    double weight = DcoPar_->entry(DcoParams::pseudoWeight);

    int numObjects = 0;
    int iObject;
    BcpsObject ** oldObject = objects_;

    for (iObject = 0; iObject < numObjects_; ++iObject) {
	DcoObjectInt * obj =
	    dynamic_cast <DcoObjectInt *>(oldObject[iObject]) ;

	if (obj) {
	    delete oldObject[iObject];
        }
	else {
	    oldObject[numObjects++] = oldObject[iObject];
        }
    }

    if (!intObjIndices_) {
        intObjIndices_ = new int [numCols];
        memset(intObjIndices_, 0, sizeof(int) * numCols);
    }

    objects_ = new BcpsObject * [(numIntObjects_ + numObjects)];
    intColIndices_ = new int [numIntObjects_];
    numObjects_ = numIntObjects_ + numObjects;

    // Walk the variables again, filling in the indices and creating objects
    // for the integer variables. Initially, the objects hold the indices,
    // variable bounds and pseudocost.
    numIntObjects_ = 0;
    for (iCol = 0; iCol < numCols; ++iCol) {
	if(lpSolver_->isInteger(iCol)) {

	    intObject = new DcoObjectInt(numIntObjects_,
                                          iCol,
                                          colLB[iCol],
                                          colUB[iCol]);
            intObject->pseudocost().setWeight(weight);

	    intObjIndices_[iCol] = numIntObjects_;
            objects_[numIntObjects_] = intObject;
	    intColIndices_[numIntObjects_++] = iCol;
	}
    }

    if (numIntObjects_) {
        sharedObjectMark_ = new char [numIntObjects_];
        memset(sharedObjectMark_, 0, sizeof(char) * numIntObjects_);
    }

    // Now append other objects
    memcpy(objects_ + numIntObjects_, oldObject, numObjects*sizeof(BcpsObject *));

    // Delete old array (just array)
    delete [] oldObject;
}

//#############################################################################

bool
DcoModel::resolve()
{
    lpSolver_->resolve();
    numIterations_ += lpSolver_->getIterationCount();
    bool feasible = (lpSolver_->isProvenOptimal() &&
                     !lpSolver_->isDualObjectiveLimitReached());

    return feasible;
}

//#############################################################################

// Delete all object information
void
DcoModel::deleteObjects()
{
    delete [] priority_;
    priority_ = NULL;
    int i;
    for (i = 0; i < numObjects_; ++i) delete objects_[i];
    delete [] objects_;
    objects_ = NULL;
    numObjects_ = 0;
    createIntgerObjects(true);
}

//#############################################################################

DcoModel::~DcoModel()
{
    gutsOfDestructor();
}

//#############################################################################

void
DcoModel::gutsOfDestructor()
{
    int i;
    bool doPresolve = DcoPar_->entry(DcoParams::presolve);


//    delete [] savedLpSolution_;
//    savedLpSolution_ = NULL;

    delete [] intObjIndices_;
    intObjIndices_ = NULL;

    delete [] intColIndices_;
    intColIndices_ = NULL;

    for (i = 0; i < numObjects_; ++i) delete objects_[i];
    delete [] objects_;
    objects_ = NULL;

    delete [] priority_;
    priority_ = NULL;

    delete [] colType_;
    delete colMatrix_;

    delete [] varLB_;
    delete [] varUB_;

    delete [] conLB_;
    delete [] conUB_;

    delete [] startVarLB_;
    delete [] startVarUB_;

    delete [] startConLB_;
    delete [] startConUB_;

    delete [] tempVarLBPos_;
    delete [] tempVarUBPos_;

    delete [] tempConLBPos_;
    delete [] tempConUBPos_;

    delete [] objCoef_;
    delete [] incumbent_;

    // free cone data
    for (int i=0; i<numCoreCones_; ++i) {
      delete[] coneMembers_[i];
    }
    if (coneMembers_) {
      delete[] coneMembers_;
    }
    if (coneTypes_) {
      delete[] coneTypes_;
    }
    if (coneSizes_) {
      delete[] coneSizes_;
    }

    delete presolve_;

    if (numHeuristics_ > 0) {
#ifdef DISCO_DEBUG
        std::cout << "MODEL: distructor: numHeuristics =  "
                  << numHeuristics_ << std::endl;
#endif
        DcoHeuristic *tempH = NULL;

        for (i = 0; i < numHeuristics_; ++i) {
            tempH = heuristics_[i];
            delete tempH;
        }

        delete [] heuristics_;
        heuristics_ = NULL;
    }

    if (generators_ != NULL) {
#ifdef DISCO_DEBUG
        std::cout << "MODEL: distructor: numCutGenerators = "
                  << numCutGenerators_ << std::endl;
#endif
        DcoConGeneratorBase *temp = NULL;
        for (i = 0; i < numCutGenerators_; ++i) {
            temp = generators_[i];
            delete temp;
        }
        delete [] generators_;
        generators_ = NULL;
    }

    delete constraintPool_;
    delete constraintPoolReceive_;
    delete constraintPoolSend_;
    delete [] oldConstraints_;
    delete branchStrategy_;
    delete rampUpBranchStrategy_;

    delete [] conRandoms_;

    delete DcoPar_;
    delete blisMessageHandler_;
    // AT  - delete structure only if problem has been really presolved, for parallel code
    //	std::cout<<" DEleteing "<<std::endl;
    if (doPresolve && presolved) {
	//delete presolvedLpSolver_;
	delete lpSolver_;
    }

    delete [] sharedObjectMark_;
    sharedObjectMark_ = NULL;
}

//#############################################################################

DcoSolution *
DcoModel::feasibleSolutionHeur(const double *solution)
{
    int j, ind;

    bool feasible = true;
    bool userFeasible = true;

    double tol = 1.0e-6;
    double value, nearest, objValue = 0.0;
    double *rowAct = NULL;

    DcoSolution *blisSol = NULL;

    // Check if within column bounds
    for (j = 0; j < numCols_; ++j) {
        value = solution[j];
        if (varLB_[j] > -ALPS_INFINITY) {
            if (value < varLB_[j] - tol) {
                feasible = false;
                goto TERM_FEAS_HEUR;
            }
        }
        if (varUB_[j] < ALPS_INFINITY) {
            if (value > varUB_[j] + tol) {
                feasible = false;
                goto TERM_FEAS_HEUR;
            }
        }
        objValue += value * objCoef_[j];
    }

    if (broker_->getMsgLevel() > 200) {
        std::cout << "FEASIBLE HEUR: numCols = " << numCols_
                  << " ; pass column bound check." << std::endl;
    }

    // Check if integeral
    for (j = 0; j < numIntObjects_; ++j) {
        ind = intColIndices_[j];
        value = solution[ind];
        //std::cout << "ind = " << ind << " ; value = " << value << std::endl;

        nearest = static_cast<int>(value + 0.5);
        if (fabs(value - nearest) > integerTol_) {
            feasible = false;
            goto TERM_FEAS_HEUR;
        }
    }

    if (broker_->getMsgLevel() > 200) {
        std::cout << "FEASIBLE HEUR: numInts = " << numIntObjects_
                  << " ; pass integral check." << std::endl;
    }

    // Check if within row bounds
    rowAct = new double [numRows_];
    colMatrix_->times(solution, rowAct);
    for (j = 0; j < numRows_; ++j) {
        value = rowAct[j];
        if (conLB_[j] > -ALPS_INFINITY) {
            if (value < conLB_[j] - tol) {
                feasible = false;
                goto TERM_FEAS_HEUR;
            }
        }
        if (conUB_[j] < ALPS_INFINITY) {
            if (value > conUB_[j] + tol) {
                feasible = false;
                goto TERM_FEAS_HEUR;
            }
        }
    }

    if (broker_->getMsgLevel() > 200) {
        std::cout << "FEASIBLE HEUR: numRows = " << numRows_
                  << " ; pass row bounds check." << std::endl;
    }

TERM_FEAS_HEUR:

    // Check if satisfy user criteria
    if (feasible) {
        blisSol = userFeasibleSolution(solution, userFeasible);

        if (broker_->getMsgLevel() > 200 && blisSol && userFeasible) {
            std::cout << "FEASIBLE HEUR: pass user feasibility check." << std::endl;
        }

	if (!blisSol && userFeasible) {
	    // User doesn't provide feasibility check.
	    numBranchResolve_ = 10;
	}
    }

    if (feasible && userFeasible && !blisSol) {
        // User doesn't overload feasible solution function.
        blisSol = new DcoSolution(getNumCols(), solution, objValue);
    }

    if (rowAct) delete [] rowAct;
    return blisSol;
}

//#############################################################################

DcoSolution *
DcoModel::feasibleSolution(int & numIntegerInfs, int & numObjectInfs)
{
    int preferredWay, j;
    int numUnsatisfied = 0;

    bool userFeasible = true;
    double sumUnsatisfied = 0.0;

    DcoSolution* sol = NULL;

    for (j = 0; j < numIntObjects_; ++j) {
	const BcpsObject * object = objects_[j];

	double infeasibility = object->infeasibility(this, preferredWay);
	if (infeasibility) {
	    assert (infeasibility>0);
	    numUnsatisfied++;
	    sumUnsatisfied += infeasibility;
	}
    }
    numIntegerInfs = numUnsatisfied;
    for (; j < numObjects_; ++j) {
	const BcpsObject * object = objects_[j];
	double infeasibility = object->infeasibility(this, preferredWay);
	if (infeasibility) {
	    assert (infeasibility > 0);
	    numUnsatisfied++;
	    sumUnsatisfied += infeasibility;
	}
    }

    numObjectInfs = numUnsatisfied - numIntegerInfs;

    // check feasibility of conic constraints.
    int start = numRows_;
    int end = numRows_ + numCoreCones_;
    for (j=start; j<end; j++) {
      double infeas;
      DcoConicConstraint const * con;
      con = dynamic_cast<DcoConicConstraint*>(constraints_[j]);
      infeas = con->infeasibility(this, preferredWay);
      if (infeas) {
	sumUnsatisfied += infeas;
	numUnsatisfied++;
      }
    }
    if (broker_->getMsgLevel() > 200) {
	std::cout << "FEASIBLE SOL: numUnsatisfied = "
		  << numUnsatisfied << std::endl;
    }

    if (!numUnsatisfied) {
        sol = userFeasibleSolution(getLpSolution(), userFeasible);
        if (broker_->getMsgLevel() > 200 && sol && userFeasible) {
            std::cout << "FEASIBLE HEUR: pass user feasibility check." << std::endl;
        }
	if (!sol && userFeasible) {
	    // User doesn't provide feasibility check.
	    numBranchResolve_ = 10;
	}
    }

    if (!numUnsatisfied && userFeasible && !sol) {
        // User doesn't overload feasible solution function.
        sol = new DcoSolution(getNumCols(),
                               getLpSolution(),
                               objSense_ * getLpObjValue());
    }

    return sol;
}

//#############################################################################
/*
  Set branching priorities.

  Setting integer priorities looks pretty robust; the call to createIntgerObjects
  makes sure that integer objects are in place. Setting priorities for
  other objects is entirely dependent on their existence, and the routine may
  quietly fail in several directions.
*/
void
DcoModel::passInPriorities (const int * priorities,
			     bool ifObject,
			     int defaultValue)
{

    // FIXME: not completed.
    int i;

    createIntgerObjects(false);

    if (!priority_) {
	priority_ = new int[numObjects_];
	for (i = 0; i < numObjects_; ++i) {
            priority_[i] = defaultValue;
        }
    }

    if (priorities) {
	if (ifObject) {
	    memcpy(priority_ + numIntObjects_,
                   priorities,
		   (numObjects_ - numIntObjects_) * sizeof(int));
        }
	else {
	    memcpy(priority_, priorities, numIntObjects_ * sizeof(int));
        }
    }
}

//#############################################################################

// todo(aykut) we need conic constraints here.
// a node should track of conic constraints
AlpsTreeNode *
DcoModel::createRoot() {

    //-------------------------------------------------------------
    // NOTE: Root will be deleted by ALPS. Root is an explicit node.
    //-------------------------------------------------------------

    DcoTreeNode* root = new DcoTreeNode;
    DcoNodeDesc* desc = new DcoNodeDesc(this);
    root->setDesc(desc);

    //-------------------------------------------------------------
    // NOTE: Although original data are stored in model when reading.
    //   Root desc still store a full copy of col and row bounds when creating.
    //   It will store soft differences after finding a branching object.
    //   The soft difference are due to reduced cost fixing and probing.
    //   Also the added cols and rows will be stored.
    //-------------------------------------------------------------
    int k;

    std::vector<BcpsVariable *> vars = getVariables();
    std::vector<BcpsConstraint *> cons = getConstraints();

    int numVars = static_cast<int> (vars.size());
    int numCons = static_cast<int> (cons.size());

#ifdef DISCO_DEBUG
    std::cout << "DISCO: createRoot(): numVars=" << numVars
              << ", numCons=" << numCons
              << "; numCoreVariables_=" << numCoreVariables_
              << ", numCoreConstraints_=" << numCoreConstraints_ << std::endl;
#endif

    int *varIndices1 = new int [numVars];
    int *varIndices2 = new int [numVars];
    int *varIndices3 = NULL;
    int *varIndices4 = NULL;
    double *vlhe = new double [numVars];
    double *vuhe = new double [numVars];
    double *vlse = NULL;
    double *vuse = NULL;

    int *conIndices1 = new int [numCons];
    int *conIndices2 = new int [numCons];
    int *conIndices3 = NULL;
    int *conIndices4 = NULL;
    double *clhe = new double [numCons];
    double *cuhe = new double [numCons];
    double *clse = NULL;
    double *cuse = NULL;

    //-------------------------------------------------------------
    // Get var bounds and indices.
    //-------------------------------------------------------------

    for (k = 0; k < numVars; ++k) {
        vlhe[k] = vars[k]->getLbHard();
        vuhe[k] = vars[k]->getUbHard();
        varIndices1[k] = k;
        varIndices2[k] = k;

#ifdef DISCO_DEBUG_MORE
        std::cout << "DISCO: createRoot(): var "<< k << ": hard: lb=" << vlhe[k]
                  << ", ub=" << vuhe[k] << std::endl;
#endif
    }

    //-------------------------------------------------------------
    // Get con bounds and indices.
    //-------------------------------------------------------------

    for (k = 0; k < numCons; ++k) {
        clhe[k] = cons[k]->getLbHard();
        cuhe[k] = cons[k]->getUbHard();
        conIndices1[k] = k;
        conIndices2[k] = k;
    }

    int *tempInd = NULL;
    BcpsObject **tempObj = NULL;

    desc->assignVars(0 /*numRem*/, tempInd,
                     0 /*numAdd*/, tempObj,
                     false, numVars, varIndices1, vlhe, /*Var hard lb*/
                     false, numVars, varIndices2, vuhe, /*Var hard ub*/
                     false, 0, varIndices3, vlse,       /*Var soft lb*/
                     false, 0, varIndices4, vuse);      /*Var soft ub*/
    desc->assignCons(0 /*numRem*/, tempInd,
                     0 /*numAdd*/, tempObj,
                     false, numCons, conIndices1, clhe, /*Con hard lb*/
                     false, numCons, conIndices2, cuhe, /*Con hard ub*/
                     false, 0,conIndices3,clse,         /*Con soft lb*/
                     false, 0,conIndices4,cuse);        /*Con soft ub*/

    //-------------------------------------------------------------
    // Mark it as an explicit node.
    //-------------------------------------------------------------

    root->setExplicit(1);

    return root;
}

//#############################################################################

void
DcoModel::addHeuristic(DcoHeuristic * heuristic)
{
    DcoHeuristic ** temp = heuristics_;
    heuristics_ = new DcoHeuristic * [numHeuristics_ + 1];

    memcpy(heuristics_, temp, numHeuristics_ * sizeof(DcoHeuristic *));
    delete [] temp;

    heuristics_[numHeuristics_++] = heuristic;
}

//#############################################################################

void
DcoModel::addCutGenerator(CglCutGenerator * generator,
			   const char * name,
			   DcoCutStrategy strategy,
			   int freq,
			   bool normal,
			   bool atSolution,
			   bool whenInfeasible)
{
    DcoConGeneratorBase ** temp = generators_;

    generators_ = new DcoConGeneratorBase * [(numCutGenerators_ + 1)];
    memcpy(generators_, temp, numCutGenerators_ * sizeof(DcoConGeneratorBase *));

    generators_[numCutGenerators_++] =
        new DcoConGenerator(this, generator, name, strategy, freq,
                             normal, atSolution, whenInfeasible);
    delete [] temp;
    temp = NULL;
#if 0
    std::cout << "Added " << name << " cut generator, strategy is "
	      << strategy << ", freq is " << freq << std::endl;
#endif
}

//#############################################################################

void
DcoModel::addCutGenerator(DcoConGeneratorBase * generator)
{
    DcoConGeneratorBase ** temp = generators_;

    generators_ = new DcoConGeneratorBase * [(numCutGenerators_ + 1)];
    memcpy(generators_, temp, numCutGenerators_ * sizeof(DcoConGeneratorBase *));

    generators_[numCutGenerators_++] = generator;

    delete [] temp;
    temp = NULL;
}

//#############################################################################

AlpsReturnStatus
DcoModel::encodeDco(AlpsEncoded *encoded) const
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    DcoPar_->pack(*encoded);

    encoded->writeRep(objSense_);

    return status;
}

//#############################################################################

AlpsReturnStatus
DcoModel::decodeDco(AlpsEncoded &encoded)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    DcoPar_->unpack(encoded);

    encoded.readRep(objSense_);

    //------------------------------------------------------
    // Load problem to LP solver.
    //------------------------------------------------------

    std::vector<DcoVariable *> vars;
    std::vector<DcoConstraint *> cons;

    int k;
    int size = static_cast<int> (variables_.size());
    for (k = 0; k < size; ++k) {
        DcoVariable * aVar = dynamic_cast<DcoVariable *>(variables_[k]);
        vars.push_back(aVar);
    }

    size = static_cast<int> (constraints_.size());
    for (k = 0; k < size; ++k) {
        DcoConstraint * aCon = dynamic_cast<DcoConstraint *>(constraints_[k]);
        cons.push_back(aCon);
    }

    variables_.clear();
    constraints_.clear();

    // Import model from vars and cons.
    importModel(vars, cons);

    return status;
}

//#############################################################################

AlpsEncoded*
DcoModel::encode() const
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    // NOTE: "AlpsKnowledgeTypeModel" is the type name.
    AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeModel);

    status = encodeAlps(encoded);
    status = encodeBcps(encoded);
    status = encodeDco(encoded);

    return encoded;
}

//#############################################################################

void
DcoModel::decodeToSelf(AlpsEncoded& encoded)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    status = decodeAlps(encoded);
    status = decodeBcps(encoded);
    status = decodeDco(encoded);
}

//#############################################################################

AlpsEncoded*
DcoModel::packSharedKnowlege()
{
    AlpsEncoded* encoded = NULL;

    bool share = false;
    bool sharePseudo = false;
    bool shareCon = false;
    bool shareVar = false;

    int k = 0;
    int numPseudcosts = 0;
    int numCons = 0;
    //int numVars = 0;
    int frequency = -1, depth = -1;
    int phase = broker_->getPhase();

    //------------------------------------------------------
    // Decide if share pseudocosts.
    //------------------------------------------------------

    if (phase == AlpsPhaseRampup) {
        sharePseudo = DcoPar_->entry(DcoParams::sharePseudocostRampUp);
    }
    else if (phase == AlpsPhaseSearch) {
        sharePseudo = DcoPar_->entry(DcoParams::sharePseudocostSearch);
        if (sharePseudo) {
            // Depth and frequency
            depth = DcoPar_->entry(DcoParams::sharePcostDepth);
            frequency =  DcoPar_->entry(DcoParams::sharePcostFrequency);
            if ( /*(numNodes_ % frequency != 0) ||*/
		(broker_->getTreeDepth() >  depth) ) {
                sharePseudo = false;
            }
        }
    }

    if (sharePseudo) {
        for (k = 0; k < numIntObjects_; ++k) {
            if (sharedObjectMark_[k]){
                ++numPseudcosts;
            }
        }
        if (numPseudcosts) share = true;
    }

    //------------------------------------------------------
    // TODO: constraints, and variables, etc.
    //------------------------------------------------------

    shareCon = DcoPar_->entry(DcoParams::shareConstraints);
    numCons = constraintPoolSend_->getNumConstraints();
    if (shareCon && numCons >= DISCO_MIN_SHARE_CON) {
	share = true;
    }

    shareVar = DcoPar_->entry(DcoParams::shareVariables);
    //numVars = constraintPoolSend_->getNumVariables();
    //if (shareVar && numVars > 10) share = true;

#if 0
    std::cout << "++++ sharePseudo =  " << sharePseudo
              << ", shareCon = " << shareCon
              << ", shareVar = " << shareVar
              << ", numPseudcosts = " << numPseudcosts
              << ", numNodes_ = " << numNodes_
              << ", frequency = " << frequency
              << ", depth = " << depth
              << ", treeDepth = " << broker_->getTreeDepth()
              << ", share = " << share << std::endl;
#endif

    if (share) {
	encoded = new AlpsEncoded(AlpsKnowledgeTypeModelGen);
	packSharedPseudocost(encoded, numPseudcosts);
	packSharedConstraints(encoded);
	packSharedVariables(encoded);
    }

    return encoded;
}

//#############################################################################

void
DcoModel::unpackSharedKnowledge(AlpsEncoded& encoded)
{
    unpackSharedPseudocost(encoded);
    unpackSharedConstraints(encoded);
    unpackSharedVariables(encoded);
}

//#############################################################################

/** Register knowledge. */
void
DcoModel::registerKnowledge() {
    // Register model, solution, and tree node
    assert(broker_);
    broker_->registerClass(AlpsKnowledgeTypeModel, new DcoModel);
    if (broker_->getMsgLevel() > 100) {
	std::cout << "DISCO: Register Alps model." << std::endl;
    }

    broker_->registerClass(AlpsKnowledgeTypeNode, new DcoTreeNode(this));
    if (broker_->getMsgLevel() > 100) {
	std::cout << "DISCO: Register Alps node." << std::endl;
    }

    broker_->registerClass(AlpsKnowledgeTypeSolution, new DcoSolution);
    if (broker_->getMsgLevel() > 100) {
	std::cout << "DISCO: Register Alps solution." << std::endl;
    }

    broker_->registerClass(BcpsKnowledgeTypeConstraint, new DcoConstraint);
    if (broker_->getMsgLevel() > 100) {
	std::cout << "DISCO: Register Bcps constraint." << std::endl;
    }

    broker_->registerClass(BcpsKnowledgeTypeVariable, new DcoVariable);
    if (broker_->getMsgLevel() > 100) {
	std::cout << "DISCO: Register Bcps variable." << std::endl;
    }
}

//#############################################################################

/** Log of specific models. */
void
DcoModel::modelLog()
{
    int logFileLevel = AlpsPar_->entry(AlpsParams::logFileLevel);
    int msgLevel = AlpsPar_->entry(AlpsParams::msgLevel);

    if (broker_->getProcType() == AlpsProcessTypeSerial) {
        /* Only valid for serial code */
        if (logFileLevel > 0) {
            std::string logfile = AlpsPar_->entry(AlpsParams::logFile);
            std::ofstream logFout(logfile.c_str(), std::ofstream::app);
            writeParameters(logFout);
        }

        if (msgLevel > 0) {
            int k;
            for (k = 0; k < numCutGenerators_; ++k) {
                if (cutGenerators(k)->calls() > 0) {
                    blisMessageHandler()->message(DISCO_CUT_STAT_FINAL,
                                                  blisMessages())
                        << cutGenerators(k)->name()
                        << cutGenerators(k)->calls()
                        << cutGenerators(k)->numConsGenerated()
                        << cutGenerators(k)->time()
                        << cutGenerators(k)->strategy()
                        << CoinMessageEol;
                }
            }
            for (k = 0; k < numHeuristics_; ++k) {
                if (heuristics(k)->calls() > 0) {
                    blisMessageHandler()->message(DISCO_HEUR_STAT_FINAL,
                                                  blisMessages())
                        << heuristics(k)->name()
                        << heuristics(k)->calls()
                        << heuristics(k)->numSolutions()
                        << heuristics(k)->time()
                        << heuristics(k)->strategy()
                        << CoinMessageEol;
                }
            }

            // Print gap
            if (currRelGap_ > ALPS_OBJ_MAX_LESS) {
                blisMessageHandler()->message(DISCO_GAP_NO, blisMessages())
                    << CoinMessageEol;
            }
            else {
                blisMessageHandler()->message(DISCO_GAP_YES, blisMessages())
                    << currRelGap_ << CoinMessageEol;
            }

        }
    }
    else if (broker_->getProcType() == AlpsProcessTypeMaster) {
        if (msgLevel > 0) {
            // Print gap
            if (currRelGap_ > ALPS_OBJ_MAX_LESS) {
                blisMessageHandler()->message(DISCO_GAP_NO, blisMessages())
                    << CoinMessageEol;
            }
            else {
                blisMessageHandler()->message(DISCO_GAP_YES, blisMessages())
                    << currRelGap_ << CoinMessageEol;
            }
        }
    }
}

//#############################################################################

void
DcoModel::nodeLog(AlpsTreeNode *node, bool force)
{
    int nodeInterval =
	broker_->getModel()->AlpsPar()->entry(AlpsParams::nodeLogInterval);

    int numNodesProcessed = broker_->getNumNodesProcessed();
    int numNodesLeft = broker_->updateNumNodesLeft();
    int msgLevel = broker_->getMsgLevel();

    bool printLog = false;
    double feasBound = ALPS_OBJ_MAX;
    double relBound = ALPS_OBJ_MAX;
    double gapVal = ALPS_OBJ_MAX;

    AlpsTreeNode *bestNode = NULL;

    if (broker_->getProcType() == AlpsProcessTypeSerial) {
        /* For serial code only */
#if 0
        std::cout << "==== Process " << broker_->getProcRank()
                  << ": printLog = " << printLog
                  << ", msgLevel = " << msgLevel
                  << ", proc type = " << broker_->getProcType()
                  << std::endl;
#endif

        // Get gap
        if (broker_->getNumKnowledges(AlpsKnowledgeTypeSolution) > 0) {
            feasBound = (broker_->getBestKnowledge(AlpsKnowledgeTypeSolution)).second;
        }

        bestNode = broker_->getBestNode();

        if (bestNode) {
            relBound = bestNode->getQuality();
        }

        if (relBound > ALPS_OBJ_MAX_LESS) {
            currAbsGap_ = currRelGap_ = 0.0;
        }
        else if (feasBound < ALPS_OBJ_MAX_LESS) {
            gapVal = ALPS_MAX(0, feasBound - relBound);
            currAbsGap_ = ALPS_MAX(0, gapVal);
            currRelGap_ = 100 * gapVal / (ALPS_FABS(relBound) + 1.0);
        }

        // print node log
        if ((msgLevel >= 1) && (force||(numNodesProcessed % nodeInterval == 0))) {
            printLog = true;
        }
        if (msgLevel > 200) {
            printLog = true;
        }
        if (printLog) {
#if 0
            int numCols = getNumCols();
            int numRows = getNumRows();
#endif

            // Print header
            if (numNodesProcessed == 0 ||
                (numNodesProcessed % (nodeInterval * 30)) == 0 ||
		msgLevel >= 10) {
                /* Print header. */
                std::cout << std::endl;
                std::cout << "Nodes Done";         /*8 spaces*/
#if 0
                std::cout << "      ObjValue";   /*14 spaces*/
                if (msgLevel > 2) {
                    std::cout << "     Row";    /*8 Spaces*/
                    std::cout << "  Column";    /*8 Spaces*/
                    std::cout << "   Index";     /*8 spaces*/
                    std::cout << "  Parent";     /*8 spaces*/
                    std::cout << "   Depth";     /*8 spaces*/
                }
#endif
                std::cout << "   Upper Bound";
                std::cout << "   Lower Bound";
                std::cout << "      Gap";         /*9 spaces*/
                std::cout << "   Time";
                std::cout << " Nodes Left";
                std::cout << std::endl;
            }

            // Print log
            if (numNodesProcessed < 10000000) {
                printf("%10d", numNodesProcessed);
            }
            else {
                printf("%9dK", numNodesProcessed/1000);
            }

            /* Quality */
#if 0
            if (node->getStatus() == AlpsNodeStatusFathomed) {
                printf("      Fathomed");
            }
            else {
                printf(" %13g", node->getQuality());
            }

            if (msgLevel > 2) {
                if (numRows < 10000000) {
                    printf("%8d", numRows);
                }
                else {
                    printf("%7dK", numRows/1000);
                }
                if (numCols < 10000000) {
                    printf("%8d", numCols);
                }
                else {
                    printf("%7dK", numCols/1000);
                }
                /* This index */
                printf("%8d", node->getIndex());
                /* Paraent index */
                if (node->getParent()) {
                    printf(" %7d", node->getParent()->getIndex());
                }
                else {
                    printf("        ");
                }
                /* Depth */
                printf(" %7d", node->getDepth());
            }
#endif

            if (feasBound > ALPS_OBJ_MAX_LESS) {
                printf("              ");
            }
            else {
                printf(" %13g", feasBound);
            }

            if (relBound > ALPS_OBJ_MAX_LESS) {
                printf("              ");
            }
            else {
                printf(" %13g", relBound);
            }

            /* Gap */
            if (currRelGap_ > ALPS_OBJ_MAX_LESS) {
                printf("         "); /* 9 spaces*/
            }
            else {
                if (currRelGap_ < 1.0e4) {
                    printf(" %7.2f%%", currRelGap_);
                }
                else {
                    printf("% 8g", currAbsGap_);
                }
            }

            int solTime = static_cast<int>(broker_->timer().getCpuTime());
            if (solTime < 1000000) {
                printf("%7d", solTime);
            }
            else {
                solTime = static_cast<int>(solTime/3600.0);
                printf("%6d", solTime);
                printf("H");
            }

            /* Number of left nodes */
            if (numNodesLeft < 10000000) {
                printf(" %10d", numNodesLeft);
            }
            else {
                printf(" %8dK", numNodesLeft/1000);
            }

	    std::cout << std::endl;
	    if (msgLevel >= 10){
		std::cout << std::endl;
	    }
        }
    }
    else if (broker_->getProcType() == AlpsProcessTypeMaster){
        /* For parallel code, only master print node log. */

        double feasBound = broker_->getIncumbentValue();
        double relBound = broker_->getBestEstimateQuality();
        double gapVal = ALPS_OBJ_MAX;
        int numNodeLog = broker_->getNumNodeLog();
        numNodesProcessed = broker_->getNumNodesProcessedSystem();

        /* Compute gap */
        if (relBound > ALPS_OBJ_MAX_LESS) {
            currAbsGap_ = currRelGap_ = 0.0;
        }
        else if (feasBound < ALPS_OBJ_MAX_LESS) {
            gapVal = ALPS_MAX(0, feasBound - relBound);
            currAbsGap_ = ALPS_MAX(0, gapVal);
            currRelGap_ = 100 * gapVal / (ALPS_FABS(relBound) + 1.0);
        }

        if (msgLevel < 1) {
            return;
        }
        if ( (numNodeLog == 0) ||
             (numNodesProcessed - numNodeLog >= nodeInterval) ) {
            printLog = true;
        }

        if (printLog) {

            numNodesLeft = broker_->getNumNodeLeftSystem();

            // Print header
            if (numNodeLog == 0) {
                std::cout << "\n";
                std::cout << "    Node";         /*8 spaces*/
                std::cout << "  BestFeasible";
                std::cout << "     BestBound";
                std::cout << "      Gap";         /*9 spaces*/
                std::cout << "   Time";
                std::cout << "    Left";
                std::cout << std::endl;
            }

            // Print log
            if (numNodesProcessed < 10000000) {
                printf("%8d", numNodesProcessed);
            }
            else {
                printf("%7dK", numNodesProcessed/1000);
            }

            if (feasBound > ALPS_OBJ_MAX_LESS) {
                printf("              ");
            }
            else {
                printf(" %13g", feasBound);
            }

            if (relBound > ALPS_OBJ_MAX_LESS) {
                printf("              ");
            }
            else {
                printf(" %13g", relBound);
            }

            if (currRelGap_ > ALPS_OBJ_MAX_LESS) {
                printf("         "); /* 9 spaces*/
            }
            else {
                if (currRelGap_ < 1.0e4) {
                    printf(" %7.2f%%", currRelGap_);
                }
                else {
                    printf("% 8g", currAbsGap_);
                }
            }

            int solTime = static_cast<int>(broker_->timer().getCpuTime());
            if (solTime < 1000000) {
                printf("%7d", solTime);
            }
            else {
                solTime = static_cast<int>(solTime/3600.0);
                printf("%6d", solTime);
                printf("H");
            }

            /* Number of left nodes */
            if (numNodesLeft < 10000000) {
                printf(" %7d", numNodesLeft);
            }
            else {
                printf(" %6dK", numNodesLeft/1000);
            }

            printf("\n");
            broker_->setNumNodeLog(numNodesProcessed);
        }
    }
}



//#############################################################################

bool
DcoModel::fathomAllNodes()
{
    double feasBound = ALPS_OBJ_MAX;
    double relBound = ALPS_OBJ_MAX;
    double gapVal = ALPS_OBJ_MAX;

    AlpsTreeNode *bestNode = NULL;

    // Compute gap
    feasBound = broker_->getIncumbentValue();
    bestNode = broker_->getBestNode();

    if (bestNode) {
        relBound = bestNode->getQuality();
    }

    if (relBound > ALPS_OBJ_MAX_LESS) {
        currAbsGap_ = currRelGap_ = 0.0;
    }
    else if (feasBound < ALPS_OBJ_MAX_LESS) {
        gapVal = ALPS_MAX(0, feasBound - relBound);
        currAbsGap_ = ALPS_MAX(0, gapVal);
        currRelGap_ = 100 * gapVal / (ALPS_FABS(relBound) + 1.0);
    }
#if 0
    printf("+++ Process %d: currAbsGap_ %g, currRelGap_%g\n",
           broker_->getProcRank(), currAbsGap_,  currRelGap_);
#endif

    if ( (currAbsGap_ <= optimalAbsGap_ + ALPS_ZERO) ||
         (currRelGap_ <= optimalRelGap_ + ALPS_ZERO) ) {
        return true;
    }
    else {
        return false;
    }
}

//#############################################################################

// Borrow from CBC
void
DcoModel::analyzeObjective()
{
    const double *objective = getObjCoef();
    const double *lower = getColLower();
    const double *upper = getColUpper();

    /*------------------------------------------------------
     * Take a first scan to see if there are unfixed continuous
     * variables in the objective.  If so, the minimum objective change
     * could be arbitrarily small. Also pick off the maximum
     * coefficient of an unfixed integer variable.
     * If the objective is found to contain only integer variables,
     * set the fathoming discipline to strict.
     *------------------------------------------------------*/

    double maximumCost = 0.0;
    bool possibleMultiple = true;
    int iColumn;
    int numberColumns = getNumCols();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
	if (upper[iColumn] > lower[iColumn]+1.0e-8) {
	    if (colType_[iColumn] == 'I' || colType_[iColumn] == 'B') {
		maximumCost = CoinMax(maximumCost,
				      fabs(objective[iColumn]));
	    }
	    else if (objective[iColumn]) {
		possibleMultiple = false;
	    }
	}
    }

    //setIntParam(CbcModel::CbcFathomDiscipline,possibleMultiple);

    /*------------------------------------------------------
     * If a nontrivial increment is possible, try and figure
     * it out. We're looking for gcd(c<j>) for all c<j> that
     * are coefficients of unfixed integer variables. Since
     * the c<j> might not be integers, try and inflate them
     * sufficiently that they look like integers (and we'll
     * deflate the gcd later).
     * 2520.0 is used as it is a nice multiple of 2,3,5,7
     *-----------------------------------------------------*/

    if (possibleMultiple && maximumCost) {
	int increment = 0;
	double multiplier = 2520.0;
	while (10.0*multiplier * maximumCost < 1.0e8) {
	    multiplier *= 10.0;
	}
	for (iColumn = 0; iColumn < numberColumns; iColumn++) {
	    if (upper[iColumn] > lower[iColumn]+1.0e-8) {
		if ( (colType_[iColumn] == 'I' || colType_[iColumn] == 'B')
		     && objective[iColumn] ) {
		    double value = fabs(objective[iColumn])*multiplier;
		    int nearest = (int) floor(value+0.5);
		    if (fabs(value-floor(value+0.5)) > 1.0e-8) {
			increment = 0;
			break;
		    }
		    else if (!increment) {
			increment = nearest;
		    }
		    else {
			increment = gcd(increment,nearest);
		    }
		}
	    }
	}

	/*--------------------------------------------------
	 * If the increment beats the current value for objective
	 * change, install it.
	 *--------------------------------------------------*/

	if (increment) {
	    double value = increment;
	    double cutoffInc = DcoPar_->entry(DcoParams::cutoffInc);

	    value /= multiplier;
	    if (value * 0.999 > cutoffInc) {
		if (broker_->getProcRank() == broker_->getMasterRank()) {
		    blisMessageHandler()->message(DISCO_CUTOFF_INC,
						  blisMessages())
			<< value << CoinMessageEol;
		}
		DcoPar_->setEntry(DcoParams::cutoffInc, -value*0.999);
	    }
	}
    }
}

//#############################################################################

void
DcoModel::packSharedPseudocost(AlpsEncoded *encoded, int numToShare)
{
    int k;

    DcoObjectInt *intObj = NULL;
    if (numToShare > 0) {
	// Record how many can be shared.
	encoded->writeRep(numToShare);
	for (k = 0; k < numIntObjects_; ++k) {
	    if (sharedObjectMark_[k]) {
		// Recored which variable.
		encoded->writeRep(k);
		intObj = dynamic_cast<DcoObjectInt*>(objects_[k]);
		(intObj->pseudocost()).encodeTo(encoded);
	    }
	}

	// Clear the mark for next round of sharing.
	clearSharedObjectMark();
    }
    else {
	// numToShare is 0.
	encoded->writeRep(numToShare);
    }

    // Make sure don't exceed large size.
    assert(encoded->size() < broker_->getLargeSize());

}

//#############################################################################

void
DcoModel::unpackSharedPseudocost(AlpsEncoded &encoded)
{
    int k, objIndex, size = 0;

    // Encode and store pseudocost
    DcoObjectInt *intObj = NULL;
    encoded.readRep(size);
    for (k = 0; k < size; ++k) {
        encoded.readRep(objIndex);
        intObj = dynamic_cast<DcoObjectInt *>(objects_[objIndex]);
        (intObj->pseudocost()).decodeFrom(encoded);
    }
}

//#############################################################################

// TODO: Limit msg size
void
DcoModel::packSharedConstraints(AlpsEncoded *encoded)
{
    int numCons = constraintPoolSend_->getNumConstraints();

    if (numCons < DISCO_MIN_SHARE_CON) {
	// Do not send
	int zero = 0; // Let writeRep know it is an integer.
	encoded->writeRep(zero);
	//std::cout << "Don't send " << numCons << " constraints"<< std::endl;
    }
    else {
	// Send constraints.
	if (numCons > DISCO_MAX_SHARE_CON) {
	    numCons = DISCO_MAX_SHARE_CON;
	}
	encoded->writeRep(numCons);
	for (int k = 0; k < numCons; ++k) {
	    AlpsKnowledge *know = constraintPoolSend_->getConstraint(k);
	    know->encode(encoded);
	}
	// Delete all constraints since they are sent.
	constraintPoolSend_->freeGuts();
	//std::cout << "Send " << numCons << " constraints"
	//  << "; msg size = " << encoded->size() << std::endl;
    }
}

//#############################################################################

void
DcoModel::unpackSharedConstraints(AlpsEncoded &encoded)
{
    int numCons = 0;

    encoded.readRep(numCons);
    //std::cout << "Received " << numCons << " constraints"<< std::endl;

    for (int k = 0; k < numCons; ++k) {
	// Unpack and store constraints
	const AlpsKnowledge* know =
            broker_->decoderObject(BcpsKnowledgeTypeConstraint);
	BcpsConstraint *con = NULL;
        con = static_cast<BcpsConstraint *>(know->decode(encoded));
	constraintPoolReceive_->addConstraint(con);
	con = NULL;
    }
}

//#############################################################################

void
DcoModel::packSharedVariables(AlpsEncoded *encoded)
{

}

//#############################################################################

void
DcoModel::unpackSharedVariables(AlpsEncoded &encoded)
{

}

//#############################################################################
