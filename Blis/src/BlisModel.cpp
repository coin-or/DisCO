/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include <cstdio>

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

#include "BlisBranchObjectInt.h"
#include "BlisBranchStrategyMaxInf.h"
#include "BlisBranchStrategyPseudo.h"
#include "BlisBranchStrategyRel.h"
#include "BlisBranchStrategyStrong.h"

#include "BlisConstraint.h"
#include "BlisHeurRound.h"
#include "BlisModel.h"
#include "BlisObjectInt.h"
#include "BlisSolution.h"
#include "BlisTreeNode.h"
#include "BlisVariable.h"

//#############################################################################

void 
BlisModel::init() 
{
    processType_ = AlpsProcessTypeMaster;
    origLpSolver_ = NULL;
    presolvedLpSolver_ = NULL;
    lpSolver_ = NULL;
    numCols_ = 0;
    numRows_ = 0;
    numElems_ = 0;
    colMatrix_ = 0;

    origVarLB_ = NULL;
    origVarUB_ = NULL;
    origConLB_ = NULL;
    origConUB_ = NULL;

    startVarLB_ = NULL;
    startVarUB_ = NULL;
    startConLB_ = NULL;
    startConUB_ = NULL;

    tempVarLBPos_ = NULL;
    tempVarUBPos_ = NULL;
    tempConLBPos_ = NULL;
    tempConUBPos_ = NULL;
    
    objSense_ = 1.0;
    objCoef_ = NULL;

    colType_ = 0;
    numIntObjects_ = 0;
    intColIndices_ = NULL;
    intObjIndices_ = NULL;
    
    presolve_ = new BlisPresolve();
    
    numSolutions_ = 0;
    numHeurSolutions_ = 0;
    
    //savedLpSolution_ = NULL;
    incumbent_ = NULL;
    hotstartStrategy_ = 0;
    
    activeNode_ = NULL;
    numStrong_ = 0;
    
    cutoff_ = COIN_DBL_MAX;
    incObjValue_ = COIN_DBL_MAX;
    blisMessageHandler_ = new CoinMessageHandler();
    blisMessages_ = BlisMessage();

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
    BlisPar_ = new BlisParams;

    /// Timing
    startTime_ = 0.0;
    /// Max solution time allowed
    timeLimit_ = 1.0e75;

    /// Tolerance 
    integerTol_ = 1.0e-5;
    optimalRelGap_ = 1.0e-4;
    optimalAbsGap_ = 1.0e-6;

    /// Heuristic
    heurStrategy_ = true;
    numHeuristics_ = 0;
    heuristics_ = NULL;

    /// Cons related
    cutStrategy_ = 0;
    numCutGenerators_ = 0;
    generators_ = NULL;
    constraintPool_ = NULL;
    oldConstraints_ = NULL;
    oldConstraintsSize_ = 0;
    numOldConstraints_ = 0;
    conRandoms_ = NULL;
}

//#############################################################################

// Read from file (currently MPS format. TODO: LP format).
// Load to solver.
void
BlisModel::readInstance(const char* dataFile)
{
    int j;
    
    int msgLevel =  AlpsPar_->entry(AlpsParams::msgLevel);

    //------------------------------------------------------
    // Read in data from MPS file.
    //------------------------------------------------------
    
    CoinMpsIO *mps = new CoinMpsIO;
    
    int rc = mps->readMps(dataFile, "");
    if(rc) {
        delete mps;
        throw CoinError("Unable to read in instance",
                        "readInstance",
                        "BlisModel");
    }
    mps->messageHandler()->setLogLevel(msgLevel);

    //------------------------------------------------------
    // Get problem data.
    //------------------------------------------------------

    numCols_ = mps->getNumCols();
    numRows_ = mps->getNumRows();
    numElems_ = mps->getNumElements();

    colMatrix_ = new CoinPackedMatrix();    
    *colMatrix_ = *(mps->getMatrixByCol());
    
    origVarLB_ = new double [numCols_];
    origVarUB_ = new double [numCols_];

    origConLB_ = new double [numRows_];
    origConUB_ = new double [numRows_];
    
    memcpy(origVarLB_, mps->getColLower(), sizeof(double) * numCols_);
    memcpy(origVarUB_, mps->getColUpper(), sizeof(double) * numCols_);
    
    memcpy(origConLB_, mps->getRowLower(), sizeof(double) * numRows_);
    memcpy(origConUB_, mps->getRowUpper(), sizeof(double) * numRows_);
    
    //memcpy(startVarLB_, mps->getColLower(), sizeof(double) * numCols_);
    //memcpy(startVarUB_, mps->getColUpper(), sizeof(double) * numCols_);
    
    //memcpy(startConLB_, mps->getRowLower(), sizeof(double) * numRows_);
    //memcpy(startConUB_, mps->getRowUpper(), sizeof(double) * numRows_);
    
    objSense_ = 1.0; /* Default from MPS is minimization */
    
    objCoef_ = new double [numCols_];
    memcpy(objCoef_, mps->getObjCoefficients(), sizeof(double) * numCols_);
    
    //------------------------------------------------------
    // Classify variable type.
    //------------------------------------------------------

    colType_ = new char [numCols_];
    intObjIndices_ = new int [numCols_];
    memset(intObjIndices_, 0, sizeof(int) * numCols_);
    
    intColIndices_ = new int [numCols_];
    numIntObjects_ = 0;
    for(j = 0; j < numCols_; ++j) {
	if (mps->isContinuous(j)) {
	    colType_[j] = 'C';
	}
	else {
	    intColIndices_[numIntObjects_++] = j;
	    if (origVarLB_[j] == 0 && origVarUB_[j] == 1.0) {
		colType_[j] = 'B';
	    }
	    else {
		colType_[j] = 'I';
	    }
	}
    }
    
    //------------------------------------------------------
    // load problem to lp solver.
    //------------------------------------------------------
    
    if (!origLpSolver_) {
        origLpSolver_ = new OsiClpSolverInterface();
    }

    origLpSolver_->loadProblem(*colMatrix_,
			       origVarLB_, origVarUB_,   
			       objCoef_,
			       origConLB_, origConUB_);
    
    origLpSolver_->setObjSense(objSense_);
    origLpSolver_->setInteger(intColIndices_, numIntObjects_);
    
    delete mps;
    
    if (numIntObjects_ == 0) {
        
	// solve lp and throw error.
	//lpSolver_->initialSolve();
	//throw CoinError("The problem does not have integer variables", 
        //	"readInstance", "BlisModel");

        if (broker_->getMsgLevel() > 0) {
            bcpsMessageHandler_->message(BLIS_W_LP, blisMessages())
                << CoinMessageEol;
        }
    }
}

//############################################################################ 

/** Read in Alps parameters. */
void 
BlisModel::readParameters(const int argnum, const char * const * arglist)
{
    //std::cout << "Reading in ALPS parameters ..." << std::endl;
    AlpsPar_->readFromArglist(argnum, arglist);
    //std::cout << "Reading in BLIS parameters ..." << std::endl;
    BlisPar_->readFromArglist(argnum, arglist);
} 

//##############################################################################

/** Write out parameters. */
void 
BlisModel::writeParameters(std::ostream& outstream) const
{
    outstream << "\n================================================"
              <<std::endl;
    outstream << "ALPS Parameters: " << std::endl;
    AlpsPar_->writeToStream(outstream);
    outstream << "\n================================================"
              <<std::endl;
    outstream << "BLIS Parameters: " << std::endl;
    BlisPar_->writeToStream(outstream);
}

//############################################################################ 

/** Do necessary work to make model usable and presolve.
    Return success or not. 
    This function is called when constructing knowledge broker. */
bool 
BlisModel::setupSelf()
{
    int j;

    if (broker_->getMsgLevel() > 0) {

      //std::cout << "**** getProcType = " << broker_->getProcType() << std::endl;
      //if (broker_->getProcType() == AlpsProcessTypeMaster) {
      if (broker_->getProcRank() == broker_->getMasterRank()) {
	bcpsMessageHandler_->message(BCPS_S_VERSION, bcpsMessages())
	  << CoinMessageEol;
	blisMessageHandler()->message(BLIS_S_VERSION, blisMessages())
	  << CoinMessageEol;
      }
    }
    
    //------------------------------------------------------
    // Starting time.
    //------------------------------------------------------
    
    startTime_ = CoinCpuTime();

    processType_ = broker_->getProcType();

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
    
    timeLimit_ = AlpsPar_->entry(AlpsParams::timeLimit);
    
    integerTol_ = BlisPar_->entry(BlisParams::integerTol);
    optimalRelGap_ = BlisPar_->entry(BlisParams::optimalRelGap);
    optimalAbsGap_ = BlisPar_->entry(BlisParams::optimalAbsGap);

    int relibility = BlisPar_->entry(BlisParams::pseudoRelibility);
    
    //------------------------------------------------------
    // Modify parameters.
    //------------------------------------------------------
    
    // Disable Alps message
    //AlpsPar()->setEntry(AlpsParams::msgLevel, 1);
    
    //------------------------------------------------------
    // Create core variables and constraints.
    //------------------------------------------------------

#ifdef BLIS_DEBUG
    std::cout << "setupSelf: numCols_ " << numCols_ 
	      << ", numRows_" << numRows_ 
	      << std::endl;
    std::cout << "Create core ..." << std::endl;
#endif

    numCoreVariables_ = numCols_;
    numCoreConstraints_ = numRows_;
    
    // BlisVariable ** tempVars = new BlisVariable* [numCols_];
    // BlisConstraint ** tempCons = new BlisConstraint* [numRows_];

    coreVariables_ = new BcpsVariable* [numCols_];
    coreConstraints_ = new BcpsConstraint* [numRows_];
    
    for (j = 0; j < numCols_; ++j) {
	BlisVariable * var = new BlisVariable(origVarLB_[j],
					      origVarUB_[j], 
					      origVarLB_[j], 
					      origVarUB_[j]);
	coreVariables_[j] = var;
	var = NULL;
	coreVariables_[j]->setObjectIndex(j);
	coreVariables_[j]->setRepType(BCPS_CORE);
	coreVariables_[j]->setIntType(colType_[j]);
	coreVariables_[j]->setStatus(BCPS_NONREMOVALBE);
    }

    for (j = 0; j < numRows_; ++j) {
        BlisConstraint *con = new BlisConstraint(origConLB_[j], 
                                                 origConUB_[j], 
                                                 origConLB_[j], 
                                                 origConUB_[j]);
        coreConstraints_[j] = con;
        con = NULL;
        coreConstraints_[j]->setObjectIndex(j);
        coreConstraints_[j]->setRepType(BCPS_CORE);
        //coreContraints_[j]->setIntType(colType_[j]);
        coreConstraints_[j]->setStatus(BCPS_NONREMOVALBE);
    }
    
    //------------------------------------------------------
    // Identify integers.
    //------------------------------------------------------
    
    findIntegers(true);
    
    // lpSolver_->initialSolve();
    
#ifdef BLIS_DEBUG_MORE
    std::string problemName;
    lpSolver_->getStrParam(OsiProbName, problemName);
    printf("BLIS: setupSelf: Problem name - %s\n", problemName.c_str());
    lpSolver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
#endif

    //------------------------------------------------------
    // Set branch strategy.
    //------------------------------------------------------

    int brStrategy = BlisPar_->entry(BlisParams::branchStrategy);

    if (brStrategy == BLIS_BS_MAXINFEAS) {
        // Max inf
        branchStrategy_ = new BlisBranchStrategyMaxInf(this);
    }
    else if (brStrategy == BLIS_BS_PSEUDOCOST) {
        // Pseudocost
        branchStrategy_ = new BlisBranchStrategyPseudo(this, 1);
    }
    else if (brStrategy == BLIS_BS_RELIBILITY) {
        // Relibility
        branchStrategy_ = new BlisBranchStrategyRel(this, relibility);
    }
    else if (brStrategy == BLIS_BS_STRONG) {
        // Strong
        branchStrategy_ = new BlisBranchStrategyStrong(this);
    }
    else {
        throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
    }

    brStrategy = BlisPar_->entry(BlisParams::branchStrategyRampUp);

    if (brStrategy == BLIS_BS_MAXINFEAS) {
        // Max inf
      rampUpBranchStrategy_ = new BlisBranchStrategyMaxInf(this);
    }
    else if (brStrategy == BLIS_BS_PSEUDOCOST) {
        // Pseudocost
        rampUpBranchStrategy_ = new BlisBranchStrategyPseudo(this, 1);
    }
    else if (brStrategy == BLIS_BS_RELIBILITY) {
        // Relibility
        rampUpBranchStrategy_ = new BlisBranchStrategyRel(this, relibility);
    }
    else if (brStrategy == BLIS_BS_STRONG) {
        // Strong
        rampUpBranchStrategy_ = new BlisBranchStrategyStrong(this);
    }
    else {
        throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
    }

    //------------------------------------------------------
    // Add heuristics.
    //------------------------------------------------------

    heurStrategy_ = BlisPar_->entry(BlisParams::heurStrategy);
    int useRound = BlisPar_->entry(BlisParams::heurRound); 

    if (useRound == BLIS_NOT_SET) {
        useRound = heurStrategy_;
    }
    if (useRound) {
        // Add rounding heuristic
        BlisHeurRound *heurRound = new BlisHeurRound(this, 
                                                     "Rounding",
                                                     useRound);
        addHeuristic(heurRound);
    }

    // Adjust heurStrategy
    for (j = 0; j < numHeuristics_; ++j) {
        if (heuristics_[j]->strategy() != BLIS_NONE) {
            // Doesn't matter what's the strategy, we just want to 
            // call heuristics.
	  heurStrategy_ = 1;//BLIS_AUTO;
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
    double denseConFactor = BlisPar_->entry(BlisParams::denseConFactor);

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
    denseConCutoff_ = static_cast<int>(aveLen + denseConFactor*stdLen);    
    denseConCutoff_ = ALPS_MIN(numCols_/2, denseConCutoff_);
    denseConCutoff_ = ALPS_MAX(10, denseConCutoff_);
    
#ifdef BLIS_DEBUG
    std::cout << "aveLen=" << aveLen << ", minLen=" << minLen
	      << ", maxLen=" << maxLen << ", stdLen=" << stdLen
	      << ", denseConCutoff_=" << denseConCutoff_ << std::endl;
#endif
    
    // NOTE: maxNumCons is valid only for automatic strategy.
    double cutFactor = BlisPar_->entry(BlisParams::cutFactor);

    maxNumCons_ = (int)((cutFactor - 1.0) * numCoreConstraints_);
    
    constraintPool_ = new BcpsConstraintPool();
    oldConstraints_ = new BlisConstraint* [maxNumCons_];
    oldConstraintsSize_ = maxNumCons_;
    
    cutStrategy_ = BlisPar_->entry(BlisParams::cutStrategy); 

#ifdef BLIS_DEBUG
    std::cout << "cutStrategy_ = " << cutStrategy_ << std::endl;
#endif

    int clique = BlisPar_->entry(BlisParams::cutClique);
    int fCover = BlisPar_->entry(BlisParams::cutFlowCover);
    int gomory = BlisPar_->entry(BlisParams::cutGomory); 
    int knap = BlisPar_->entry(BlisParams::cutKnapsack); 
    int mir = BlisPar_->entry(BlisParams::cutMir); 
    int oddHole = BlisPar_->entry(BlisParams::cutOddHole);
    int probe = BlisPar_->entry(BlisParams::cutProbing);
    int twoMir = BlisPar_->entry(BlisParams::cutTwoMir); 

    //------------------------------------------------------
    // Add cut generators.
    //------------------------------------------------------

    if (probe == BLIS_NOT_SET) {
        // Only at root by default
        probe = BLIS_ROOT;
    }
    if (probe != BLIS_NONE) {
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
        addCutGenerator(probing, "Probing", probe);
    }

    if (clique == BLIS_NOT_SET) {
        // Only at root by default
        clique = cutStrategy_;
    }
    if (clique != BLIS_NONE) {
        CglClique *cliqueCut = new CglClique ;
        cliqueCut->setStarCliqueReport(false);
        cliqueCut->setRowCliqueReport(false);
        addCutGenerator(cliqueCut, "Clique", clique);
    }

    if (oddHole == BLIS_NOT_SET) {
        // Disable by default
        oddHole = BLIS_NONE;
    }
    if (oddHole != BLIS_NONE) {
        CglOddHole *oldHoleCut = new CglOddHole;
        oldHoleCut->setMinimumViolation(0.005);
        oldHoleCut->setMinimumViolationPer(0.00002);
        // try larger limit
        oldHoleCut->setMaximumEntries(200);
        addCutGenerator(oldHoleCut, "OddHole", oddHole);
    }

    if (fCover == BLIS_NOT_SET) {
         fCover = cutStrategy_;
    }
    if (fCover != BLIS_NONE) {
        CglFlowCover *flowGen = new CglFlowCover;
        addCutGenerator(flowGen, "Flow Cover", fCover);
    }

    if (knap == BLIS_NOT_SET) {
        knap = cutStrategy_;
    }
    if (knap != BLIS_NONE) {
        CglKnapsackCover *knapCut = new CglKnapsackCover;
        addCutGenerator(knapCut, "Knapsack", knap);
    }

    if (mir == BLIS_NOT_SET) {
        // Disable by default
        twoMir = BLIS_NONE;
    }
    if (mir != BLIS_NONE) {
        CglMixedIntegerRounding2 *mixedGen = new CglMixedIntegerRounding2;
        addCutGenerator(mixedGen, "MIR", mir);
    }

    if (gomory == BLIS_NOT_SET) {
        gomory = cutStrategy_;
    }
    if (gomory != BLIS_NONE) {
        CglGomory *gomoryCut = new CglGomory;
        // try larger limit
        gomoryCut->setLimit(300);
        addCutGenerator(gomoryCut, "Gomory", gomory);
    }

    if (twoMir == BLIS_NOT_SET) {
        // Disable by default
        twoMir = BLIS_NONE;
    }
    if (twoMir != BLIS_NONE) {
        CglTwomir *twoMirCut =  new CglTwomir;
        addCutGenerator(twoMirCut, "Two MIR", twoMir);
    }
    
    // Random vector
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

    // Adjust cutstrategy
    for (j = 0; j < numCutGenerators_; ++j) {
        int strategy = cutGenerators(j)->strategy();
        if (strategy != BLIS_NONE) {
            // Doesn't matter what's the strategy, we just want to 
            // Generate cuts.
            cutStrategy_ =  BLIS_AUTO;
            break;
        }
    }
     
#ifdef BLIS_DEBUG_MORE
    std::cout << "AFTER: cutStrategy_ = " << cutStrategy_ << std::endl;
#endif

    return true;
}

//############################################################################ 

void 
BlisModel::preprocess()
{
    int numPasses = 10;
    double feaTol = 1.0e-6;
    bool keepIntegers = true;
    char *prohibited = 0;

    bool doPresolve = BlisPar_->entry(BlisParams::presolve);
    
    //std::cout << "Presolve = "<< doPresolve << std::endl;
    
    doPresolve = false;
    
    if (doPresolve) {
	presolvedLpSolver_ = presolve_->preprocess(*lpSolver_,
						   feaTol,
						   keepIntegers,
						   numPasses,
						   prohibited);
	lpSolver_ = presolvedLpSolver_->clone();
    }
    else {
	lpSolver_ = origLpSolver_;
    }
}

//############################################################################ 

void 
BlisModel::postprocess()
{
    
}

//############################################################################ 

bool 
BlisModel::feasibleSolution(int & numIntegerInfs)
{
    bool feasible = true;
    numIntegerInfs = 0;
    int i = -1;
    //const int numCols = lpSolver_->getNumCols();
    const double *savedLpSolution = lpSolver_->getColSolution();

#if 0
    if (savedLpSolution_ != 0) {
      delete [] savedLpSolution_;
      savedLpSolution_ = 0;
    }
    
    savedLpSolution_ = new double [numCols];
    memcpy(savedLpSolution_, 
	   lpSolver_->getColSolution(), 
	   sizeof(double) * numCols);
#endif
    
    for (i = 0; i < numIntObjects_; ++i) {
      if ( ! checkInteger(savedLpSolution[intColIndices_[i]]) ) {
	++numIntegerInfs;
	feasible = false;
      }
    }

    return feasible;
}

//############################################################################ 

bool
BlisModel::setBestSolution(BLIS_SOL_TYPE how,
			   double & objectiveValue, 
			   const double * solution, 
			   bool fixVariables)
{
    double cutoff = getCutoff();
    
    // Double check the solution to catch pretenders.
    if (objectiveValue >= cutoff) {  // Bad news
        return false;
    }
    else {  // Better solution
	incObjValue_ = objectiveValue;

	int numColumns = lpSolver_->getNumCols();
	if (incumbent_ == 0) {
	    incumbent_ = new double[numColumns];
	}
	
	memcpy(incumbent_, solution, numColumns*sizeof(double));

        // Update cutoff value in lp solver.
	setCutoff(incObjValue_);
        ++numSolutions_;

	switch (how) {
	case BLIS_SOL_BOUNDING:
#ifdef BLIS_DEBUG
	  std::cout << "Rounding heuristics found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_BRANCHING:
#ifdef BLIS_DEBUG
	  std::cout << "Branching found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_DIVING:
            ++numHeurSolutions_;
#ifdef BLIS_DEBUG
	  std::cout << "Branching found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_ROUNDING:
            ++numHeurSolutions_;
#ifdef BLIS_DEBUG
	  std::cout << "Rounding heuristics found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	case BLIS_SOL_STRONG:
#ifdef BLIS_DEBUG
	  std::cout << "Strong branching found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	default:
#ifdef BLIS_DEBUG
	  std::cout << "Nowhere found a better solution" 
		    <<", old cutoff = " << cutoff 
		    << ", new cutoff = " << getCutoff()  << std::endl;
#endif
	  break;
	}


	//setMinimizationObjValue(objectiveValue * lpSolver_->getObjSense());
	
	return true;
    }
}

//############################################################################ 

void 
BlisModel::findIntegers(bool startAgain)
{
    assert(lpSolver_);

    if (numIntObjects_ && !startAgain && objects_) return;

    int iCol;    
    int numCols = getNumCols();
    
    const double *colLB = lpSolver_->getColLower();
    const double *colUB = lpSolver_->getColUpper();
    BlisObjectInt *intObject = NULL;

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

    double weight = BlisPar_->entry(BlisParams::pseudoWeight);
    
    int numObjects = 0;
    int iObject;
    BcpsObject ** oldObject = objects_;    

    for (iObject = 0; iObject < numObjects_; ++iObject) {
	BlisObjectInt * obj =
	    dynamic_cast <BlisObjectInt *>(oldObject[iObject]) ;
        
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
            
	    intObject = new BlisObjectInt(numIntObjects_,
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
BlisModel::resolve()
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
BlisModel::deleteObjects()
{
    delete [] priority_;
    priority_ = NULL;
    int i;
    for (i = 0; i < numObjects_; ++i) delete objects_[i];
    delete [] objects_;
    objects_ = NULL;
    numObjects_ = 0;
    findIntegers(true);
}

//#############################################################################

BlisModel::~BlisModel()
{
    gutsOfDestructor();
}

//#############################################################################

void 
BlisModel::gutsOfDestructor()
{
    int i;
    bool doPresolve = BlisPar_->entry(BlisParams::presolve);
    

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

    delete [] origVarLB_;
    delete [] origVarUB_;

    delete [] origConLB_;
    delete [] origConUB_;

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

    delete presolve_;
    
    if (numHeuristics_ > 0) {
#ifdef BLIS_DEBUG
        std::cout << "MODEL: distructor: numHeuristics =  " 
                  << numHeuristics_ << std::endl;
#endif
        BlisHeuristic *tempH = NULL;
        
        for (i = 0; i < numHeuristics_; ++i) {
            tempH = heuristics_[i];
            delete tempH;
        }
        
        delete [] heuristics_;
        heuristics_ = NULL;
    }

    if (generators_ != NULL) {
#ifdef BLIS_DEBUG
        std::cout << "MODEL: distructor: numCutGenerators = " 
                  << numCutGenerators_ << std::endl;
#endif
        BlisConGenerator *temp = NULL;
        for (i = 0; i < numCutGenerators_; ++i) {
            temp = generators_[i];
            delete temp;
        }
        delete [] generators_;
        generators_ = NULL;
    }

    delete constraintPool_;
    delete [] oldConstraints_;
    delete branchStrategy_;
    delete rampUpBranchStrategy_;

    delete [] conRandoms_;
    
    delete BlisPar_;
    delete blisMessageHandler_;
    if (doPresolve) {
	delete presolvedLpSolver_;
	delete lpSolver_;
    }
}

//#############################################################################

bool 
BlisModel::feasibleSolution(int & numIntegerInfs, int & numObjectInfs)
{
    int numUnsatisfied = 0;
    double sumUnsatisfied = 0.0;
    int preferredWay;
    int j;
    
#if 0
    if (savedLpSolution_ != 0) {
        delete [] savedLpSolution_;
        savedLpSolution_ = 0;
    }
    
    savedLpSolution_ = new double [lpSolver_->getNumCols()];
    
    // Put current solution in safe place
    memcpy(savedLpSolution_, 
           lpSolver_->getColSolution(),
	   lpSolver_->getNumCols() * sizeof(double));
#endif

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

    //printf("numUnsatisfied = %d\n",numUnsatisfied);
    
    return (!numUnsatisfied);
}

//#############################################################################
/*
  Set branching priorities.

  Setting integer priorities looks pretty robust; the call to findIntegers
  makes sure that integer objects are in place. Setting priorities for
  other objects is entirely dependent on their existence, and the routine may
  quietly fail in several directions.
*/
void 
BlisModel::passInPriorities (const int * priorities,
			     bool ifObject, 
			     int defaultValue)
{

    // FIXME: not completed.
    int i;

    findIntegers(false);

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

AlpsTreeNode * 
BlisModel::createRoot() {
  
  //-------------------------------------------------------------
  // NOTE: Root will be deleted by ALPS. Root is an explicit node.
  //-------------------------------------------------------------
    
  BlisTreeNode* root = new BlisTreeNode;
  BlisNodeDesc* desc = new BlisNodeDesc(this);
  root->setDesc(desc);

  //-------------------------------------------------------------
  // NOTE: Although original data are stored in model when reading. 
  //   Root desc still store a full copy of col and row bounds when creating.
  //   It will store soft differences after finding a branching object. 
  //   The soft difference are due to reduced cost fixing and probing.
  //   Also the added cols and rows will be stored.
  //-------------------------------------------------------------
  int k;

  BcpsVariable ** vars = getCoreVariables();
  BcpsConstraint ** cons = getCoreConstraints();

#ifdef BLIS_DEBUG
  std::cout << "BLIS: createRoot(): numCoreVariables_=" << numCoreVariables_
	    << ", numCoreConstraints_=" << numCoreConstraints_ << std::endl;
#endif  

  int *varIndices1 = new int [numCoreVariables_];
  int *varIndices2 = new int [numCoreVariables_];
  int *varIndices3 = NULL; //new int [numCoreVariables_];
  int *varIndices4 = NULL; //new int [numCoreVariables_];
  double *vlhe = new double [numCoreVariables_];
  double *vuhe = new double [numCoreVariables_];
  double *vlse = NULL; //new double [numCoreVariables_];
  double *vuse = NULL; //new double [numCoreVariables_];

  int *conIndices1 = new int [numCoreConstraints_];
  int *conIndices2 = new int [numCoreConstraints_];
  int *conIndices3 = NULL; //new int [numCoreConstraints_];
  int *conIndices4 = NULL; //new int [numCoreConstraints_];
  double *clhe = new double [numCoreConstraints_];
  double *cuhe = new double [numCoreConstraints_];
  double *clse = NULL; //new double [numCoreConstraints_];
  double *cuse = NULL; //new double [numCoreConstraints_];

  //-------------------------------------------------------------
  // Get var bounds and indices.
  //-------------------------------------------------------------

  for (k = 0; k < numCoreVariables_; ++k) {
    vlhe[k] = vars[k]->getLbHard();
    vuhe[k] = vars[k]->getUbHard();
    //vlse[k] = vars[k]->getLbSoft();
    //vuse[k] = vars[k]->getUbSoft();
    varIndices1[k] = k;
    varIndices2[k] = k;
    
    //varIndices3[k] = k;
    //varIndices4[k] = k;

#ifdef BLIS_DEBUG_MORE
    std::cout << "BLIS: createRoot(): var "<< k << ": hard: lb=" << vlhe[k]
	      << ", ub=" << vuhe[k] << std::endl;
#endif  
    
  }

  //-------------------------------------------------------------  
  // Get con bounds and indices.
  //-------------------------------------------------------------

  for (k = 0; k < numCoreConstraints_; ++k) {
    clhe[k] = cons[k]->getLbHard();
    cuhe[k] = cons[k]->getUbHard();
    //clse[k] = cons[k]->getLbSoft();
    //cuse[k] = cons[k]->getUbSoft();
    conIndices1[k] = k;
    conIndices2[k] = k;
    //conIndices3[k] = k;
    //conIndices4[k] = k;
  }

  int *tempInd = NULL;
  BcpsObject **tempObj = NULL;
  
  desc->assignVars(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numCoreVariables_, varIndices1, vlhe, /*Var hard lb*/
		   false, numCoreVariables_, varIndices2, vuhe, /*Var hard ub*/
		   false, 0, varIndices3, vlse, /*Var soft lb*/
		   false, 0, varIndices4, vuse);/*Var soft ub*/
  desc->assignCons(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numCoreConstraints_,conIndices1,clhe, /*Con hard lb*/
		   false, numCoreConstraints_,conIndices2,cuhe, /*Con hard ub*/
		   false, 0,conIndices3,clse, /*Con soft lb*/
		   false, 0,conIndices4,cuse);/*Con soft ub*/

  //-------------------------------------------------------------  
  // Mark it as an explicit node.
  //-------------------------------------------------------------

  root->setExplicit(1);
  
  return root;
}

//#############################################################################

void 
BlisModel::addHeuristic(BlisHeuristic * heuristic)
{
    BlisHeuristic ** temp = heuristics_;
    heuristics_ = new BlisHeuristic * [numHeuristics_ + 1];

    memcpy(heuristics_, temp, numHeuristics_ * sizeof(BlisHeuristic *));
    delete [] temp;
    
    heuristics_[numHeuristics_++] = heuristic;
}

//#############################################################################

void 
BlisModel::addCutGenerator(CglCutGenerator * generator,
			   const char * name,
			   int strategy,
			   bool normal, 
			   bool atSolution,
			   bool whenInfeasible)
{

#if 0
    if (!generators_) {        
        generators_ = new BlisCutGenerator* [50];
    }
    generators_[numCutGenerators_++] = new BlisConGenerator(this, generator, 
                                                            strategy, name,
                                                            normal, atSolution,
                                                            whenInfeasible);
#else
    BlisConGenerator ** temp = generators_;

    generators_ = new BlisConGenerator * [(numCutGenerators_ + 1)];
    memcpy(generators_, temp, numCutGenerators_ * sizeof(BlisConGenerator *));
    
    generators_[numCutGenerators_++] = 
        new BlisConGenerator(this, generator, name, strategy,
                             normal, atSolution, whenInfeasible);
    delete [] temp;
    temp = NULL;
#endif

}

//#############################################################################

AlpsEncoded* 
BlisModel::encode() const 
{ 
    AlpsReturnCode status = ALPS_OK;

    // NOTE: "ALPS_MODEL" is the type name.
    AlpsEncoded* encoded = new AlpsEncoded("ALPS_MODEL");

    //------------------------------------------------------
    // Encode Alps part. 
    // NOTE: Nothing to do for Bcps part.
    //------------------------------------------------------

    status = encodeAlps(encoded);
    
    //------------------------------------------------------
    // Encode Blis part. 
    //------------------------------------------------------

    //------------------------------------------------------
    // Blis parameter.
    //------------------------------------------------------

    BlisPar_->pack(*encoded);

    //------------------------------------------------------
    // Get a column matrix.
    //------------------------------------------------------
    
    const CoinPackedMatrix *matrixByCol = lpSolver_->getMatrixByCol();
    
    const int numRows = lpSolver_->getNumRows();
    encoded->writeRep(numRows);
    
    const int numCols = lpSolver_->getNumCols();
    encoded->writeRep(numCols);

#ifdef BLIS_DEBUG
    std::cout << "BlisModel::encode()-- numRows="<< numRows << "; numCols=" 
	      << numCols << std::endl;
#endif

    //------------------------------------------------------
    // Variable bounds.
    //------------------------------------------------------

    const double* collb = lpSolver_->getColLower();
    // NOTE: when write a array to buffer, the length is written
    //       before the actual array. So don't need send numCols or
    //       numRows.
    encoded->writeRep(collb, numCols);
    const double* colub = lpSolver_->getColUpper();
    encoded->writeRep(colub, numCols);

    //------------------------------------------------------
    // Objective.
    //------------------------------------------------------

    const double* obj = lpSolver_->getObjCoefficients();
    encoded->writeRep(obj, numCols);
    const double objSense = lpSolver_->getObjSense();
    encoded->writeRep(objSense);

    //------------------------------------------------------
    // Constraint bounds.
    //------------------------------------------------------

    const double* rowlb = lpSolver_->getRowLower();
    encoded->writeRep(rowlb, numRows);
    const double* rowub = lpSolver_->getRowUpper();
    encoded->writeRep(rowub, numRows);

    //------------------------------------------------------
    // Matrix.
    //------------------------------------------------------

    int numElements = lpSolver_->getNumElements();
    encoded->writeRep(numElements);
    const double* elementValue = matrixByCol->getElements();
    encoded->writeRep(elementValue, numElements);

    const CoinBigIndex* colStart = matrixByCol->getVectorStarts();

    int numStart = numCols + 1;
    encoded->writeRep(colStart, numStart);

    const int* index = matrixByCol->getIndices();
    encoded->writeRep(index, numElements);

    //------------------------------------------------------
    // Variable type.
    //------------------------------------------------------

    encoded->writeRep(numIntObjects_);
    encoded->writeRep(intColIndices_, numIntObjects_);

    //------------------------------------------------------
    // Debug.
    //------------------------------------------------------

#if 0
    std::cout << "BlisModel::encode()-- objSense="<< objSense
	      << "; numElements="<< numElements 
	      << "; numIntObjects_=" << numIntObjects_ 
	      << "; numStart = " << numStart <<std::endl;
#endif

#if 0
    std::cout << "rowub=";
    for (int i = 0; i < numRows; ++i){
	std::cout <<rowub[i]<<" ";
    }
    std::cout << std::endl;
    std::cout << "elementValue=";
    for (int j = 0; j < numElements; ++j) {
	std::cout << elementValue[j] << " ";
    }
    std::cout << std::endl;    
#endif

    return encoded;
}

//#############################################################################

void
BlisModel::decodeToSelf(AlpsEncoded& encoded) 
{
    AlpsReturnCode status = ALPS_OK;

    //------------------------------------------------------
    // Decode Alps part. 
    // NOTE: Nothing to do for Bcps part.
    //------------------------------------------------------

    status = decodeAlps(encoded);

    //------------------------------------------------------
    // Decode Blis part. 
    //------------------------------------------------------


    //------------------------------------------------------
    // Blis Parameters.
    //------------------------------------------------------

    BlisPar_->unpack(encoded);

    encoded.readRep(numRows_);
    encoded.readRep(numCols_);    

#ifdef BLIS_DEBUG
    std::cout << "BlisModel::decode()-- numRows_="<< numRows_ 
	      << "; numCols_=" << numCols_ << std::endl;
#endif

    //------------------------------------------------------
    // Variable bounds.
    //------------------------------------------------------

    int numCols;
    encoded.readRep(origVarLB_, numCols);
    assert(numCols == numCols_);

    encoded.readRep(origVarUB_, numCols);
    assert(numCols == numCols_);
    
    //------------------------------------------------------
    // Objective.
    //------------------------------------------------------

    encoded.readRep(objCoef_, numCols);
    assert(numCols == numCols_);
    
    encoded.readRep(objSense_);
    
    //------------------------------------------------------
    // Constraint bounds.
    //------------------------------------------------------
    
    int numRows;
    
    encoded.readRep(origConLB_, numRows);
    assert(numRows == numRows_);

    encoded.readRep(origConUB_, numRows);
    assert(numRows == numRows_);
    
    //------------------------------------------------------
    // Matrix.
    //------------------------------------------------------

    encoded.readRep(numElems_);

    int numElements;
    double* elementValue = NULL;
    encoded.readRep(elementValue, numElements);
    assert(numElements == numElems_);

    CoinBigIndex* colStart;

    int numStart;
    encoded.readRep(colStart, numStart);
    assert(numStart == numCols_ + 1);

    int* index;
    encoded.readRep(index, numElements);
    assert(numElements == numElems_);

    //------------------------------------------------------
    // Variable type.
    //------------------------------------------------------

    encoded.readRep(numIntObjects_);
    int numInts;
    encoded.readRep(intColIndices_, numInts);
    assert(numInts == numIntObjects_);

    //------------------------------------------------
    // Classify variable type
    //------------------------------------------------

    colType_ = new char [numCols_];
    int j, ind = -1;
    for(j = 0; j < numCols_; ++j) {
	colType_[j] = 'C';
    }
    for(j = 0; j < numInts; ++j) {
	ind = intColIndices_[j];
	assert(ind >= 0 && ind < numCols_);
	if (origVarLB_[ind] == 0.0 && origVarUB_[ind] == 1.0) {
	    colType_[ind] = 'B';
	}
	else {
	    colType_[ind] = 'I';
	}
    }
    
    //------------------------------------------------------
    // Debug.
    //------------------------------------------------------

#if 0
    std::cout << "BlisModel::decode()-- objSense="<< objSense_
	      <<  "; numElements="<< numElements 
	      << "; numberIntegers_=" << numIntObjects_ 
	      << "; numStart = " << numStart <<std::endl;
#endif

#if 0
    int i;
    std::cout << "origConUB_=";
    for (i = 0; i < numRows; ++i){
	std::cout <<origConUB_[i]<<" ";
    }
    std::cout << std::endl;
    std::cout << "elementValue=";
    for (j = 0; j < numElements; ++j) {
	std::cout << elementValue[j] << " ";
    }
    std::cout << std::endl;  
    std::cout << "index=";
    for (j = 0; j < numElements; ++j) {
	std::cout << index[j] << " ";
    }
    std::cout << std::endl;  
    std::cout << "colStart=";
    for (j = 0; j < numCols + 1; ++j) {
	std::cout << colStart[j] << " ";
    }
    std::cout << std::endl;   
#endif

    //------------------------------------------------------
    // Check if lpSolver_ is declared in main.
    // decodeToSelf is called by processes other than master.
    //------------------------------------------------------

    lpSolver_ = origLpSolver_;
    assert(lpSolver_);

    //------------------------------------------------------
    // Load data to lp solver.
    //------------------------------------------------------

    lpSolver_->loadProblem(numCols, numRows,
			   colStart, index, elementValue,
			   origVarLB_, origVarUB_, 
			   objCoef_,
			   origConLB_, origConUB_);
    
    lpSolver_->setObjSense(objSense_);
    lpSolver_->setInteger(intColIndices_, numIntObjects_);

    //------------------------------------------------------
    // Clean up.
    //------------------------------------------------------

    delete [] colStart;
    colStart = NULL;

    delete [] index;
    index = NULL;

    delete [] elementValue;
    elementValue = NULL;
}

//#############################################################################

AlpsEncoded* 
BlisModel::encodeKnowlegeShared()
{
    AlpsEncoded* encoded = 0;

    int k = 0;
    int size = 0;
    int numShared = 0;
    int frequency = -1, depth = -1;
    
    bool sharePseudo = true;
    bool shareCon = false;
    bool shareVar = false;
    bool share = false;
    
    int phase = broker_->getPhase();
    
    if (phase == ALPS_PHASE_RAMPUP) {
        sharePseudo = BlisPar_->entry(BlisParams::sharePseudocostRampUp);
    }
    else if (phase == ALPS_PHASE_SEARCH) {
        sharePseudo = BlisPar_->entry(BlisParams::sharePseudocostSearch);
        if (sharePseudo) {
            // Depth and frequency
            depth = BlisPar_->entry(BlisParams::sharePcostDepth);
            frequency =  BlisPar_->entry(BlisParams::sharePcostFrequency);
            if ( (numNodes_ % frequency != 0) || 
                 (broker_->getTreeDepth() >  depth) ) {
                sharePseudo = false;
            }
        }
    }

    if (sharePseudo) {
        for (k = 0; k < numIntObjects_; ++k) {
            if (sharedObjectMark_[k]){
                ++numShared;
            }
        }
        if (numShared) share = true;
    }
    // TODO: cuts, etc.

#if 0
    std::cout << "++++ sharePseudo =  " << sharePseudo
              << ", shareCon = " << shareCon
              << ", shareVar = " << shareVar
              << ", numShared = " << numShared 
              << ", numNodes_ = " << numNodes_
              << ", frequency = " << frequency
              << ", depth = " << depth
              << ", treeDepth = " << broker_->getTreeDepth()
              << ", share = " << share << std::endl;
#endif        

    if (share) {
        // NOTE: "ALPS_MODEL_GEN" is the type name. We don't need to
        //       register it since ALPS_MODEL is registered.
        
        BlisObjectInt *intObj = NULL;
        
        if (numShared > 0) {
            encoded = new AlpsEncoded("ALPS_MODEL_GEN");
            // Record how many can be shared.
            encoded->writeRep(numShared);
            for (k = 0; k < numIntObjects_; ++k) {
                if (sharedObjectMark_[k]) {
                    // Recored which object.
                    encoded->writeRep(k);
                    intObj = dynamic_cast<BlisObjectInt*>(objects_[k]);
                    (intObj->pseudocost()).encodeTo(encoded);
                }    
            }
            
            // Clear the mark for next round of sharing.
            clearSharedObjectMark();
        }
        else {
            encoded->writeRep(numShared);
        }
        
	// Make sure don't exceed large size.
	assert(encoded->size() < broker_->getLargeSize());
    }

    return encoded;
}

//#############################################################################

void 
BlisModel::decodeKnowledgeShared(AlpsEncoded& encoded)
{
    int k, objIndex, size = 0;
    
    // Encode and store pseudocost
    BlisObjectInt *intObj = NULL;
    encoded.readRep(size);
    for (k = 0; k < size; ++k) {
        encoded.readRep(objIndex);
        intObj = dynamic_cast<BlisObjectInt *>(objects_[objIndex]);
        (intObj->pseudocost()).decodeFrom(encoded);
    }

#if 0        
    // Encode and store generated constraints that should be shared.
    encoded.readRep(size);
    for (k = 0; k < size; ++k) {
        // TODO
        assert(0);
    }
    
    // Encode and store variables that should be shared.
    encoded.readRep(size);
    for (k = 0; k < size; ++k) {
        // TODO
        assert(0);
    }
#endif
}

//#############################################################################

/** Register knowledge. */
void 
BlisModel::registerKnowledge() {
    // Register model, solution, and tree node
    assert(broker_);
    broker_->registerClass("ALPS_MODEL", new BlisModel);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Alps model." << std::endl;
    }
    
    broker_->registerClass("ALPS_NODE", new BlisTreeNode(this));
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Alps node." << std::endl;
    }
    
    broker_->registerClass("ALPS_SOLUTION", new BlisSolution);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Alps solution." << std::endl;
    }
    
    broker_->registerClass("BCPS_CONSTRAINT", new BlisConstraint);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Bcps constraint." << std::endl;
    }
    
    broker_->registerClass("BCPS_VARIABLE", new BlisVariable);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Bcps variable." << std::endl;
    }
}

//#############################################################################

/** Log of specific models. */
void 
BlisModel::modelLog() 
{
    if (processType_ != AlpsProcessTypeMaster) return;
    
    int logFileLevel = AlpsPar_->entry(AlpsParams::logFileLevel);
    int msgLevel = AlpsPar_->entry(AlpsParams::msgLevel);
    if (logFileLevel > 0) {
	std::string logfile = AlpsPar_->entry(AlpsParams::logFile);
	std::ofstream logFout(logfile.c_str(), std::ofstream::app);
	writeParameters(logFout);
    }

    if (msgLevel > 0) {
        int k;
        for (k = 0; k < numCutGenerators_; ++k) {
            if (cutGenerators(k)->calls() > 0) {
                blisMessageHandler()->message(BLIS_CUT_STAT_FINAL,
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
                blisMessageHandler()->message(BLIS_HEUR_STAT_FINAL,
                                              blisMessages())
                    << heuristics(k)->name()
                    << heuristics(k)->calls()
                    << heuristics(k)->numSolutions()
                    << heuristics(k)->time()
                    << heuristics(k)->strategy()
                    << CoinMessageEol;
            }   
        }
    }    
}

//#############################################################################

void 
BlisModel::nodeLog(AlpsTreeNode *node, bool force) 
{
    int nodeInterval = 
	broker_->getModel()->AlpsPar()->entry(AlpsParams::nodeLogInterval);
    
    int numNodesProcessed = broker_->getNumNodesProcessed();
    int numNodesLeft = broker_->updateNumNodesLeft();
    int msgLevel = broker_->getMsgLevel();
    bool printLog = false;
    
    //std::cout << "nodeInterval = " << nodeInterval << std::endl;
    
    AlpsTreeNode *bestNode = NULL;

    if ((msgLevel > 1) && (force||(numNodesProcessed % nodeInterval == 0))) {
        printLog = true;
    }

    if (broker_->getProcType() != AlpsProcessTypeMaster) {
        printLog = false;
    }
    
    if (msgLevel > 200) {
        printLog = true;
    }

#if 0
    std::cout << "==== Process " << broker_->getProcRank()
              << ": printLog = " << printLog 
              << ", msgLevel = " << msgLevel 
              << ", proc type = " << broker_->getProcType()
              << std::endl;
#endif
    
    if (printLog) {
        double feasBound = ALPS_OBJ_MAX;
	double relBound = ALPS_OBJ_MAX;
	double gap = ALPS_OBJ_MAX;
	double gapVal = ALPS_OBJ_MAX;
	
        if (broker_->getNumKnowledges(ALPS_SOLUTION) > 0) {
            feasBound = (broker_->getBestKnowledge(ALPS_SOLUTION)).second;
        }
	
        bestNode = broker_->getBestNode();
        
        if (bestNode) {
            relBound = bestNode->getQuality();
        }	

	if (numNodesProcessed == 0 ||
	    (numNodesProcessed % (nodeInterval * 30)) == 0) {
	    /* Print header. */
	    std::cout << "\n";
            std::cout << "    Node";         /*8 spaces*/
            std::cout << "      ObjValue";   /*14 spaces*/
            if (msgLevel > 2) {
                std::cout << "   Index";     /*8 spaces*/
                std::cout << "  Parent";     /*8 spaces*/
                std::cout << "   Depth";     /*8 spaces*/
            }

            std::cout << "  BestFeasible";
            std::cout << "     BestBound";
            std::cout << "      Gap";         /*9 spaces*/
            std::cout << "   Time";
            std::cout << "    Left";
            std::cout << std::endl;
	}
	
	if (numNodesProcessed < 10000000) {
	    printf("%8d", numNodesProcessed);
	}
	else {
	    printf("%7dK", numNodesProcessed/1000);
	}

        /* Quality */
	if (node->getStatus() == AlpsNodeStatusFathomed) {
	    printf("      Fathomed");
	}
	else {
	    printf(" %13g", node->getQuality());
	}

        if (msgLevel > 2) {
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
	if ( (feasBound < ALPS_OBJ_MAX_LESS) &&
	     (relBound < ALPS_OBJ_MAX_LESS) ) {
	    gapVal = ALPS_MAX(0, feasBound - relBound);
	    gap = 100 * gapVal / (ALPS_FABS(relBound) + 1.0);
	}
	if (gap > ALPS_OBJ_MAX_LESS) {
	    printf("         "); /* 9 spaces*/
 	}
	else {
	    if (gap < 1.0e4) {
		printf(" %7.2f%%", gap);
	    }
	    else {
		printf("% 8g", gapVal);
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
    }
}

//#############################################################################
