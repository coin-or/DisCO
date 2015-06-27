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
// This file is modified from COIN/Cbc/CbcHeuristic.hpp
//#############################################################################


#ifndef DcoHeuristic_h_
#define DcoHeuristic_h_

#include <string>
#include <vector>

#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

#include "Dco.hpp"

class DcoModel;


//#############################################################################


/** Heuristic base class */
class DcoHeuristic {

 private:

    /** Illegal Assignment operator */
    DcoHeuristic & operator = (const DcoHeuristic& rhs);
    
 protected:

    /** Pointer to the model */
    DcoModel *model_;
    
    /** Heuristics name */
    char *name_;
    
    /** When to call findSolution() routine. 
        DcoHeurStrategyNone:     disable
        DcoHeurStrategyRoot:     just root
        DcoHeurStrategyAuto:     automatically decided by DISCO
        DcoHeurStrategyPeriodic: every 't' nodes
	DcoHeurStrategyBeforeRoot: before solving first LP
    */
    DcoHeurStrategy strategy_;

    /** The frequency with which to call the heuristic */
    int heurCallFrequency_;
    
    /** Number of solutions found. */
    int numSolutions_;

    /** Used CPU/User time. */
    double time_;

    /** The times of calling this heuristic. */
    int calls_;

    /** The times of calling this heuristic and no solutions found. */
    int noSolsCalls_;
    
public:
    
    /** Default Constructor. */
    DcoHeuristic() {
        model_ = NULL;
        name_ = strdup("Unknown");
        strategy_ = DcoHeurStrategyAuto;
	heurCallFrequency_ = 1;
        numSolutions_ = 0;
        time_ = 0.0;
        calls_ = 0;
        noSolsCalls_ = 0;
    }
    
    /** Useful constructor. */
    DcoHeuristic(DcoModel *model, const char *name,
		  DcoHeurStrategy strategy, int heurCallFrequency) {
        model_ = model;
        if (name) {
            name_ = strdup(name);
        }
        else {
            name_ = strdup("Unknown");
        }
        strategy_ = strategy;
	heurCallFrequency_ = heurCallFrequency;
        numSolutions_ = 0;
        time_ = 0.0;
        calls_ = 0;
        noSolsCalls_ = 0;
    }

    /** Distructor. */
    virtual ~DcoHeuristic() { if (name_) free(name_); }

    /** Copy constructor. */
    DcoHeuristic(const DcoHeuristic & rhs) {        
        model_ = rhs.model_;
        name_ = strdup(rhs.name_);
        strategy_ = rhs.strategy_; // What if disabled?
	heurCallFrequency_ = rhs.heurCallFrequency_;
        numSolutions_ = 0;
        time_ = 0.0;
        calls_ = 0;
        noSolsCalls_ = 0;
    }
 
    /** update model (This is needed if cliques update matrix etc). */
    virtual void setModel(DcoModel * model) { model_ = model ;}
    
    /** Get/set strategy. */
    //@{
    virtual void setStrategy(DcoHeurStrategy strategy) {strategy_ = strategy;}
    virtual DcoHeurStrategy strategy() const { return strategy_; }
    //@]

    /** Get/set call frequency. */
    //@{
    virtual void setHeurCallFrequency(int freq) { heurCallFrequency_ = freq; }
    virtual int heurCallFrequency() const { return heurCallFrequency_; }
    //@]

    /** Clone a heuristic. */
    virtual DcoHeuristic * clone() const {
	DcoHeuristic *h = NULL;
	assert(0);
	return h;
    }
    
    /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value 
        This is called after cuts have been added - so can not add cuts
    */
    virtual bool searchSolution(double & objectiveValue, 
                               double * newSolution)=0;
    
    /** returns 0 if no solution, 1 if valid solution, -1 if just
        returning an estimate of best possible solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value 
        (only if nonzero code) 
        This is called at same time as cut generators - so can add cuts
        Default is do nothing
    */
    virtual bool searchSolution(double & objectiveValue,
                               double * newSolution,
                               OsiCuts & cs) { return 0; }
    
    /** return name of generator. */
    inline const char * name() const { return name_; }
    
    /** Record number of solutions found. */
    inline void addNumSolutions(int num=1) { numSolutions_ += num; }

    /** Number of solutions found. */
    inline int numSolutions() const { return numSolutions_; }
    
    /** Record Cpu time used. */
    inline void addTime(double t=0.0) { time_ += t; }

    /** Cpu time used. */
    inline double time() const { return time_; }
    
    /** Record number of times called. */
    inline void addCalls(int c=1) { calls_ += c; }
    
    /** Number of times called. */
    inline int calls() const { return calls_; } 

    /** Number called and no cons found. */
    inline int noSolCalls() const { return noSolsCalls_; }

    /** Increase the number of no cons called. */
    inline void addNoSolCalls(int n=1) { noSolsCalls_ += n; }
};

#endif

