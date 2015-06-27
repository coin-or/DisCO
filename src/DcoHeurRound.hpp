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
// This file is modified from CbcHeuristic.hpp
//#############################################################################

#ifndef DcoHeurRound_h_
#define DcoHeurRound_h_

#include <string>
#include <vector>

#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

#include "DcoHeuristic.hpp"

class DcoModel;

//#############################################################################

/** Rounding Heuristic.  */
class DcoHeurRound : public DcoHeuristic {
 private:
    /** Illegal Assignment operator. */ 
    DcoHeurRound & operator=(const DcoHeurRound& rhs);

 protected:
    /** Column majored matrix. */
    CoinPackedMatrix matrix_;
    
    /** Row majored matrix. */
    CoinPackedMatrix matrixByRow_;
    
    /** Seed for random stuff. */
    int seed_;
    
 public:
    /** Default Constructor. */
    DcoHeurRound() {}
    
    /** Constructor with model - assumed before cuts. */
    DcoHeurRound(DcoModel * model, const char *name,
		  DcoHeurStrategy strategy, int freq)
        :
        DcoHeuristic(model, name, strategy, freq)
        {
            assert(model->solver());
        }

    /** Destructor. */
    ~DcoHeurRound() {}
  
    /** Copy constructor. */
    DcoHeurRound( const DcoHeurRound &);
    
    /** Clone a rounding heuristic */
    virtual DcoHeuristic * clone() const;
    
    /** update model (This is needed if cliques update matrix etc). */
    virtual void setModel(DcoModel * model);

    //using DcoHeuristic::searchSolution ;
    /** returns 0 if no solution, 1 if valid solution
        with better objective value than one passed in
        Sets solution values if good, sets objective value (only if good)
        This is called after cuts have been added - so can not add cuts
    */
    virtual bool searchSolution(double & objectiveValue,
                               double * newSolution);
    
    /** Set seed */
    void setSeed(int value) { seed_ = value; }

};
#endif
