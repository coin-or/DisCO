/*===========================================================================*
 * This file is part of a solver for the Vehicle Routing Problem             *
 * developed using the BiCePS Linear Integer Solver (BLIS).                  *
 *                                                                           *
 * This solver is distributed under the Common Public License as part of     * 
 * the COIN-OR repository (http://www.coin-or.org).                          *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2007 Yan Xu and Ted Ralphs.                                 *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef VrpSolution_h_
#define VrpSolution_h_

#include "Alps.h"

#include "BlisSolution.h"

#include "VrpCommonTypes.h"

//#############################################################################
/** This class contains a vrp solution */
//#############################################################################

class VrpSolution : public BlisSolution {

protected:

    best_tours tour;
    
public:
    
    /** Default constructor. */
    VrpSolution() 
	: 
	BlisSolution()
	{}

    /** Useful constructor. */
    VrpSolution(int s, const double *values, double objValue)
	:
	BlisSolution(s, values, objValue)
	{}

    /** Destructor. */
    virtual ~VrpSolution() { }
    
    /** Print the solution.*/
    virtual void print(std::ostream& os) const {
        double nearInt = 0.0;
	for (int j = 0; j < size_; ++j) {
	    if (values_[j] > 1.0e-15 || values_[j] < -1.0e-15) {
                nearInt = floor(values_[j] + 0.5);
                if (ALPS_FABS(nearInt - values_[j]) < 1.0e-6) {
                    os << "x[" << j << "] = " << nearInt << std::endl;
                }
                else {
                    os << "x[" << j << "] = " << values_[j] << std::endl;
                }   
	    }
	}
    }
    
    /** The method that encodes the solution into a encoded object. */
    virtual AlpsEncoded* encode() const {
	AlpsEncoded* encoded = NULL;//new AlpsEncoded("ALPS_SOLUTION");
	//encodeBcps(encoded);
	// Nothing to do for Vrp part.
	return encoded;
    }
  
    /** The method that decodes the solution from a encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded& encoded) const {
	VrpSolution * sol = NULL;//new VrpSolution();
	//sol->decodeBcps(encoded);
	return sol;
    }
    
};

//#############################################################################
//#############################################################################

#endif
