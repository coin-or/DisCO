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

#include <vector>
#include "VrpHeurTSP.h"
#include "VrpModel.h"
#include "VrpVariable.h"

//#############################################################################

// Given current LP solution. Rouding number closed to 1 to 1
// Let those edge be in the route, then use nearest neighbor idea
// to construct a tour.
bool
VrpHeurTSP::searchSolution(double & objectiveValue, double * newSolution)
{
    bool feasible = false;

    const double * sol = model_->getLpSolution();
    

    
    return feasible;
}

//#############################################################################

void
VrpHeurTSP::createAdjList(VrpModel *model) 
{
    std::vector<VrpVariable *> edges = model->getEdgeList();
    int numVertices = model->getNumVertices();
    int numEdges = edges.size();
    
}

//#############################################################################
#if 0
void VrpHeurTSP::findNearest(VrpModel * model)
{
    std::vector<VrpVariable *> edges = model->getEdgeList();
    
    VrpVariable *var = NULL;
    int beg, end, k = 0;
    double cost;

    int numVertices = model->getNumVertices();
    int numEdges = edges.size();
    
    assert(numEdges == model->getNumEdges());
    
    // Allocate memory for nearest_ and weights_;
    nearest_ = new int [numVertices];
    weights_ = new double [numVertices]; 
    CoinZeroN(nearest_, numVertices);
    CoinFillN(weights_, numVertices, 1.0e20);
    for (k = 0; k < numEdges; ++k) {
        cost = var->getObjCoef();
        beg = var->getv0();
        end = var->getv1();
        if (cost < weights_[beg]) {
            weights_[beg] = cost;
        }
        if (cost < weights_[end]) {
            weights_[end] = cost;
        }
    }

    for (k = 0; k < numVertices; ++k) {
        assert(weights_[k] < 1.0e19);
    }
}
#endif

//#############################################################################
