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
 * Ted Ralphs. All Rights Reserved.                                          *
 *===========================================================================*/


#ifndef DcoBranchStrategyMaxInf_hpp_
#define DcoBranchStrategyMaxInf_hpp_

#include <BcpsBranchStrategy.h>

class DcoModel;

/**
   # How does branching work?

   When a node is picked from the pool and it is pregnant Alps calls
   AlpsSearchStrategy::createNewNodes(). AlpsSearchStrategy::createNewNodes()
   calls AlpsNode::branch() and Alps processes output of this function further.

   When DcoTreeNode::process() decides to branch it calls
   DcoBranchStrategy[Name]::createCandBranchObjects(). This will populate
   branchObjects_ inherited from BcpsBranchStrategy.

   When DcoTreeNode::branch() is called we create child nodes using the
   bestBranchObject_.  DcoTreeNode::branch() returns the new children nodes.
   Alps will add them into the node pool.

   DcoTreeNode::process calls DcoTreeNode::boundingLoop.
   DcoTreeNode::boundingLoop calls DcoTreeNode::branchConstrainOrPrice.
   DcoTreeNode::branchConstrainOrPrice may decide to branch. When it does
   DcoTreeNode::boundingLoop calls DcoBranchStrategy::createCandBranchObjects.
   Then sets DcoTreeNode::branchObject_ to the best branch object found by
   calling DcoBranchStrategy::bestBranchObject().

   Branching is triggered by an infeasible object. An infeasible object should
   know how to branch itself. Thus, they have member function
   createBranchObject(). This means for each infeasible object we will have a
   branch object.

   A Branch object keeps track of how many branches are created. Has detailed
   information on creating these branches.



   # Ideas/Questions

   Q1: I do not see clearly what a branch object represents. For me each branch is
   a node. I do not see how we need branch object class. Branch itself does not
   stand for something, branching is a process that outputs something (nodes).

   A1: You are right. Branching is a process and it creates new
   nodes. Ideally the branching object data should be stored in the node class.
   In the current implementatition the data needed for branching is separately
   encapsulated in a branch object and it is a member of the DcoTreeNode class.

   Branch strategy object (branchStrategy_) is a member of DcoModel. Once we
   decide to branch when processing a node, we create candidates using the
   strategy object (branchStrategy_->createCandBranchObjects()). But this
   function has no idea what node is being branched. This creates problems
   since branching process might need node data, (ie. pseudocost branching
   needs the quality of the node being processed and its parent's quality).

   # How to update scoring

   Keeping scores for most infeasible branching strategy is easy. Each time
   we will decide the branching variable, we will check the fractionality of
   all variables. Things get mroe complicated for other strategies.

   ## Pseudocost

   When pseudocost strategy is used, we need to update statistics related to
   performance of variables. When a branch decision is being made data only
   from previous decisions are available. The performance of current branching
   decison will be available only after the children's subproblems are solved.

   When the statistics should be updated???? Ideally it should be updated
   after we solve a subproblem, DcoTreeNode::bound() call.

   When a branch() is called in a node, first update its score. This way we
   update scores for each node only once.

   # BcpsBranchStrategy

   BcpsBranchStrategy is an abstract base class for defining a branch stragey.
   It provides an interface for creating a set of branching candidates and
   comparing them.

   # Maximum infeasibility branching, DcoBranchStrategyMaxInf

   This class represents maximum infeasibility branching strategy. Inherits
   BcpsBranchStrategy.

   Implements maximum infeasibility branching in interface inherited from
   BcpsBranchStrategy.

   # What is a branch object?

   A branching object keeps necessary information for a branching to be
   perfomed.

   BcpsBranchObject contains the member data required when choosing
   branching entities and excuting actual branching. It also has
   the member funtions to do branching by adjusting bounds, etc.
   in solver. Branching objects can be simple integer variables or more
   complicated objects like SOS.

   BcpsBranchObject fields are type_, model_, objectIndex_, upScore_,
   downScore_, direction_, value_, numBranchesLeft_.

   virtual functions are: clone(), numBranches(), numBranchesLeft(), branch(),
   print(), boundBranch(),

 */
class DcoBranchStrategyMaxInf: virtual public BcpsBranchStrategy {
public:
  DcoBranchStrategyMaxInf(DcoModel * model);
  virtual ~DcoBranchStrategyMaxInf() {}
  virtual int createCandBranchObjects(BcpsTreeNode * node);
  /// Compare current to other, return 1 if current is better, 0 otherwise
  virtual int betterBranchObject(BcpsBranchObject const * current,
                                 BcpsBranchObject const * other);
private:
  DcoBranchStrategyMaxInf();
  DcoBranchStrategyMaxInf & operator=(DcoBranchStrategyMaxInf const & rhs);
  DcoBranchStrategyMaxInf(DcoBranchStrategyMaxInf const & other);
};

#endif
