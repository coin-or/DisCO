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
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#ifndef DcoBranchStrategyPseudo_hpp_
#define DcoBranchStrategyPseudo_hpp_

#include <BcpsBranchStrategy.h>

#include <map>

class DcoModel;
class DcoTreeNode;
/*!  This class implements pseudocost branching strategy. In this part we use
  the notation in Achterberg's PhD dissertation.

  Define following,

  \f$f_j ^+ := \lceil \overline{x}_j \rceil - \overline{x}_j \f$<br>
  \f$f_j ^- := \overline{x}_j - \lceil \overline{x}_j \rceil \f$

  Let \f$ \Delta _j ^- \f$ be the objective value change in the down branch,
  and \f$ \Delta _j ^+ \f$ in the upper branch when we branch on variable
  \f$ x_j \f$.

  Let \f$ \varphi _j ^- \f$ be the average of \f$ \frac{\Delta _j ^-}{f_j ^-}
  \f$ considering all the times we branched on \f$ x_j \f$.

  Then we compute score for variable \f$ x_j \f$ as follows,

  \f[ score(\varphi _j ^-, \varphi _j ^+) =
        (1-u) min(\varphi _j ^-, \varphi _j^+)
      + u max(\varphi _j ^-, \varphi _j^+)
  \f]

  where u is called the scale factor and it is fixed to 1/6.

  We branch on the maximum score variable.

  # Design, keeping statistics
  How do we update \f$ \varphi \f$, since optimal value of the children are not
  known before solving the children's subproblem.

  One way to do this is when createCandBranchObjects() is called.

  // what happens in case of cuts?
  // ans: the update of statistics is carried once the node is decided to be
  // branched and the branching process is started. Hence the statistics will
  // have not only the branching affect but the cutting too. This is OK since
  // in the statistics we want to measure how much does the branching results
  // a tighter estimate of the problem. If it leads a tighter estimate with the
  // help of cuts let it be.


 */

class DcoBranchStrategyPseudo: virtual public BcpsBranchStrategy {
  /// score factor used. See class documentation.
  double score_factor_;
  ///@name Statistics
  //@{
  /// number of observations for each integer variable
  int * down_num_;
  int * up_num_;
  /// estimated improvement in the objective value per change in each variable
  /// derivative_[i] is the average of all observations for variable i.
  /// these are \f$ \varphi \f$ variables in the documentation
  double * down_derivative_;
  double * up_derivative_;
  /// reverse map of relaxed columns, rev_relaxed_[index] gives the index of
  /// the varaible in relaxed columns array.
  std::map<int,int> rev_relaxed_;
  /// update scores of the stored branch objects.
  void update_statistics(DcoTreeNode * node);
public:
  DcoBranchStrategyPseudo(DcoModel * model);
  virtual ~DcoBranchStrategyPseudo();
  virtual int createCandBranchObjects(BcpsTreeNode * node);
  /// Compare current to other, return 1 if current is better, 0 otherwise
  virtual int betterBranchObject(BcpsBranchObject const * current,
                                 BcpsBranchObject const * other);
private:
  /// Disable default constructor.
  DcoBranchStrategyPseudo();
  /// Disable copy constructor.
  DcoBranchStrategyPseudo(DcoBranchStrategyPseudo const & other);
  /// Disable copy assignment operator.
  DcoBranchStrategyPseudo & operator=(DcoBranchStrategyPseudo const & rhs);
};

#endif
