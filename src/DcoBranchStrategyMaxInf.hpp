#ifndef DcoBranchStrategyMaxInf_hpp_
#define DcoBranchStrategyMaxInf_hpp_

#include <BcpsBranchStrategy.h>

class DcoModel;

/**
   Represents maximum infeasibility branching strategy. Inherits
   BcpsBranchStrategy.

   BcpsBranchStrategy:
   BcpsBranchStrategy is an abstract base class for defining a branch stragey.
   It provides an interface for creating a set of branching candidates and
   comparing them.

   DcoBranchStrategyMaxInf:
   Implements maximum infeasibility branching in interface inherited from
   BcpsBranchStrategy.

 */
class DcoBranchStrategyMaxInf: public BcpsBranchStrategy {
public:
  DcoBranchStrategyMaxInf();
  DcoBranchStrategyMaxInf(DcoBranchStrategyMaxInf const & other);
  DcoBranchStrategyMaxInf(DcoModel * model);
  virtual ~DcoBranchStrategyMaxInf();
  virtual BcpsBranchStrategy * clone() const;
  virtual int createCandBranchObjects(int numPassesLeft, double ub);
  virtual int betterBranchObject(BcpsBranchObject * b,
				 BcpsBranchObject * bestSoFar);
private:
  DcoBranchStrategyMaxInf & operator=(DcoBranchStrategyMaxInf const & rhs);
};

#endif
