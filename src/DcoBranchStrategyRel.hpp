#ifndef DcoBranchStrategyRel_hpp_
#define DcoBranchStrategyRel_hpp_

#include <BcpsBranchStrategy.h>

class DcoModel;

/**
   Represents Reliability branching strategy. Inherits BcpsBranchStrategy.

   BcpsBranchStrategy:
   BcpsBranchStrategy is an abstract base class for defining a branch stragey.
   It provides an interface for creating a set of branching candidates and
   comparing them.

   DcoBranchStrategyRel:
   Implements reliability branching in interface inherited from
   BcpsBranchStrategy.

 */
class DcoBranchStrategyRel: public BcpsBranchStrategy {
public:
  DcoBranchStrategyRel();
  DcoBranchStrategyRel(DcoBranchStrategyRel const & other);
  DcoBranchStrategyRel(DcoModel * model);
  virtual ~DcoBranchStrategyRel();
  virtual BcpsBranchStrategy * clone() const;
  virtual int createCandBranchObjects(int numPassesLeft, double ub);
  virtual int betterBranchObject(BcpsBranchObject * b,
				 BcpsBranchObject * bestSoFar);
private:
  DcoBranchStrategyRel & operator=(DcoBranchStrategyRel const & rhs);
};

#endif
