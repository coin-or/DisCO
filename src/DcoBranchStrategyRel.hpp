#ifndef DcoBranchStrategyRel_hpp_
#define DcoBranchStrategyRel_hpp_

#include <BcpsBranchStrategy.h>

class DcoModel;

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
