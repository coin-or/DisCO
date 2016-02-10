#ifndef DcoNodeDesc_hpp_
#define DcoNodeDesc_hpp_

#include <CoinWarmStartBasis.hpp>
#include <CoinWarmStartBasis.hpp>
#include <BcpsNodeDesc.h>
#include "DcoModel.hpp"
#include "Dco.hpp"

class DcoNodeDesc: public BcpsNodeDesc {
  /** Branched direction to create it. For updating pseudocost. */
  DcoNodeBranchDir branchedDir_;
  /** Branched object index to create it. For updating pseudocost. */
  int branchedInd_;
  /** Branched value to create it. For updating pseudocost. */
  double branchedVal_;
  /** Warm start. */
  CoinWarmStartBasis * basis_;
public:
  DcoNodeDesc();
  DcoNodeDesc(DcoModel * model);
  virtual ~DcoNodeDesc();
  DcoNodeBranchDir branchedDir() const;
  int branchedInd() const;
  double branchedVal() const;
  // set fields
  void setBranchedDir(DcoNodeBranchDir dir);
  void setBranchedInd(int ind);
  void setBranchedVal(double val);
  /** Set basis. */
  void setBasis(CoinWarmStartBasis *& ws);
  /** Get warm start basis. */
  CoinWarmStartBasis * getBasis() const;
};

#endif
