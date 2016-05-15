#ifndef DcoNodeDesc_hpp_
#define DcoNodeDesc_hpp_

#include <CoinWarmStartBasis.hpp>
#include <CoinWarmStartBasis.hpp>
#include <BcpsNodeDesc.h>
#include "DcoModel.hpp"
#include "Dco.hpp"

/*!
  Represents subproblem data, ie. tree node data. This class stores a branch
  and bound tree node data. Inherits BcpsNodeDesc. BcpsNodeDesc inherits
  AlpsNodeDesc.

  # AlpsNodeDesc

  AlpsNodeDesc is an abstract base class for Alps application node data. It
  does not have any pure virtual functions but it will have (this will get
  fixed.). It has one field only and t is model_. It only provides interface
  for encode and decode functions.

  # BcpsNodeDesc

  Inherits AlpsNodeDesc. It has two fields, vars_ and cons_. Both of these
  fields are pointers to BcpsObjectListMod instances. It provides an interface
  to modify these two fields.

  # BcpsObjectListMod

  This is a simple struct with data members only. It has the following members.

  <ul>
  <li> numRemove, integer. Number of objects to be deleted.
  <li> posRemove, integer pointer.
  <li> numAdd, integer. Number of objects to be added.
  <li> objects, BcpsObject pointer pointer.
  <li> lbHard, BcpsFieldListMod<double>.
  <li> ubHard, BcpsFieldListMod<double>.
  <li> lbSoft, BcpsFieldListMod<double>.
  <li> ubSoft, BcpsFieldListMod<double>.
  <li> status, BcpsFieldListMod<int>.
  </ul>

  Struct BcpsFieldListMod explained next.

  # BcpsFieldListMod

  This is a simple struct with a template class. It has the following fields.

  <ul>
  <li> relative, boolean.
  <li> numModify, integer.
  <li> posModify, integer pointer.
  <li> entries, template class T pointer.
  </ul>

 */

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
  DcoNodeBranchDir getBranchedDir() const;
  int getBranchedInd() const;
  double getBranchedVal() const;
  // set fields
  void setBranchedDir(DcoNodeBranchDir dir);
  void setBranchedInd(int ind);
  void setBranchedVal(double val);
  /** Set basis. */
  void setBasis(CoinWarmStartBasis *& ws);
  /** Get warm start basis. */
  CoinWarmStartBasis * getBasis() const;
  /// get DcoModel this node belongs.
  DcoModel * getModel() const;
};

#endif
