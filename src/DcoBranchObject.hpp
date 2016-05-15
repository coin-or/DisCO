#ifndef DcoBranchObject_hpp_
#define DcoBranchObject_hpp_

#include <BcpsBranchObject.h>
#include "DcoModel.hpp"

/**
   Represents a DisCO branch object. DcoBranchObject inherits BcpsBranchObject.

   BcpsBranchObject:
   BcpsBranchObject is an abstract base class for describing a branch object. A
   branch object keeps data to create a branch. Data stored are type_, model_,
   objectIndex_, upScore_, downScore_, direction_, value and numBranchesLeft_.
   It has an interface that gets/sets these data. See BcpsBranchObject
   documentation for details.

   BcpsBranchObject stands for the most general branch objects, SOS, variable,
   etc. DcoBranchObject is for a simple branch on an integral
   variable. Implements BcpsBranchObject this in mind.
 */

class DcoBranchObject: virtual public BcpsBranchObject {
  /// upper bound of the down branch
  double ubDownBranch_;
  /// lower bound of the up branch
  double lbUpBranch_;
public:
  ///@name Constructor and Destructors.
  //@{
  /// Constructor.
  DcoBranchObject(DcoModel * model, int colInd, int intScore, double dblScore,
                  int direction, double value);
  /// Destructor.
  virtual ~DcoBranchObject();
  //@}

  ///@name Other functions
  //@{
  /// Clone this to a new object and return pointer to it.
  virtual BcpsBranchObject * clone() const;
  //@}

  ///@name Branch related functions.
  //@{
  /// Perform branching as specified by the branching object. Update the status
  /// of this branching object.
  virtual double branch(bool normalBranch = false);
  /// Return true if branching should fix object bounds.
  virtual bool boundBranch() const;
  //@}

  ///@name Bound getting functions.
  //@{
  /// Get upper bound of the down branch.
  double ubDownBranch() const { return ubDownBranch_; }
  /// Get lower bound of the up branch.
  double lbUpBranch() const { return lbUpBranch_; }
  //@}

};

#endif
