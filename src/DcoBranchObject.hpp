#ifndef DcoBranchObject_hpp_
#define DcoBranchObject_hpp_

#include <BcpsBranchObject.h>
#include "DcoModel.hpp"

/*!
  Represents a DisCO branch object. DcoBranchObject inherits BcpsBranchObject.

  # DcoBranchObject

  DcoBranchObject represents a branch object for a simple branch on an
  integral variable.

  Has two fields, ubDownBranch_ and lbUpBranch_.

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
  DcoBranchObject(DcoModel * model, int index, double score, double value);
  /// Copy constructor.
  DcoBranchObject(DcoBranchObject const & other);
  /// Helpful Copy constructor.
  DcoBranchObject(BcpsBranchObject const * other);
  /// Copy assignment operator
  DcoBranchObject & operator=(DcoBranchObject const & rhs);
  /// Destructor.
  virtual ~DcoBranchObject();
  //@}

  ///@name Other functions
  //@{
  /// Clone this to a new object and return pointer to it.
  //virtual BcpsBranchObject * clone() const;
  //@}

  ///@name Virtual functions inherited from BcpsBranchObject
  /// The number of branch arms created for this branch object.
  virtual int numBranches() const;
  /// The number of branch arms left to be evaluated.
  virtual int numBranchesLeft() const;
  /// Spit out a branch and, update this or superclass fields if necessary.
  virtual double branch(bool normalBranch = false);
  /// Pack to an encoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// Unpack a branching object from an encoded object.
  virtual AlpsReturnStatus decode(AlpsEncoded & encoded);
  //@}

  ///@name Bound getting functions.
  //@{
  /// Get upper bound of the down branch.
  double ubDownBranch() const { return ubDownBranch_; }
  /// Get lower bound of the up branch.
  double lbUpBranch() const { return lbUpBranch_; }
  //@}

private:
  /// Disable default constructor.
  DcoBranchObject();
};

#endif
