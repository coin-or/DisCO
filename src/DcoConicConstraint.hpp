#ifndef DcoConicConstraint_h_
#define DcoConicConstraint_h_

#include "BcpsObject.h"
#include <OsiConicSolverInterface.hpp>

class OsiRowCut;

//#############################################################################

class DcoConicConstraint : public BcpsConstraint {

protected:
  /** Size of cone */
  int size_;
  /** Cone type, 0 for Lorentz cone, 1 for rotated Lorentz cone */
  OsiLorentzConeType type_;
  /** Variable indices */
  int * members_;
public:
  /** Default constructor. */
  DcoConicConstraint();
  /** Useful constructor. */
  DcoConicConstraint(int size, OsiLorentzConeType type, int const * members);
  /** Destructor. */
  virtual ~DcoConicConstraint();
  /** Copy constructor. */
  DcoConicConstraint(const DcoConicConstraint & rhs);
  /** Return data  */
  /**@{*/
  int getSize() const       { return size_; }
  OsiLorentzConeType getType() const   { return type_; }
  int const * getMembers() const { return members_; }
  /**@}*/
  /** Virtual Functions inherited from BcpsObjcet */
  /**@{*/
  /** Infeasibility of the object
      This is some measure of the infeasibility of the object. It should be
      scaled to be in the range [0.0, 0.5], with 0.0 indicating the object
      is satisfied.

      The preferred branching direction is returned in preferredWay,

      This is used to prepare for strong branching but should also think of
      case when no strong branching

      The object may also compute an estimate of cost of going "up" or
      "down". This will probably be based on pseudo-cost ideas. */
  virtual double infeasibility(BcpsModel *m, int &preferredWay) const;
  /**@}*/

  /** Set data  */
  /**@{*/
  void setData(int size, OsiLorentzConeType type, int const * members) {
    if (size_ < size) {
      if (members_) {
	delete [] members_;
      }
      members_ = new int[size];
    }
    size_ = size;
    std::copy(members, members+size, members_);
    //memcpy(values_, val, sizeof(double) * s);
  }
  /**@}*/

protected:

  /** Pack Dco part into an encoded object. */
  AlpsReturnStatus encodeDco(AlpsEncoded *encoded);

  /** Unpack Dco part from a encode object. */
  AlpsReturnStatus decodeDco(AlpsEncoded &encoded);

public:

  /** Create a OsiRowCut based on this constraint. */
  OsiRowCut *createOsiRowCut();

  /** Compute a hash key. */
  virtual void hashing(BcpsModel *model=NULL);

  /** Check if violates a given lp solution. */
  double violation(double const * lpSolution);

  using AlpsKnowledge::encode ;
  /** Pack into a encode object. */
  virtual AlpsReturnStatus encode(AlpsEncoded *encoded);

  /** Decode a constraint from an encoded object. */
  virtual AlpsKnowledge* decode(AlpsEncoded& encoded) const;
};

//#############################################################################

#endif
