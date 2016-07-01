#ifndef DcoSolution_hpp_
#define DcoSolution_hpp_

#include <BcpsSolution.h>

/*!

  DcoSolution represents a solution for the MISOCO problem (master problem)
  represented by the DcoModel. Inherits BcpsSolution.

  # Design concerns:

  BcpsSolution can also keep objects (variables) associated with the solution.
  We are not using that for now, but it can be usefull.

 */

class DcoSolution: public BcpsSolution {
public:
  DcoSolution();
  DcoSolution(int size, double const * values, double quality);
  virtual ~DcoSolution();
  virtual BcpsSolution * selectNonzeros(const double etol=1e-5) const;
  virtual BcpsSolution * selectFractional(const double etol=1e-5) const;

  ///@name Encode and Decode functions
  //@{
  /// Get encode defined in AlpsKnowledge.
  /// Encodes the solution into AlpsEncoded object and return pointer to it.
  using AlpsKnowledge::encode;
  /// Encode this into the given AlpsEncoded object.
  virtual AlpsReturnStatus encode(AlpsEncoded * encoded) const;
  /// Decodes the given object into a new solution and returns the pointer to
  /// it.
  virtual AlpsKnowledge * decode(AlpsEncoded & encoded) const;
  /// Decode the given AlpsEncoded object into this.
  virtual AlpsReturnStatus decodeToSelf(AlpsEncoded & encoded);
  //@}
private:
  DcoSolution(DcoSolution const & other);
  DcoSolution & operator=(DcoSolution const & rhs);
};

#endif
