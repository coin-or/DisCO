#include <algorithm>
#include <cmath>

#include "DcoConicConstraint.hpp"
#include "DcoModel.hpp"

class OsiRowCut;

#define EPS1 1e-5

/** Default constructor. */
DcoConicConstraint::DcoConicConstraint() {
  members_ = 0;
  size_ = 0;
  // cone is lorentz cone by default.
  type_ = OSI_QUAD;
}

/** Useful constructor. */
DcoConicConstraint::DcoConicConstraint(int size, OsiLorentzConeType type, int const * members) {
  size_ = size;
  type_ = type;
  members_ = new int[size];
  std::copy(members, members+size, members_);
}

/** Destructor. */
DcoConicConstraint::~DcoConicConstraint() {
  if (members_) {
    delete[] members_;
  }
}

/** Copy constructor. */
DcoConicConstraint::DcoConicConstraint(const DcoConicConstraint & rhs) {
  setData(rhs.getSize(), rhs.getType(), rhs.getMembers());
}

double DcoConicConstraint::infeasibility(BcpsModel * m,
					 int & preferredWay) const {
  DcoModel * model = dynamic_cast<DcoModel*>(m);
  double const * sol = model->getLpSolution();
  double * par_sol = new double[size_];
  for (int i=0; i<size_; ++i) {
    par_sol[i] = sol[members_[i]];
  }
  double infeas;
  if (type_==OSI_QUAD) {
    double term1 = par_sol[0];
    double term2 = std::inner_product(par_sol+1, par_sol+size_,
				      par_sol+1, 0.0);
    term2 = sqrt(term2);
    infeas = term2-term1;
    if (infeas<EPS1) {
      infeas = 0.0;
    }
  }
  else if (type_==OSI_RQUAD) {
    double term1 = 2.0*par_sol[0]*par_sol[1];
    double term2 = std::inner_product(par_sol+2, par_sol+size_,
				      par_sol+2, 0.0);
    infeas = term2-term1;
    if (infeas<EPS1) {
      infeas = 0.0;
    }
  }
  else {
    std::cerr << "Unknown cone type!" << std::endl;
    throw std::exception();
  }
  delete[] par_sol;
  return infeas;
}

/** Pack Dco part into an encoded object. */
AlpsReturnStatus DcoConicConstraint::encodeDco(AlpsEncoded *encoded) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

/** Unpack Dco part from a encode object. */
AlpsReturnStatus DcoConicConstraint::decodeDco(AlpsEncoded &encoded) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

/** Create a OsiRowCut based on this constraint. */
OsiRowCut * DcoConicConstraint::createOsiRowCut() {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

/** Compute a hash key. */
void DcoConicConstraint::hashing(BcpsModel * model) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

/** Check if violates a given lp solution. */
double DcoConicConstraint::violation(double const * lpSolution) {
  double viol = 0.0;
  double * par_sol = new double[size_];
  for (int j=0; j<size_; ++j) {
    par_sol[j] = lpSolution[members_[j]];
  }
  if (type_==OSI_QUAD) {
    double term1 = par_sol[0];
    double term2 = std::inner_product(par_sol+1, par_sol+size_,
				      par_sol+1, 0.0);
    term2 = sqrt(term2);
    viol = term2 - term1;
  }
  else if (type_==OSI_RQUAD) {
    double term1 = 2.0*par_sol[0]*par_sol[1];
    double term2 = std::inner_product(par_sol+2, par_sol+size_,
				      par_sol+2, 0.0);
    viol = term2 - term1;
  }
  else {
    std::cerr << "Unknown cone type!" << std::endl;
    throw std::exception();
  }
  delete[] par_sol;
  // todo(aykut) this epsilon should be a parameter
  if (viol<EPS1) {
    return 0.0;
  }
  return viol;
}

/** Pack into a encode object. */
AlpsReturnStatus DcoConicConstraint::encode(AlpsEncoded *encoded) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

/** Decode a constraint from an encoded object. */
AlpsKnowledge * DcoConicConstraint::decode(AlpsEncoded& encoded) const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}
