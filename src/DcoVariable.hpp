#ifndef DcoVariable_hpp_
#define DcoVariable_hpp_

#include <BcpsObject.h>

class DcoVariable: public BcpsVariable {
public:
  DcoVariable();
  DcoVariable(double lbh, double ubh, double lbs, double ubs);
  DcoVariable(DcoVariable const & other);
  virtual ~DcoVariable();
};

#endif
