#ifndef DcoNodeDesc_hpp_
#define DcoNodeDesc_hpp_

#include <BcpsNodeDesc.h>
#include "DcoModel.hpp"


class DcoNodeDesc: public BcpsNodeDesc {
public:
  DcoNodeDesc();
  DcoNodeDesc(DcoModel * model);
  virtual ~DcoNodeDesc();
};

#endif
