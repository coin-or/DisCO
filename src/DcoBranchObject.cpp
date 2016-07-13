#include "DcoBranchObject.hpp"
#include "DcoMessage.hpp"

DcoBranchObject::DcoBranchObject(int index, double score, double value)
  : BcpsBranchObject(DcoBranchingObjectTypeInt, index, score, value) {
  ubDownBranch_ = floor(value);
  lbUpBranch_ = ceil(value);
}

/// Copy constructor.
DcoBranchObject::DcoBranchObject(DcoBranchObject const & other)
  : BcpsBranchObject(other) {
  ubDownBranch_ = other.ubDownBranch();
  lbUpBranch_ = other.lbUpBranch();
}

/// Helpful Copy constructor.
DcoBranchObject::DcoBranchObject(BcpsBranchObject const * other)
  : AlpsKnowledge(other->getType(), other->broker_),
    BcpsBranchObject(*other) {
  DcoBranchObject const * dco_other = dynamic_cast<DcoBranchObject const *>(other);
  if (dco_other==NULL) {
    std::cerr << "Fatal error!" << std::endl;
    throw std::exception();
  }
  ubDownBranch_ = dco_other->ubDownBranch();
  lbUpBranch_ = dco_other->lbUpBranch();
}

/// Copy assignment operator
DcoBranchObject & DcoBranchObject::operator=(DcoBranchObject const & rhs) {
  BcpsBranchObject::operator=(rhs);
  ubDownBranch_ = rhs.ubDownBranch();
  lbUpBranch_ = rhs.lbUpBranch();
  return *this;
}

DcoBranchObject::~DcoBranchObject() {
}

// BcpsBranchObject * DcoBranchObject::clone() const {
//   DcoModel * dco_model = dynamic_cast<DcoModel*> (model());
//   int index = getObjectIndex();
//   double up_score = getUpScore();
//   double down_score = getDownScore();
//   int direction = getDirection();
//   double value = getValue();
//   BcpsBranchObject * nbo = new DcoBranchObject(dco_model, index, up_score,
//                                                down_score, direction, value);
//   return nbo;
// }

// double DcoBranchObject::branch(bool normalBranch) {
//   DcoModel * model = dynamic_cast<DcoModel*> (model_);
//   int col_index = model->relaxedCols()[objectIndex_];
//   // Decrement number of branches left by 1.
//   --numBranchesLeft_;
//   if (direction_<0) {
//     model->solver()->setColUpper(col_index, ubDownBranch_);
//     direction_ = 1;
//   }
//   else {
//     model->solver()->setColLower(col_index, lbUpBranch_);
//     direction_ = -1;          // Swap direction
//   }
//   return 0.0;
// }

// bool DcoBranchObject::boundBranch() const {
//   // return True if branching fixes the bound.
//   // return True if variable we are branching is binary
//   return false;
// }


/// The number of branch arms created for this branch object.
int DcoBranchObject::numBranches() const {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
}

/// The number of branch arms left to be evaluated.
int DcoBranchObject::numBranchesLeft() const {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
}

/// Spit out a branch and, update this or superclass fields if necessary.
double DcoBranchObject::branch(bool normalBranch) {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  message_handler->message(DISCO_NOT_IMPLEMENTED, *messages)
    << __FILE__ << __LINE__ << CoinMessageEol;
}

/// Encode the content of this into the given AlpsEncoded object.
AlpsReturnStatus DcoBranchObject::encode(AlpsEncoded * encoded) const {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  AlpsReturnStatus status;
  status = BcpsBranchObject::encode(encoded);
  assert(status==AlpsReturnStatusOk);
  encoded->writeRep(ubDownBranch_);
  encoded->writeRep(lbUpBranch_);
  return status;
}

/// Decode the given AlpsEncoded object into a new AlpsKnowledge object and
/// return a pointer to it.
AlpsKnowledge * DcoBranchObject::decode(AlpsEncoded & encoded) const {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  AlpsReturnStatus status;
  // create a new object with default values,
  // Bcps decode will decode right values into them.
  AlpsKnowledge * new_bo = new DcoBranchObject(-1, 0.0, 0.0);
  status = new_bo->decodeToSelf(encoded);
  assert(status=AlpsReturnStatusOk);
  return new_bo;;
}

/// Decode the given AlpsEncoded object into this.
AlpsReturnStatus DcoBranchObject::decodeToSelf(AlpsEncoded & encoded) {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  AlpsReturnStatus status;
  // decode Bcps part.
  status = BcpsBranchObject::decodeToSelf(encoded);
  assert(status=AlpsReturnStatusOk);
  // decode fields of DcoBranchObject
  encoded.readRep(ubDownBranch_);
  encoded.readRep(lbUpBranch_);
  return status;
}
