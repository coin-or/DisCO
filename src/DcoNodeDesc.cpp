/*===========================================================================*
 * This file is part of the Discrete Conic Optimization (DisCO) Solver.      *
 *                                                                           *
 * DisCO is distributed under the Eclipse Public License as part of the      *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *          Aykut Bulut, Lehigh University                                   *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2001-2018, Lehigh University, Aykut Bulut, Yan Xu, and      *
 *                          Ted Ralphs.                                      *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#include "DcoNodeDesc.hpp"
#include "DcoMessage.hpp"

DcoNodeDesc::DcoNodeDesc() {
  // set if as down branch by default
  branchedDir_ = DcoNodeBranchDirectionDown;
  branchedInd_ = -1;
  // set to 0.0 by default.
  branchedVal_ = 0.0;
  basis_ = NULL;
}

DcoNodeDesc::DcoNodeDesc(DcoModel * model): BcpsNodeDesc(model) {
  // set if as down branch by default
  branchedDir_ = DcoNodeBranchDirectionDown;
  branchedInd_ = -1;
  // set to 0.0 by default.
  branchedVal_ = 0.0;
  basis_ = NULL;
}

DcoNodeDesc::~DcoNodeDesc() {
  if (basis_) {
    delete basis_;
  }
}

DcoNodeBranchDir DcoNodeDesc::getBranchedDir() const {
  return branchedDir_;
}

int DcoNodeDesc::getBranchedInd() const {
  return branchedInd_;
}

double DcoNodeDesc::getBranchedVal() const {
  return branchedVal_;
}

// set fields
void DcoNodeDesc::setBranchedDir(DcoNodeBranchDir dir) {
  branchedDir_ = dir;
}

void DcoNodeDesc::setBranchedInd(int ind) {
  branchedInd_ = ind;
}

void DcoNodeDesc::setBranchedVal(double val) {
  branchedVal_ = val;
}

void DcoNodeDesc::setBasis(CoinWarmStartBasis *& ws) {
  if (basis_) {
    delete basis_;
  }
  basis_= ws;
  ws = NULL;
}

/** Get warm start basis. */
CoinWarmStartBasis * DcoNodeDesc::getBasis() const {
  return basis_;
}

/// Encode this to an AlpsEncoded object.
AlpsReturnStatus DcoNodeDesc::encode(AlpsEncoded * encoded) const {
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  // return value
  AlpsReturnStatus status;
  status = AlpsNodeDesc::encode(encoded);
  // todo(aykut) rename this function in Bcps level?
  status = BcpsNodeDesc::encodeBcps(encoded);
  assert(status==AlpsReturnStatusOk);
  if (status != AlpsReturnStatusOk) {
    message_handler->message(DISCO_UNEXPECTED_ENCODE_STATUS,
                             *messages)
      << __FILE__ << __LINE__ << CoinMessageEol;
  }
  encoded->writeRep(branchedDir_);
  encoded->writeRep(branchedInd_);
  encoded->writeRep(branchedVal_);
  // Encode basis if available
  int available = 0;
  if (basis_) {
    available = 1;
    encoded->writeRep(available);
    int numCols = basis_->getNumStructural();
    int numRows = basis_->getNumArtificial();
    encoded->writeRep(numCols);
    encoded->writeRep(numRows);
    // Pack structural.
    int nint = (basis_->getNumStructural() + 15) >> 4;
    encoded->writeRep(basis_->getStructuralStatus(), nint * 4);
    // Pack artificial.
    nint = (basis_->getNumArtificial() + 15) >> 4;
    encoded->writeRep(basis_->getArtificialStatus(), nint * 4);
  }
  else {
    encoded->writeRep(available);
  }
  return status;
}

/// Decode a given AlpsEncoded object to an AlpsKnowledge object and return a
/// pointer to it.
AlpsNodeDesc * DcoNodeDesc::decode(AlpsEncoded & encoded) const {
  // get pointers for message logging
  DcoModel * dco_model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = dco_model->dcoMessageHandler_;
  CoinMessages * messages = dco_model->dcoMessages_;
  AlpsReturnStatus status;
  DcoNodeDesc * new_desc = new DcoNodeDesc();
  status = new_desc->decodeToSelf(encoded);
  if (status != AlpsReturnStatusOk) {
    message_handler->message(DISCO_UNEXPECTED_DECODE_STATUS,
                             *messages)
      << __FILE__ << __LINE__ << CoinMessageEol;
  }
  return new_desc;
}

/// Decode a given AlpsEncoded object into self.
AlpsReturnStatus DcoNodeDesc::decodeToSelf(AlpsEncoded & encoded) {
  AlpsReturnStatus status = AlpsReturnStatusOk;
  status = AlpsNodeDesc::decodeToSelf(encoded);
  // todo(aykut) rename this function in Bcps level?
  status = BcpsNodeDesc::decodeBcps(encoded);

  assert(status==AlpsReturnStatusOk);
  encoded.readRep(branchedDir_);
  encoded.readRep(branchedInd_);
  encoded.readRep(branchedVal_);
  // decode basis if available
  int available;
  encoded.readRep(available);
  if (available==1) {
    if (basis_) {
      delete basis_;
    }
    int numCols;
    int numRows;
    encoded.readRep(numCols);
    encoded.readRep(numRows);
    int tempInt;
    // Structural
    int nint = (numCols + 15) >> 4;
    char * structuralStatus = new char[4 * nint];
    encoded.readRep(structuralStatus, tempInt);
    assert(tempInt == nint*4);
    // Artificial
    nint = (numRows + 15) >> 4;
    char * artificialStatus = new char[4 * nint];
    encoded.readRep(artificialStatus, tempInt);
    assert(tempInt == nint*4);
    basis_ = new CoinWarmStartBasis();
    if (!basis_) {
      // todo(aykut) We should do this through messages as in the rest of the code.
        throw CoinError("Out of memory", "DcoDecodeWarmStart", "HELP");
    }
    basis_->assignBasisStatus(numCols, numRows,
                              structuralStatus, artificialStatus);
    assert(!structuralStatus);
    assert(!artificialStatus);
  }
  else {
    basis_ = NULL;
  }
  return status;
}
