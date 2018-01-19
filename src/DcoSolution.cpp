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


#include "DcoSolution.hpp"
#include "DcoMessage.hpp"
#include "DcoModel.hpp"

DcoSolution::DcoSolution() {
}

DcoSolution::DcoSolution(int size, double const * values, double quality):
  BcpsSolution(size, values, quality) {
}

DcoSolution::~DcoSolution() {
}

BcpsSolution * DcoSolution::selectNonzeros(const double etol) const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}

BcpsSolution * DcoSolution::selectFractional(const double etol) const {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
  return NULL;
}

/// Encodes the solution into AlpsEncoded object and return pointer to it.
AlpsReturnStatus DcoSolution::encode(AlpsEncoded * encoded) const {
  // get pointers for message logging
  assert(broker_);
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  //CoinMessages * messages = model->dcoMessages_;

  // return value
  AlpsReturnStatus status;
  status = AlpsSolution::encode(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected encode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  // encode Bcps fields
  status = BcpsSolution::encode(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected encode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }

  // Nothing to do for DisCO part.

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "Proc[" << broker_->getProcRank() << "]"
            << " solution " << this << " encoded, quality "
            << quality_ << std::endl;
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                              'G', DISCO_DLOG_MPI)
    << CoinMessageEol;
  // end of debug stuff


  return status;
}

/// Decodes the given object into a new solution and returns the pointer to
/// it.
AlpsKnowledge * DcoSolution::decode(AlpsEncoded & encoded) const {
  AlpsReturnStatus status;
  DcoSolution * sol = new DcoSolution();
  sol->setBroker(broker_);
  status = sol->decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected decode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  return sol;
}

/// Decode the given AlpsEncoded object into this.
AlpsReturnStatus DcoSolution::decodeToSelf(AlpsEncoded & encoded) {
  // get pointers for message logging
  assert(broker_);
  DcoModel * model = dynamic_cast<DcoModel*>(broker_->getModel());
  CoinMessageHandler * message_handler = model->dcoMessageHandler_;
  //CoinMessages * messages = model->dcoMessages_;

  AlpsReturnStatus status;
  status = AlpsSolution::decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected decode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }
  status = BcpsSolution::decodeToSelf(encoded);
  if (status!=AlpsReturnStatusOk) {
    std::cerr << "Unexpected decode status, "
              << "file: " <<  __FILE__
              << "line: " << __LINE__
              << std::endl;
    throw std::exception();
  }

  // debug stuff
  std::stringstream debug_msg;
  debug_msg << "Proc[" << broker_->getProcRank() << "]"
            << " solution decoded into " << this << ". quality "
            << quality_ << std::endl;
  message_handler->message(0, "Dco", debug_msg.str().c_str(),
                              'G', DISCO_DLOG_MPI)
    << CoinMessageEol;
  // end of debug stuff

  return status;
}
