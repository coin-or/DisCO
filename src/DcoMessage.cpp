/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           *
 *                                                                           *
 * Copyright (C) 2001-2015, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "DcoMessage.hpp"
#include <cstring>

//#############################################################################
typedef struct {
    DISCO_Message internalNumber;
    int externalNumber;              // or continuation
    char detail;
    const char * message;
} Dco_message;

//#############################################################################

//
// 100-200 DcoModel::readInstance messages
// 200-300 DcoTreeNode messages
//
//
static Dco_message us_english[]=
{
    {DISCO_CUTOFF_INC, 43, 1, "Objective coefficients are multiples of %g"},
    {DISCO_CUT_STAT_FINAL, 53, 1, "Called %s cut generator %d times, generated %d cuts, CPU time %.4f seconds, current strategy %d"},
    {DISCO_CUT_STAT_NODE, 55, 1, "Node %d, called %s cut generator %d times, generated %d cuts, CPU time %.4f seconds, current strategy %d"},
    {DISCO_GAP_NO, 57, 1, "Relative optimality gap is infinity because no solution was found"},
    {DISCO_GAP_YES, 58, 1, "Relative optimality gap is %.2f%%"},
    {DISCO_HEUR_BEFORE_ROOT, 60, 1, "%s heuristic found a solution; quality is %g"},
    {DISCO_HEUR_STAT_FINAL, 63, 1, "Called %s heuristic %d times, found %d solutions, CPU time %.4f seconds, current strategy %d"},
    {DISCO_HEUR_STAT_NODE, 65, 1, "Node %d, called %s heuristic %d times, found %d solutions, CPU time %.4f seconds, current strategy %d"},
    {DISCO_ROOT_PROCESS, 30, 1, "Processing the root node (%d rows, %d columns)"},
    {DISCO_ROOT_TIME, 35, 1, "Processing the first root relaxation took CPU time %.4f seconds"},
    // reading mps files
    {DISCO_READ_NOINTS, 20, 1, "Problem does not have integer variables"},
    {DISCO_READ_NOCONES, 21, 1, "Problem does not have conic constraints."},
    {DISCO_READ_MPSERROR, 9001, 1, "Reading conic mps file failed with code %d." },
    {DISCO_READ_MPSFILEONLY,9002, 1, "Mosek style conic mps files only."},
    {DISCO_READ_CONEERROR, 9002, 1, "Invalid cone type."},
    {DISCO_READ_ROTATEDCONESIZE, 9002, 1, "Rotated cones should have at least 3 members."},
    {DISCO_READ_CONESTATS1,101,1,"Problem has %d cones."},
    {DISCO_READ_CONESTATS2,102,1, "Cone %d has %d entries (type %d)"},
    {DISCO_NODE_BRANCHONINT,9201,1, "Branched on integer variable. Variable index %d."},
    {DISCO_NODE_UNEXPECTEDSTATUS,9202,1, "Unexpected node status %d"},
    // relaxation solver messages
    {DISCO_SOLVER_UNKNOWN_STATUS,9301,1, "Unknown relaxation solver status."},
    // general messages
    {DISCO_OUT_OF_MEMORY,9901,1, "Out of memory, file: %s, line: %d."},
    {DISCO_NOT_IMPLEMENTED,9902,1, "Not implemented yet, file: %s, line: %d."},
    {DISCO_DUMMY_END, 9999, 0, ""}
};

//#############################################################################

/* Constructor */
DcoMessage::DcoMessage(Language language):
  CoinMessages(sizeof(us_english) / sizeof(Dco_message)) {
  language_ = language;
  strcpy(source_, "Dco");
  Dco_message * message = us_english;
  while (message->internalNumber != DISCO_DUMMY_END) {
    CoinOneMessage oneMessage(message->externalNumber, message->detail,
			      message->message);
    addMessage(message->internalNumber, oneMessage);
    message++;
  }
  // now override any language ones
  switch (language) {
  default:
    message = NULL;
    break;
  }
  // replace if any found
  if (message) {
    while (message->internalNumber != DISCO_DUMMY_END) {
      replaceMessage(message->internalNumber, message->message);
      message++;
    }
  }
}
